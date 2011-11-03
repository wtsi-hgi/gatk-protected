/*
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.walkers.genotyper.MultiallelicGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.indels.ConstrainedMateFixingManager;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

public class GenotypingEngine {

    // Smith-Waterman parameters originally copied from IndelRealigner
    private final double SW_MATCH = 5.0;      // 1.0;
    private final double SW_MISMATCH = -8.0;  //-1.0/3.0;
    private final double SW_GAP;       //-1.0-1.0/3.0;
    private final double SW_GAP_EXTEND; //-1.0/.0;

    private final boolean DEBUG;

    private final static double LOG_ONE_HALF = -Math.log10(2.0);
    private final static List<Allele> noCall = new ArrayList<Allele>(); // used to noCall all genotypes until the exact model is applied

    final HashMap<Integer, ArrayList<Event>> allEventDictionary;

    private class Event {
        public VariantContext vc;
        public Allele refAllele;
        public Allele altAllele;
        public int index;
    }

    public GenotypingEngine( final boolean DEBUG, final double gop, final double gcp ) {
        this.DEBUG = DEBUG;
        SW_GAP = -1.0 * gop;
        SW_GAP_EXTEND = -1.0 * gcp;
        noCall.add(Allele.NO_CALL);
        allEventDictionary = new HashMap<Integer, ArrayList<Event>>();
    }

    public ArrayList<VariantContext> alignAndAssignGenotypeLikelihoods( final GenomeLocParser genomeLocParser, final ArrayList<Haplotype> allHaplotypes, final Set<Haplotype> bestHaplotypes, final byte[] ref, final GenomeLoc loc, final GenomeLoc window, final HashMap<String,Double[][]> haplotypeLikelihoodMatrixMap ) {

        final HashMap<Integer, ArrayList<Event>> bestEventDictionary = new HashMap<Integer, ArrayList<Event>>(); // These are the events we will actually be genotyping
        final ArrayList<VariantContext> returnCallContexts = new ArrayList<VariantContext>();

        // Create the dictionary of all genotype-able events sorted by start location and annotated with the originating haplotype index
        if( DEBUG ) { System.out.println(" ========  Top Haplotypes ======== "); }
        populateEventDictionary(bestEventDictionary, bestHaplotypes, ref, loc, window, false);

        // walk along the haplotype in genomic order, genotyping the events as the come up
        final ArrayList<Integer> sortedKeySet = new ArrayList<Integer>();
        sortedKeySet.addAll(bestEventDictionary.keySet());
        Collections.sort(sortedKeySet);
        for( final Integer key : sortedKeySet ) {
            final ArrayList<Event> allEventList = allEventDictionary.get(key);
            final ArrayList<Event> bestEventList = bestEventDictionary.get(key);

            // Gather together all the VCs at this start location and use VariantContextUtils to merge the alleles together
            final ArrayList<VariantContext> vcsToGenotype = new ArrayList<VariantContext>();
            for( final Event e : bestEventList ) {
                vcsToGenotype.add(e.vc);
            }

            // For multi-allelic deletion records the alleles need to be expanded to match the length of the longest allele
            // Also add reference events to the dictionary for every haplotype that doesn't have any event to make the for loop below very easy to write
            final VariantContext mergedVC = VariantContextUtils.simpleMerge(genomeLocParser, vcsToGenotype, null, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, VariantContextUtils.GenotypeMergeType.UNSORTED, false, false, null, false, false);
            correctAndExpandEventListWithRefEvents(allEventList, mergedVC, haplotypeLikelihoodMatrixMap.values().iterator().next()[0].length);
            vcsToGenotype.clear();
            vcsToGenotype.add( mergedVC );

            for( final VariantContext vcToGenotype : vcsToGenotype ) {
                if( DEBUG ) { System.out.println("Genotyping event at " + key + " with alleles: " + vcToGenotype.getAlleles()); }
                final Map<String, Genotype> genotypes = new LinkedHashMap<String, Genotype>();
                // Grab the genotype likelihoods from the appropriate places in the haplotype likelihood matrix -- calculation performed independently per sample
                for( final String sample : haplotypeLikelihoodMatrixMap.keySet() ) {
                    final double[] genotypeLikelihoods = new double[(vcToGenotype.getAlleles().size() * (vcToGenotype.getAlleles().size()+1)) / 2];
                    final Double[][] haplotypeLikelihoodMatrix = haplotypeLikelihoodMatrixMap.get( sample );
                    int glIndex = 0;
                    for( int iii = 0; iii < vcToGenotype.getAlleles().size(); iii++ ) {
                        for( int jjj = 0; jjj <= iii; jjj++ ) {
                            double likelihood = Double.NEGATIVE_INFINITY;
                            final Pair<Allele, Allele> allelePair1 = new Pair<Allele, Allele>(vcToGenotype.getReference(), vcToGenotype.getAlleles().get(jjj));
                            final Pair<Allele, Allele> allelePair2 = new Pair<Allele, Allele>(vcToGenotype.getReference(), vcToGenotype.getAlleles().get(iii));

                            // Loop through all haplotype pairs and find the max likelihood that has this given combination of events on the pair of haplotypes
                            for( final Event e1 : allEventList ) {
                                if( allelePair1.equals( new Pair<Allele, Allele>(e1.refAllele, e1.altAllele) ) ) {
                                    for( final Event e2 : allEventList ) {
                                        if( allelePair2.equals( new Pair<Allele, Allele>(e2.refAllele, e2.altAllele) ) ) {
                                            likelihood = Math.max( likelihood, haplotypeLikelihoodMatrix[e1.index][e2.index] );
                                        }
                                    }
                                }
                            }

                            if( Double.isInfinite(likelihood) ) {
                                throw new ReviewedStingException("Infinite likelihood detected. Maybe the correct event wasn't found in the event dictionary.");
                            }

                            genotypeLikelihoods[glIndex++] = likelihood;
                            if( DEBUG ) { System.out.println(iii + ", " + jjj + ": " + likelihood); }
                        }
                    }
                    final HashMap<String, Object> attributes = new HashMap<String, Object>();
                    attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods((
                            new MultiallelicGenotypeLikelihoods(sample, vcToGenotype.getAlleles(), genotypeLikelihoods, 40)).getLikelihoods()));
                    genotypes.put(sample, new Genotype(sample, noCall, Genotype.NO_NEG_LOG_10PERROR, null, attributes, false));
                }
                returnCallContexts.add( VariantContext.modifyGenotypes(vcToGenotype, genotypes) );
            }
        }

        return returnCallContexts;
    }

    public void createEventDictionaryAndFilterBadHaplotypes( final ArrayList<Haplotype> allHaplotypes, final byte[] ref, final GenomeLoc loc, final GenomeLoc window ) {
        allEventDictionary.clear();
        populateEventDictionary(allEventDictionary, allHaplotypes, ref, loc, window, true);
    }

    private void populateEventDictionary(final HashMap<Integer, ArrayList<Event>> eventDictionary, final Collection<Haplotype> haplotypes, final byte[] ref, final GenomeLoc loc, final GenomeLoc window, final boolean filterBadHaplotypes ) {
        int hIndex = 0;
        final HashSet<Haplotype> haplotypesToRemove = new HashSet<Haplotype>();
        for( final Haplotype h : haplotypes ) {

            // Align the haplotype to the reference
            final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( ref, h.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
            if( DEBUG ) {
                System.out.println( h.toString() );
                System.out.println( "Cigar = " + swConsensus.getCigar() );
            }
            if( swConsensus.getCigar().getReadLength() < 10 ) {
                if( DEBUG ) { System.out.println("Filtered!"); }
                if( filterBadHaplotypes ) { haplotypesToRemove.add(h); }
                continue; // Protection against SW failures
            }

            // Walk along the alignment and turn any difference from the reference into an event
            final ArrayList<VariantContext> vcs = generateVCsFromAlignment(swConsensus, ref, h.bases, loc);

            if( vcs == null || tooManyClusteredVariantsOnHaplotype( vcs ) ) { // too many variants on this haplotype means it wasn't assembled very well
                if( DEBUG ) { System.out.println("Filtered!"); }
                if( filterBadHaplotypes ) { haplotypesToRemove.add(h); }
                continue; // Protection against SW failures
            }

            // Separate all events into dictionaries partitioned by start location
            for( final VariantContext vc : vcs ) {
                if( vc.getStart() >= window.getStart() && vc.getStart() <= window.getStop() ) {
                    ArrayList<Event> eventList = eventDictionary.get(vc.getStart());
                    if(eventList == null) { // haven't seen this start location yet, so need to create a new list
                        eventList = new ArrayList<Event>();
                        eventDictionary.put(vc.getStart(), eventList);
                    }
                    final Event e = new Event();
                    e.vc = vc; // BUGBUG: probably don't need to keep all the VC's around anymore or even create them in the first place
                    e.refAllele = vc.getReference();
                    e.altAllele = vc.getAlternateAlleles().get(0); // This vc is guaranteed to have only one alternate allele
                    if( vc.getAlternateAlleles().size() != 1 ) {
                        throw new ReviewedStingException("BUG: smith-waterman derived variant context has more than one alternate allele!");
                    }
                    e.index = hIndex;
                    eventList.add(e);
                }
            }

            hIndex++;
        }
        if( filterBadHaplotypes ) { haplotypes.removeAll(haplotypesToRemove); }
    }

    private boolean tooManyClusteredVariantsOnHaplotype( final ArrayList<VariantContext> vcs ) {
        final int clusterSize = 60;
        final int threshold = 4;

        for(int iii = 0; iii < vcs.size() - threshold + 1; iii++) {
            final int size = vcs.get(iii + threshold - 1).getStart() - vcs.get(iii).getStart();
            if( size <= clusterSize ) {
                return true;
            }
        }

        return false;
    }

    private void correctAndExpandEventListWithRefEvents( final ArrayList<Event> inputEvents, final VariantContext mergedVC, final int maxHaplotypeIndex ) {
        if( DEBUG ) { System.out.println( "! Merged event ref = " + mergedVC.getReference() ); }
        for(int iii = 0; iii < maxHaplotypeIndex; iii++) {
            Event myEvent = null;
            for( final Event e : inputEvents) {
                if( e.index == iii ) { myEvent = e; }
            }
            if( myEvent == null ) { // this is a ref haplotype so add a ref event
                final Event e = new Event();
                e.vc = null;
                e.refAllele = mergedVC.getReference();
                e.altAllele = mergedVC.getReference();
                e.index = iii;
                inputEvents.add(e);
            } else { // might need to correct the alleles because of the potential merging of multiallelic records
                if( DEBUG ) { System.out.println( "My Event = " + myEvent.refAllele + "/" + myEvent.altAllele ); }
                if( mergedVC.getAlternateAlleles().size() > 1 && !myEvent.refAllele.equals(mergedVC.getReference()) ) {
                    final int suffixSize = mergedVC.getReference().getBases().length - myEvent.refAllele.length();
                    if( suffixSize > 0 ) {
                        if( (myEvent.vc.isSimpleInsertion() || myEvent.vc.isSimpleDeletion()) && mergedVC.isMixed() ) { // the special case of combining a SNP and an in/del (one has padded reference but the other doesn't)
                            myEvent.altAllele = Allele.extend(myEvent.altAllele, Arrays.copyOfRange(mergedVC.getReference().getBases(), myEvent.refAllele.getBases().length + 1, myEvent.refAllele.getBases().length + suffixSize));
                            myEvent.altAllele = Allele.create(mergedVC.getReference().getBaseString().charAt(0) + myEvent.altAllele.getBaseString());
                        } else {
                            myEvent.altAllele = Allele.extend(myEvent.altAllele, Arrays.copyOfRange(mergedVC.getReference().getBases(), myEvent.refAllele.getBases().length, myEvent.refAllele.getBases().length + suffixSize));
                        }
                        myEvent.refAllele = mergedVC.getReference();
                    }
                }
                if( DEBUG ) { System.out.println( "--> Updated Event = " + myEvent.refAllele + "/" + myEvent.altAllele ); }
            }
        }
    }

    private ArrayList<VariantContext> generateVCsFromAlignment( final SWPairwiseAlignment swConsensus, final byte[] ref, final byte[] read, final GenomeLoc loc ) {
        final ArrayList<VariantContext> vcs = new ArrayList<VariantContext>();

        int refPos = swConsensus.getAlignmentStart2wrt1();
        if( refPos==0 ) { return null; } // Protection against SW failures
        if( swConsensus.getCigar().toString().contains("S") ) { return null; } // Protection against SW failures
        int readPos = 0;
        final int lookAhead = 3;

        for( final CigarElement ce : swConsensus.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
                case I:
                {
                    final byte[] insertionBases = Arrays.copyOfRange( read, readPos, readPos + elementLength);
                    boolean allN = true;
                    for( byte b : insertionBases ) {
                        if( b != (byte) 'N' ) {
                            allN = false;
                            break;
                        }
                    }
                    if( !allN ) {
                        final ArrayList<Allele> alleles = new ArrayList<Allele>();
                        alleles.add( Allele.create(Allele.NULL_ALLELE_STRING, true));
                        alleles.add( Allele.create(insertionBases, false));
                        if( DEBUG ) { System.out.println("@ " + (loc.getStart() + refPos - 1) + " > Insertion: " + alleles); }
                        vcs.add(new VariantContext("HaplotypeCaller", loc.getContig(), loc.getStart() + refPos - 1, loc.getStart() + refPos - 1,
                                alleles, VariantContext.NO_GENOTYPES, VariantContext.NO_NEG_LOG_10PERROR, null, null, ref[refPos-1]));
                    }
                    readPos += elementLength;
                    break;
                }
                case S:
                {
                    readPos += elementLength;
                    break;
                }
                case D:
                {
                    final byte[] deletionBases = Arrays.copyOfRange( ref, refPos, refPos + elementLength);
                    final ArrayList<Allele> alleles = new ArrayList<Allele>();
                    alleles.add( Allele.create(deletionBases, true) );
                    alleles.add( Allele.create(Allele.NULL_ALLELE_STRING, false) );
                    if( DEBUG ) { System.out.println( "@ " + (loc.getStart() + refPos - 1) + " Deletion: " + alleles); }
                    vcs.add( new VariantContext("HaplotypeCaller", loc.getContig(), loc.getStart() + refPos - 1, loc.getStart() + refPos + elementLength - 1,
                            alleles, VariantContext.NO_GENOTYPES, VariantContext.NO_NEG_LOG_10PERROR, null, null, ref[refPos-1]) );
                    refPos += elementLength;
                    break;
                }
                case M:
                {
                    int numSinceMismatch = -1;
                    int stopOfMismatch = -1;
                    int startOfMismatch = -1;
                    int refPosStartOfMismatch = -1;
                    for( int iii = 0; iii < elementLength; iii++ ) {
                        if( ref[refPos] != read[readPos] && read[readPos] != ((byte) 'N') ) {
                            // SNP or start of possible MNP
                            if( stopOfMismatch == -1 ) {
                                startOfMismatch = readPos;
                                stopOfMismatch = readPos;
                                numSinceMismatch = 0;
                                refPosStartOfMismatch = refPos;
                            } else {
                                stopOfMismatch = readPos;
                            }
                        }

                        if( stopOfMismatch != -1) {
                            numSinceMismatch++;
                        }

                        if( numSinceMismatch > lookAhead || (iii == elementLength - 1 && stopOfMismatch != -1) ) {
                            final byte[] refBases = Arrays.copyOfRange( ref, refPosStartOfMismatch, refPosStartOfMismatch + (stopOfMismatch - startOfMismatch) + 1 );
                            final byte[] mismatchBases = Arrays.copyOfRange( read, startOfMismatch, stopOfMismatch + 1 );
                            final ArrayList<Allele> alleles = new ArrayList<Allele>();
                            alleles.add( Allele.create( refBases, true ) );
                            alleles.add( Allele.create( mismatchBases, false ) );
                            if( DEBUG ) { System.out.println( "@ " + (loc.getStart() + refPosStartOfMismatch) + " > SNP/MNP: " + alleles); }
                            vcs.add( new VariantContext("HaplotypeCaller", loc.getContig(), loc.getStart() + refPosStartOfMismatch,
                                    loc.getStart() + refPosStartOfMismatch + (stopOfMismatch - startOfMismatch), alleles,
                                    VariantContext.NO_GENOTYPES, VariantContext.NO_NEG_LOG_10PERROR, null, null) );
                            numSinceMismatch = -1;
                            stopOfMismatch = -1;
                            startOfMismatch = -1;
                            refPosStartOfMismatch = -1;
                        }

                        refPos++;
                        readPos++;
                    }
                    break;
                }

                case N:
                case H:
                case P:
                default:
                    throw new ReviewedStingException( "Unsupported cigar operator: " + ce.getOperator() );
            }
        }

        if( DEBUG && vcs.size() == 0 ) {
            System.out.println("> Reference!");
        }

        return vcs;
    }







    //////////////////////////////////////////
    //
    // Code for debug output follows
    //
    //////////////////////////////////////////        

    public void alignAllHaplotypes( final ArrayList<Haplotype> haplotypes, final byte[] ref, final GenomeLoc loc, final StingSAMFileWriter writer, final SAMRecord exampleRead ) {

        int iii = 0;
        for( final Haplotype h : haplotypes ) {
            final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( ref, h.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
            exampleRead.setReadName("Haplotype" + iii);
            exampleRead.setReadBases(h.bases);
            exampleRead.setAlignmentStart(loc.getStart() + swConsensus.getAlignmentStart2wrt1());
            exampleRead.setCigar(swConsensus.getCigar());
            final byte[] quals = new byte[h.bases.length];
            Arrays.fill(quals, (byte) 25);
            exampleRead.setBaseQualities(quals);
            writer.addAlignment(exampleRead);
            iii++;
        }
                
    }

    public void alignAllReads( final Pair<Haplotype,Haplotype> bestTwoHaplotypes, final byte[] ref, final GenomeLoc loc, final ConstrainedMateFixingManager manager, final ArrayList<SAMRecord> reads, final double[][] likelihoods ) {

        final SWPairwiseAlignment swConsensus0 = new SWPairwiseAlignment( ref, bestTwoHaplotypes.first.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        final SWPairwiseAlignment swConsensus1 = new SWPairwiseAlignment( ref, bestTwoHaplotypes.second.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        final Consensus consensus0 = new Consensus(bestTwoHaplotypes.first.bases, swConsensus0.getCigar(), swConsensus0.getAlignmentStart2wrt1());
        final Consensus consensus1 = new Consensus(bestTwoHaplotypes.second.bases, swConsensus1.getCigar(), swConsensus1.getAlignmentStart2wrt1());

        int iii = 0;
        for( final SAMRecord read : reads ) {
            final Consensus bestConsensus = ( likelihoods[iii][0] > likelihoods[iii][1] ? consensus0 : consensus1 );
            final AlignedRead aRead = new AlignedRead( read );
            bestConsensus.cigar = AlignmentUtils.leftAlignIndel(bestConsensus.cigar, ref, bestConsensus.str, bestConsensus.positionOnReference, bestConsensus.positionOnReference);
            Pair<Integer, Integer> altAlignment = findBestOffset(bestConsensus.str, aRead, loc.getStart());
            updateRead(bestConsensus.cigar, bestConsensus.positionOnReference, altAlignment.second, aRead, loc.getStart());
            aRead.finalizeUpdate();
            manager.addRead(aRead.getRead(), true);
            iii++;
        }
    }

    // private functions copied from IndelRealigner
    private Pair<Integer, Integer> findBestOffset(final byte[] ref, final AlignedRead read, final int leftmostIndex) {

        // optimization: try the most likely alignment first (to get a low score to beat)
        int originalAlignment = read.getOriginalAlignmentStart() - leftmostIndex;
        int bestScore = mismatchQualitySumIgnoreCigar(read, ref, originalAlignment, Integer.MAX_VALUE);
        int bestIndex = originalAlignment;

        // optimization: we can't get better than 0, so we can quit now
        if ( bestScore == 0 )
            return new Pair<Integer, Integer>(bestIndex, 0);

        // optimization: the correct alignment shouldn't be too far from the original one (or else the read wouldn't have aligned in the first place)
        for ( int i = 0; i < originalAlignment; i++ ) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if ( bestScore == 0 )
                return new Pair<Integer, Integer>(bestIndex, 0);
        }

        final int maxPossibleStart = ref.length - read.getReadLength();
        for ( int i = originalAlignment + 1; i <= maxPossibleStart; i++ ) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if ( bestScore == 0 )
                return new Pair<Integer, Integer>(bestIndex, 0);
        }

        return new Pair<Integer, Integer>(bestIndex, bestScore);
    }

    private static int mismatchQualitySumIgnoreCigar(final AlignedRead aRead, final byte[] refSeq, int refIndex, int quitAboveThisValue) {
        final byte[] readSeq = aRead.getReadBases();
        final byte[] quals = aRead.getBaseQualities();
        int sum = 0;
        for (int readIndex = 0 ; readIndex < readSeq.length ; refIndex++, readIndex++ ) {
            if ( refIndex >= refSeq.length ) {
                sum += 99;
                // optimization: once we pass the threshold, stop calculating
                if ( sum > quitAboveThisValue )
                    return sum;
            } else {
                byte refChr = refSeq[refIndex];
                byte readChr = readSeq[readIndex];
                if ( !BaseUtils.isRegularBase(readChr) || !BaseUtils.isRegularBase(refChr) )
                    continue; // do not count Ns/Xs/etc ?
                if ( readChr != refChr ) {
                    sum += (int)quals[readIndex];
                    // optimization: once we pass the threshold, stop calculating
                    if ( sum > quitAboveThisValue )
                        return sum;
                }
            }
        }
        return sum;
    }

    private boolean updateRead(final Cigar altCigar, final int altPosOnRef, final int myPosOnAlt, final AlignedRead aRead, final int leftmostIndex) {
        Cigar readCigar = new Cigar();

        // special case: there is no indel
        if ( altCigar.getCigarElements().size() == 1 ) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
            readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            aRead.setCigar(readCigar);
            return true;
        }

        CigarElement altCE1 = altCigar.getCigarElement(0);
        CigarElement altCE2 = altCigar.getCigarElement(1);

        int leadingMatchingBlockLength = 0; // length of the leading M element or 0 if the leading element is I

        CigarElement indelCE;
        if ( altCE1.getOperator() == CigarOperator.I  ) {
            indelCE=altCE1;
            if ( altCE2.getOperator() != CigarOperator.M  ) {
                return false;
            }
        }
        else {
            if ( altCE1.getOperator() != CigarOperator.M  ) {
                return false;
            }
            if ( altCE2.getOperator() == CigarOperator.I  || altCE2.getOperator() == CigarOperator.D ) {
                indelCE=altCE2;
            } else {
                return false;
            }
            leadingMatchingBlockLength = altCE1.getLength();
        }

        // the easiest thing to do is to take each case separately
        int endOfFirstBlock = altPosOnRef + leadingMatchingBlockLength;
        boolean sawAlignmentStart = false;

        // for reads starting before the indel
        if ( myPosOnAlt < endOfFirstBlock) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
            sawAlignmentStart = true;

            // for reads ending before the indel
            if ( myPosOnAlt + aRead.getReadLength() <= endOfFirstBlock) {
                //readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
                //aRead.setCigar(readCigar);
                aRead.setCigar(null); // reset to original alignment
                return true;
            }
            readCigar.add(new CigarElement(endOfFirstBlock - myPosOnAlt, CigarOperator.M));
        }

        // forward along the indel
        //int indelOffsetOnRef = 0, indelOffsetOnRead = 0;
        if ( indelCE.getOperator() == CigarOperator.I ) {
            // for reads that end in an insertion
            if ( myPosOnAlt + aRead.getReadLength() < endOfFirstBlock + indelCE.getLength() ) {
                int partialInsertionLength = myPosOnAlt + aRead.getReadLength() - endOfFirstBlock;
                // if we also started inside the insertion, then we need to modify the length
                if ( !sawAlignmentStart )
                    partialInsertionLength = aRead.getReadLength();
                readCigar.add(new CigarElement(partialInsertionLength, CigarOperator.I));
                aRead.setCigar(readCigar);
                return true;
            }

            // for reads that start in an insertion
            if ( !sawAlignmentStart && myPosOnAlt < endOfFirstBlock + indelCE.getLength() ) {
                aRead.setAlignmentStart(leftmostIndex + endOfFirstBlock);
                readCigar.add(new CigarElement(indelCE.getLength() - (myPosOnAlt - endOfFirstBlock), CigarOperator.I));
                //indelOffsetOnRead = myPosOnAlt - endOfFirstBlock;
                sawAlignmentStart = true;
            } else if ( sawAlignmentStart ) {
                readCigar.add(indelCE);
                //indelOffsetOnRead = indelCE.getLength();
            }
        } else if ( indelCE.getOperator() == CigarOperator.D ) {
            if ( sawAlignmentStart )
                readCigar.add(indelCE);
            //indelOffsetOnRef = indelCE.getLength();
        }

        // for reads that start after the indel
        if ( !sawAlignmentStart ) {
            //aRead.setAlignmentStart(leftmostIndex + myPosOnAlt + indelOffsetOnRef - indelOffsetOnRead);
            //readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            //aRead.setCigar(readCigar);
            aRead.setCigar(null); // reset to original alignment
            return true;
        }

        int readRemaining = aRead.getReadBases().length;
        for ( CigarElement ce : readCigar.getCigarElements() ) {
            if ( ce.getOperator() != CigarOperator.D )
                readRemaining -= ce.getLength();
        }
        if ( readRemaining > 0 )
            readCigar.add(new CigarElement(readRemaining, CigarOperator.M));
        aRead.setCigar(readCigar);

        return true;
    }


    // private classes copied from IndelRealigner
    private class AlignedRead {
        private final SAMRecord read;
        private byte[] readBases = null;
        private byte[] baseQuals = null;
        private Cigar newCigar = null;
        private int newStart = -1;
        private int mismatchScoreToReference = 0;
        private long alignerMismatchScore = 0;

        public AlignedRead(SAMRecord read) {
            this.read = read;
            mismatchScoreToReference = 0;
        }

        public SAMRecord getRead() {
               return read;
        }

        public int getReadLength() {
            return readBases != null ? readBases.length : read.getReadLength();
        }

        public byte[] getReadBases() {
            if ( readBases == null )
                getUnclippedBases();
            return readBases;
        }

        public byte[] getBaseQualities() {
            if ( baseQuals == null )
                getUnclippedBases();
            return baseQuals;
        }

        // pull out the bases that aren't clipped out
        private void getUnclippedBases() {
            readBases = new byte[getReadLength()];
            baseQuals = new byte[getReadLength()];
            byte[] actualReadBases = read.getReadBases();
            byte[] actualBaseQuals = read.getBaseQualities();
            int fromIndex = 0, toIndex = 0;

            for ( CigarElement ce : read.getCigar().getCigarElements() ) {
                int elementLength = ce.getLength();
                switch ( ce.getOperator() ) {
                    case S:
                        fromIndex += elementLength;
                        break;
                    case M:
                    case I:
                        System.arraycopy(actualReadBases, fromIndex, readBases, toIndex, elementLength);
                        System.arraycopy(actualBaseQuals, fromIndex, baseQuals, toIndex, elementLength);
                        fromIndex += elementLength;
                        toIndex += elementLength;
                    default:
                        break;
                }
            }

            // if we got clipped, trim the array
            if ( fromIndex != toIndex ) {
                byte[] trimmedRB = new byte[toIndex];
                byte[] trimmedBQ = new byte[toIndex];
                System.arraycopy(readBases, 0, trimmedRB, 0, toIndex);
                System.arraycopy(baseQuals, 0, trimmedBQ, 0, toIndex);
                readBases = trimmedRB;
                baseQuals = trimmedBQ;
            }
        }

        public Cigar getCigar() {
            return (newCigar != null ? newCigar : read.getCigar());
        }

        public void setCigar(Cigar cigar) {
            setCigar(cigar, true);
        }

        // tentatively sets the new Cigar, but it needs to be confirmed later
        public void setCigar(Cigar cigar, boolean fixClippedCigar) {
            if ( cigar == null ) {
                newCigar = null;
                return;
            }

            if ( fixClippedCigar && getReadBases().length < read.getReadLength() )
                cigar = reclipCigar(cigar);

            // no change?
            if ( read.getCigar().equals(cigar) ) {
                newCigar = null;
                return;
            }

            // no indel?
            String str = cigar.toString();
            if ( !str.contains("D") && !str.contains("I") ) {
                //    newCigar = null;
                //    return;
            }

            newCigar = cigar;
        }

        // pull out the bases that aren't clipped out
        private Cigar reclipCigar(Cigar cigar) {
            return reclipCigar(cigar, read);
        }

        private boolean isClipOperator(CigarOperator op) {
            return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
        }

        protected Cigar reclipCigar(Cigar cigar, SAMRecord read) {
            ArrayList<CigarElement> elements = new ArrayList<CigarElement>();

            int i = 0;
            int n = read.getCigar().numCigarElements();
            while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
                elements.add(read.getCigar().getCigarElement(i++));

            elements.addAll(cigar.getCigarElements());

            i++;
            while ( i < n && !isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
                i++;

            while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
                elements.add(read.getCigar().getCigarElement(i++));

            return new Cigar(elements);
        }

        // tentatively sets the new start, but it needs to be confirmed later
        public void setAlignmentStart(int start) {
            newStart = start;
        }

        public int getAlignmentStart() {
            return (newStart != -1 ? newStart : read.getAlignmentStart());
        }

        public int getOriginalAlignmentStart() {
            return read.getAlignmentStart();
        }

        // finalizes the changes made.
        // returns true if this record actually changes, false otherwise
        public boolean finalizeUpdate() {
            // if we haven't made any changes, don't do anything
            if ( newCigar == null )
                return false;
            if ( newStart == -1 )
                newStart = read.getAlignmentStart();

            read.setCigar(newCigar);
            read.setAlignmentStart(newStart);

            return true;
        }

        public void setMismatchScoreToReference(int score) {
            mismatchScoreToReference = score;
        }

        public int getMismatchScoreToReference() {
            return mismatchScoreToReference;
        }

        public void setAlignerMismatchScore(long score) {
            alignerMismatchScore = score;
        }

        public long getAlignerMismatchScore() {
            return alignerMismatchScore;
        }
    }

    private static class Consensus {
        public final byte[] str;
        public final ArrayList<Pair<Integer, Integer>> readIndexes;
        public final int positionOnReference;
        public int mismatchSum;
        public Cigar cigar;

        public Consensus(byte[] str, Cigar cigar, int positionOnReference) {
            this.str = str;
            this.cigar = cigar;
            this.positionOnReference = positionOnReference;
            mismatchSum = 0;
            readIndexes = new ArrayList<Pair<Integer, Integer>>();
        }

        @Override
        public boolean equals(Object o) {
            return ( this == o || (o instanceof Consensus && Arrays.equals(this.str,(((Consensus)o).str)) ) );
        }

        public boolean equals(Consensus c) {
            return ( this == c || Arrays.equals(this.str,c.str) ) ;
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(this.str);
        }
    }


}