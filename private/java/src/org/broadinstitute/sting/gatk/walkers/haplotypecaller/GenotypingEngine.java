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

import net.sf.samtools.CigarElement;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

public class GenotypingEngine {

    // Smith-Waterman parameters originally copied from IndelRealigner
    private final double SW_MATCH = 5.0;      // 1.0;
    private final double SW_MISMATCH = -8.0;  //-1.0/3.0;
    private final double SW_GAP;       //-1.0-1.0/3.0;
    private final double SW_GAP_EXTEND; //-1.0/.0;

    private final boolean DEBUG;

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

    public ArrayList<VariantContext> alignAndAssignGenotypeLikelihoods( final GenomeLocParser genomeLocParser, final ArrayList<Haplotype> allHaplotypes, final ArrayList<Haplotype> bestHaplotypes, final byte[] ref, final GenomeLoc loc, final GenomeLoc window, final HashMap<String,Double[][]> haplotypeLikelihoodMatrixMap, final ArrayList<VariantContext> allelesToGenotype ) {

        final HashMap<Integer, ArrayList<Event>> bestEventDictionary = new HashMap<Integer, ArrayList<Event>>(); // These are the events we will actually be genotyping
        final ArrayList<VariantContext> returnCallContexts = new ArrayList<VariantContext>();

        // Create the dictionary of all genotype-able events sorted by start location and annotated with the originating haplotype index
        if( DEBUG ) { System.out.println(" ========  Top Haplotypes ======== "); }
        if( allelesToGenotype != null && !allelesToGenotype.isEmpty() ) { // GENOTYPE_GIVEN_ALLELES mode
            bestHaplotypes.clear();
            bestHaplotypes.add( allHaplotypes.get(0) ); // the reference haplotype
            populateEventDictionary(bestEventDictionary, bestHaplotypes, ref, loc, window, true, allelesToGenotype);
        } else {
            populateEventDictionary(bestEventDictionary, bestHaplotypes, ref, loc, window, false, null);
        }

        // walk along the haplotype in genomic order, genotyping the events as they come up
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

            if( DEBUG ) { System.out.println("Genotyping event at " + key + " with alleles: " + mergedVC.getAlleles()); }

            final GenotypesContext genotypes = GenotypesContext.create(haplotypeLikelihoodMatrixMap.keySet().size());
            // Grab the genotype likelihoods from the appropriate places in the haplotype likelihood matrix -- calculation performed independently per sample
            for( final String sample : haplotypeLikelihoodMatrixMap.keySet() ) {
                final double[] genotypeLikelihoods = new double[(mergedVC.getAlleles().size() * (mergedVC.getAlleles().size()+1)) / 2];
                final Double[][] haplotypeLikelihoodMatrix = haplotypeLikelihoodMatrixMap.get( sample );
                int glIndex = 0;
                for( int iii = 0; iii < mergedVC.getAlleles().size(); iii++ ) {
                    for( int jjj = 0; jjj <= iii; jjj++ ) {
                        double likelihood = Double.NEGATIVE_INFINITY;
                        final Pair<Allele, Allele> allelePair1 = new Pair<Allele, Allele>(mergedVC.getReference(), mergedVC.getAlleles().get(jjj));
                        final Pair<Allele, Allele> allelePair2 = new Pair<Allele, Allele>(mergedVC.getReference(), mergedVC.getAlleles().get(iii));

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
                    }
                }
                final HashMap<String, Object> attributes = new HashMap<String, Object>();
                attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods((genotypeLikelihoods)));
                genotypes.add(new Genotype(sample, noCall, Genotype.NO_LOG10_PERROR, null, attributes, false));
            }
            returnCallContexts.add( new VariantContextBuilder(mergedVC).genotypes(genotypes).make() );
        }

        return returnCallContexts;
    }

    public void createEventDictionaryAndFilterBadHaplotypes( final ArrayList<Haplotype> allHaplotypes, final byte[] ref, final GenomeLoc loc, final GenomeLoc window, final ArrayList<VariantContext> allelesToGenotype ) {
        allEventDictionary.clear();
        populateEventDictionary( allEventDictionary, allHaplotypes, ref, loc, window, true, allelesToGenotype );
    }

    private void populateEventDictionary( final HashMap<Integer, ArrayList<Event>> eventDictionary, final Collection<Haplotype> haplotypes, final byte[] paddedRef, 
                                          final GenomeLoc paddedLoc, final GenomeLoc window, final boolean filterBadHaplotypes, final ArrayList<VariantContext> allelesToGenotype ) {
        int hIndex = 0;
        int sizeRefHaplotype = 0;
        final ArrayList<Haplotype> haplotypesToRemove = new ArrayList<Haplotype>();
        final ArrayList<Haplotype> haplotypesToAdd = new ArrayList<Haplotype>();
        if( filterBadHaplotypes ) {
            for( final Haplotype h : haplotypes ) {

                // Align the haplotype to the reference
                final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( paddedRef, h.getBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
                if( DEBUG ) {
                    System.out.println();
                    System.out.println( h.toString() );
                    System.out.println( "Cigar = " + swConsensus.getCigar() );
                }
                // Now check for smith-waterman failures and remove the haplotype from the candidate list
                if( swConsensus.getCigar().getReadLength() < 10 || swConsensus.getCigar().toString().contains("S") ) {
                    if( DEBUG ) { System.out.println("Filtered! -SW failure"); }
                    haplotypesToRemove.add(h);
                    continue; // Protection against SW failures
                }
                if( hIndex == 0 ) { // book-keeping, the first haplotype is always the full reference
                    sizeRefHaplotype = swConsensus.getCigar().getReadLength();
                } else if ( Math.max(swConsensus.getCigar().getReadLength(), swConsensus.getCigar().getReferenceLength() ) < 0.65 * sizeRefHaplotype ) {
                    // Sometimes the assembly produces singleton, short paths that should be filtered out, need to perform "tip clipping" on the assembly graph
                    if( DEBUG ) { System.out.println("Filtered! -too short"); }
                   haplotypesToRemove.add(h);
                    continue; // Protection against assembly failures
                }
                hIndex++;

                // GENOTYPE_GIVEN_ALLELES mode
                // Artificially insert the provided alternate alleles into the discovered haplotype to genotype them
                if( allelesToGenotype != null && !allelesToGenotype.isEmpty() ) {
                    for( final VariantContext vc : allelesToGenotype ) {
                        for( final Allele a : vc.getAlternateAlleles() ) {
                            haplotypesToAdd.add( new Haplotype( h.insertAllele(vc.getReference(), a, vc.getStart() - paddedLoc.getStart(), swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar()) ) );
                        }
                    }                
                }
            }
            haplotypes.removeAll( haplotypesToRemove );
            for( final Haplotype h : haplotypesToAdd ) {
                if( !haplotypes.contains(h) ) { haplotypes.add(h); }
            }
        }        

        haplotypesToRemove.clear();
        hIndex = 0;
        for( final Haplotype h : haplotypes ) {

            final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( paddedRef, h.getBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
            if( DEBUG ) {
                System.out.println();
                System.out.println( h.toString() );
                System.out.println( "Cigar = " + swConsensus.getCigar() );
            }
            // Extending the haplotype to include extra reference bases might have cause a smith-waterman failure so remove the haplotype from the candidate list
            if( swConsensus.getCigar().getReadLength() < 10 || swConsensus.getCigar().toString().contains("S") ) {
                if( DEBUG ) { System.out.println("Filtered! -SW failure"); }
                if( filterBadHaplotypes ) { haplotypesToRemove.add(h); }
                continue; // Protection against SW failures
            }

            // Walk along the alignment and turn any difference from the reference into an event
            final ArrayList<VariantContext> vcs = generateVCsFromAlignment(swConsensus, paddedRef, h.getBases(), paddedLoc);

            if( vcs == null || tooManyClusteredVariantsOnHaplotype( vcs ) ) { // too many variants on this haplotype means it wasn't assembled very well
                if( DEBUG ) { System.out.println("Filtered! -too complex"); }
                if( filterBadHaplotypes ) { haplotypesToRemove.add(h); }
                continue; // Protection against SW failures
            }

            // Separate all events into a dictionary partitioned by start location
            for( final VariantContext vc : vcs ) {
                if( DEBUG ) { System.out.println( ">> " + vc); }
                if( vc.getStart() >= window.getStart() && vc.getStart() <= window.getStop() ) {
                    ArrayList<Event> eventList = eventDictionary.get(vc.getStart());
                    if( eventList == null ) { // haven't seen this start location yet, so need to create a new list
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
        if( filterBadHaplotypes ) {
            haplotypes.removeAll( haplotypesToRemove );
        }
    }

    private boolean tooManyClusteredVariantsOnHaplotype( final ArrayList<VariantContext> vcs ) {
        // Turning off clustered variants intrinsic filter for now
        /*
        final int clusterSize = 60;
        final int threshold = 4;

        for(int iii = 0; iii < vcs.size() - threshold + 1; iii++) {
            final int size = vcs.get(iii + threshold - 1).getStart() - vcs.get(iii).getStart();
            if( size <= clusterSize ) {
                return true;
            }
        }
        */

        return false;
    }

    private void correctAndExpandEventListWithRefEvents( final ArrayList<Event> inputEvents, final VariantContext mergedVC, final int maxHaplotypeIndex ) {
        for(int iii = 0; iii < maxHaplotypeIndex; iii++) {
            Event myEvent = null;
            for( final Event e : inputEvents ) {
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
                for( final Event e : inputEvents ) {
                    if( e.index == iii ) {
                        myEvent = e;
                        //if( DEBUG ) { System.out.println( "My Event = " + myEvent.refAllele + "/" + myEvent.altAllele ); }
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
                        //if( DEBUG ) { System.out.println( "--> Updated Event = " + myEvent.refAllele + "/" + myEvent.altAllele ); }
                    }
                }
            }
        }
    }

    private ArrayList<VariantContext> generateVCsFromAlignment( final SWPairwiseAlignment swConsensus, final byte[] ref, final byte[] read, final GenomeLoc loc ) {
        final ArrayList<VariantContext> vcs = new ArrayList<VariantContext>();

        int refPos = swConsensus.getAlignmentStart2wrt1();
        if( refPos < 0 ) { return null; } // Protection against SW failures
        int readPos = 0;
        final int lookAhead = 0; // BUGBUG: There is trouble with multiallelics when one allele is an MNP and the other allele is one of the component SNPs

        for( final CigarElement ce : swConsensus.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
                case I:
                {
                    final byte[] insertionBases = Arrays.copyOfRange( read, readPos, readPos + elementLength);
                    boolean allN = true;
                    for( final byte b : insertionBases ) {
                        if( b != (byte) 'N' ) {
                            allN = false;
                            break;
                        }
                    }
                    if( !allN ) {
                        final ArrayList<Allele> alleles = new ArrayList<Allele>();
                        alleles.add( Allele.create(Allele.NULL_ALLELE_STRING, true));
                        alleles.add( Allele.create(insertionBases, false));
                        vcs.add(new VariantContextBuilder("HaplotypeCaller", loc.getContig(), loc.getStart() + refPos - 1, loc.getStart() + refPos - 1,
                                alleles).referenceBaseForIndel(ref[refPos-1]).make());
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
                    vcs.add( new VariantContextBuilder("HaplotypeCaller", loc.getContig(), loc.getStart() + refPos - 1, loc.getStart() + refPos + elementLength - 1,
                            alleles).referenceBaseForIndel(ref[refPos-1]).make());
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
                            vcs.add( new VariantContextBuilder("HaplotypeCaller", loc.getContig(), loc.getStart() + refPosStartOfMismatch,
                                    loc.getStart() + refPosStartOfMismatch + (stopOfMismatch - startOfMismatch), alleles).make());
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
}