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

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.CigarElement;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

public class GenotypingEngine {

    private final boolean DEBUG;
    private final static List<Allele> noCall = new ArrayList<Allele>(); // used to noCall all genotypes until the exact model is applied

    public GenotypingEngine( final boolean DEBUG ) {
        this.DEBUG = DEBUG;
        noCall.add(Allele.NO_CALL);
    }

    @Requires({"refLoc.containsP(activeRegionWindow)", "haplotypes.size() > 0"})
    public ArrayList<VariantContext> assignGenotypeLikelihoodsAndCallEvents( final UnifiedGenotyperEngine UG_engine, final ArrayList<Haplotype> haplotypes, final byte[] ref, final GenomeLoc refLoc, 
                                                                             final GenomeLoc activeRegionWindow, final GenomeLocParser genomeLocParser ) {

        // Prepare the list of haplotype indices to genotype
        final ArrayList<Allele> allelesToGenotype = new ArrayList<Allele>();

        if( DEBUG ) { System.out.println("=== Best Haplotypes ==="); }
        for( final Haplotype h : haplotypes ) {
            allelesToGenotype.add( Allele.create(h.getBases(), h.isReference()) );
        }
        final int numHaplotypes = haplotypes.size();

        // Grab the genotype likelihoods from the appropriate places in the haplotype likelihood matrix -- calculation performed independently per sample
        final GenotypesContext genotypes = GenotypesContext.create(haplotypes.get(0).getSampleKeySet().size());
        for( final String sample : haplotypes.get(0).getSampleKeySet() ) { // BUGBUG: assume all haplotypes saw the same samples
            final double[] genotypeLikelihoods = new double[numHaplotypes * (numHaplotypes+1) / 2];
            final double[][] haplotypeLikelihoodMatrix = LikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods(haplotypes, sample);
            int glIndex = 0;
            for( int iii = 0; iii < numHaplotypes; iii++ ) {
                for( int jjj = 0; jjj <= iii; jjj++ ) {
                    genotypeLikelihoods[glIndex++] = haplotypeLikelihoodMatrix[iii][jjj]; // for example: AA,AB,BB,AC,BC,CC
                }
            }
            final HashMap<String, Object> attributes = new HashMap<String, Object>();
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods(genotypeLikelihoods));
            genotypes.add(new Genotype(sample, noCall, Genotype.NO_LOG10_PERROR, null, attributes, false));
        }
        final VariantCallContext call = UG_engine.calculateGenotypes(new VariantContextBuilder().loc(activeRegionWindow).alleles(allelesToGenotype).genotypes(genotypes).make(), UG_engine.getUAC().GLmodel);
        if( call == null ) { return new ArrayList<VariantContext>(); } // exact model says that the call confidence is below the specified confidence threshold so nothing to do here

        // Prepare the list of haplotypes that need to be run through Smith-Waterman for output to VCF
        final ArrayList<Haplotype> haplotypesToRemove = new ArrayList<Haplotype>();
        for( final Haplotype h : haplotypes ) {
            if( call.getAllele(h.getBases()) == null ) { // exact model removed this allele from the list so no need to run SW and output to VCF
                haplotypesToRemove.add(h);
            }
        }
        haplotypes.removeAll(haplotypesToRemove);

        return callEventsFromHaplotypes(call, haplotypes, ref, refLoc, activeRegionWindow, genomeLocParser);
    }

    @Requires({"refLoc.containsP(activeRegionWindow)", "haplotypes.size() > 0"})
    public ArrayList<VariantContext> callEventsFromHaplotypes( final VariantCallContext call, final ArrayList<Haplotype> haplotypes, final byte[] ref, final GenomeLoc refLoc,
                                                               final GenomeLoc activeRegionWindow, final GenomeLocParser genomeLocParser ) {

        final ArrayList<VariantContext> returnVCs = new ArrayList<VariantContext>();
        // Smith-Waterman parameters originally copied from IndelRealigner
        final double SW_MATCH = 5.0;
        final double SW_MISMATCH = -8.0;
        final double SW_GAP = -30.0;
        final double SW_GAP_EXTEND = -1.4;

        // Run Smith-Waterman on each called haplotype to figure out what events need to be written out in a VCF file
        final TreeSet<Integer> startPosKeySet = new TreeSet<Integer>();
        final ArrayList<HashMap<Integer,VariantContext>> fullEventDictionary = new ArrayList<HashMap<Integer, VariantContext>>();
        int count = 0;
        for( final Haplotype h : haplotypes ) {
            final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( ref, h.getBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
            if( DEBUG ) {
                System.out.println( h.toString() );
                System.out.println( "> Cigar = " + swConsensus.getCigar() );
            }
            // Walk along the alignment and turn any difference from the reference into an event
            final HashMap<Integer,VariantContext> eventMap = generateVCsFromAlignment( swConsensus, ref, h.getBases(), refLoc, "HC" + count++ );
            fullEventDictionary.add(eventMap);
            startPosKeySet.addAll(eventMap.keySet());
        }
        
        // Create the VC merge priority list
        final ArrayList<String> priorityList = new ArrayList<String>();
        for( int iii = 0; iii < haplotypes.size(); iii++ ) {
            priorityList.add("HC" + iii);
        }
        
        // Walk along each position in the key set and create each event to be outputted
        for( final int loc : startPosKeySet ) {
            if( loc >= activeRegionWindow.getStart() && loc <= activeRegionWindow.getStop() ) {
                final ArrayList<VariantContext> eventsAtThisLoc = new ArrayList<VariantContext>();
                for( final HashMap<Integer,VariantContext> eventMap : fullEventDictionary ) {
                    final VariantContext vc = eventMap.get(loc);
                    if( vc != null && !containsVCWithMatchingAlleles(eventsAtThisLoc, vc) ) {
                        eventsAtThisLoc.add(vc);
                    }
                }
                
                // Create the allele mapping object which maps the original haplotype alleles to the alleles present in just this event
                final ArrayList<ArrayList<Integer>> alleleMapper = createAlleleMapper( loc, eventsAtThisLoc, fullEventDictionary );

                // Merge the event to find a common reference representation
                final VariantContext mergedVC = VariantContextUtils.simpleMerge(genomeLocParser, eventsAtThisLoc, priorityList, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, VariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);

                if( DEBUG ) {
                    System.out.println("Genotyping event at " + loc + " with alleles = " + mergedVC.getAlleles());
                    System.out.println("Event/haplotype allele mapping = " + alleleMapper);
                }

                // Grab the genotype likelihoods from the appropriate places in the haplotype likelihood matrix -- calculation performed independently per sample
                final GenotypesContext genotypes = GenotypesContext.create(haplotypes.get(0).getSampleKeySet().size());
                for( final String sample : haplotypes.get(0).getSampleKeySet() ) { // BUGBUG: assume all haplotypes saw the same samples
                    final int numHaplotypes = alleleMapper.size();
                    final double[] genotypeLikelihoods = new double[numHaplotypes * (numHaplotypes+1) / 2];
                    final double[][] haplotypeLikelihoodMatrix = LikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods(haplotypes, sample, alleleMapper);
                    int glIndex = 0;
                    for( int iii = 0; iii < numHaplotypes; iii++ ) {
                        for( int jjj = 0; jjj <= iii; jjj++ ) {
                            genotypeLikelihoods[glIndex++] = haplotypeLikelihoodMatrix[iii][jjj]; // for example: AA,AB,BB,AC,BC,CC
                        }
                    }
                    final HashMap<String, Object> attributes = new HashMap<String, Object>();
                    attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods(genotypeLikelihoods));

                    // using the allele mapping object translate the haplotype allele into the event allele
                    final ArrayList<Allele> eventAllelesForSample = findEventAllelesInSample( mergedVC.getAlleles(), call.getAlleles(), call.getGenotype(sample).getAlleles(), alleleMapper );
                    genotypes.add(new Genotype(sample, eventAllelesForSample, Genotype.NO_LOG10_PERROR, null, attributes, false));
                }
                returnVCs.add(new VariantContextBuilder(mergedVC).log10PError(call.getLog10PError()).genotypes(genotypes).make());
            }
        }

        return returnVCs;
    }

    @Requires({"fullEventDictionary.size() <= eventsAtThisLoc.size() + 1"})
    @Ensures({"result.size() == eventsAtThisLoc.size() + 1"})
    protected static ArrayList<ArrayList<Integer>> createAlleleMapper( final int loc, final ArrayList<VariantContext> eventsAtThisLoc, final ArrayList<HashMap<Integer,VariantContext>> fullEventDictionary ) {
        final ArrayList<ArrayList<Integer>> alleleMapper = new ArrayList<ArrayList<Integer>>();
        final ArrayList<Integer> refList = new ArrayList<Integer>();
        for( int jjj = 0; jjj < fullEventDictionary.size(); jjj++ ) {
            if( fullEventDictionary.get(jjj).get(loc) == null ) {
                refList.add(jjj);
            }
        }
        alleleMapper.add(refList);
        for( int iii = 0; iii < eventsAtThisLoc.size(); iii++ ) {
            final ArrayList<Integer> list = new ArrayList<Integer>();
            for( int jjj = 0; jjj < fullEventDictionary.size(); jjj++ ) {
                if( fullEventDictionary.get(jjj).get(loc) != null && fullEventDictionary.get(jjj).get(loc).hasSameAllelesAs(eventsAtThisLoc.get(iii)) ) {
                    list.add(jjj);
                }
            }
            alleleMapper.add(list);
        }
        return alleleMapper;
    }

    @Ensures({"result.size() == haplotypeAllelesForSample.size()"})
    protected static ArrayList<Allele> findEventAllelesInSample( final List<Allele> eventAlleles, final List<Allele> haplotypeAlleles, final List<Allele> haplotypeAllelesForSample, final ArrayList<ArrayList<Integer>> alleleMapper ) {
        final ArrayList<Allele> eventAllelesForSample = new ArrayList<Allele>();
        for( final Allele a : haplotypeAllelesForSample ) {
            final int haplotypeIndex = haplotypeAlleles.indexOf(a);
            for( int iii = 0; iii < alleleMapper.size(); iii++ ) {
                final ArrayList<Integer> haplotypeIndicesList = alleleMapper.get(iii);
                if( haplotypeIndicesList.contains(haplotypeIndex) ) {
                    eventAllelesForSample.add(eventAlleles.get(iii));
                    break;
                }
            }
        }
        return eventAllelesForSample;
    }

    protected static boolean containsVCWithMatchingAlleles( final List<VariantContext> list, final VariantContext vcToTest ) {
        for( final VariantContext vc : list ) {
            if( vc.hasSameAllelesAs(vcToTest) ) {
                return true;
            }
        }
        return false;
    }

    protected static HashMap<Integer,VariantContext> generateVCsFromAlignment( final SWPairwiseAlignment swConsensus, final byte[] ref, final byte[] alignment, final GenomeLoc refLoc, final String sourceNameToAdd ) {
        final HashMap<Integer,VariantContext> vcs = new HashMap<Integer,VariantContext>();

        int refPos = swConsensus.getAlignmentStart2wrt1();
        if( refPos < 0 ) { return null; } // Protection against SW failures
        int alignmentPos = 0;
        final int lookAhead = 5; // Used to create MNPs out of nearby SNPs on the same haplotype

        for( final CigarElement ce : swConsensus.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
                case I:
                    final byte[] insertionBases = Arrays.copyOfRange( alignment, alignmentPos, alignmentPos + elementLength );
                    boolean allN = true;
                    for( final byte b : insertionBases ) {
                        if( b != (byte) 'N' ) {
                            allN = false;
                            break;
                        }
                    }
                    if( !allN ) {
                        final ArrayList<Allele> insertionAlleles = new ArrayList<Allele>();
                        insertionAlleles.add( Allele.create(Allele.NULL_ALLELE_STRING, true));
                        insertionAlleles.add(Allele.create(insertionBases, false));
                        final int insertionStart = refLoc.getStart() + refPos - 1;
                        vcs.put(insertionStart, new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), insertionStart, insertionStart, insertionAlleles).referenceBaseForIndel(ref[refPos-1]).make());
                    }
                    alignmentPos += elementLength;
                    break;
                case S:
                    alignmentPos += elementLength;
                    break;
                case D:
                    final byte[] deletionBases = Arrays.copyOfRange( ref, refPos, refPos + elementLength );
                    final ArrayList<Allele> deletionAlleles = new ArrayList<Allele>();
                    deletionAlleles.add( Allele.create(deletionBases, true) );
                    deletionAlleles.add(Allele.create(Allele.NULL_ALLELE_STRING, false));
                    final int deletionStart = refLoc.getStart() + refPos - 1;
                    vcs.put(deletionStart, new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), deletionStart, deletionStart + elementLength, deletionAlleles).referenceBaseForIndel(ref[refPos-1]).make());
                    refPos += elementLength;
                    break;
                case M:
                    int numSinceMismatch = -1;
                    int stopOfMismatch = -1;
                    int startOfMismatch = -1;
                    int refPosStartOfMismatch = -1;
                    for( int iii = 0; iii < elementLength; iii++ ) {
                        if( ref[refPos] != alignment[alignmentPos] && alignment[alignmentPos] != ((byte) 'N') ) {
                            // SNP or start of possible MNP
                            if( stopOfMismatch == -1 ) {
                                startOfMismatch = alignmentPos;
                                stopOfMismatch = alignmentPos;
                                numSinceMismatch = 0;
                                refPosStartOfMismatch = refPos;
                            } else {
                                stopOfMismatch = alignmentPos;
                            }
                        }
                        if( stopOfMismatch != -1) {
                            numSinceMismatch++;
                        }
                        if( numSinceMismatch > lookAhead || (iii == elementLength - 1 && stopOfMismatch != -1) ) {
                            final byte[] refBases = Arrays.copyOfRange( ref, refPosStartOfMismatch, refPosStartOfMismatch + (stopOfMismatch - startOfMismatch) + 1 );
                            final byte[] mismatchBases = Arrays.copyOfRange( alignment, startOfMismatch, stopOfMismatch + 1 );
                            final ArrayList<Allele> snpAlleles = new ArrayList<Allele>();
                            snpAlleles.add( Allele.create( refBases, true ) );
                            snpAlleles.add(Allele.create(mismatchBases, false));
                            final int snpStart = refLoc.getStart() + refPosStartOfMismatch;
                            vcs.put(snpStart, new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), snpStart, snpStart + (stopOfMismatch - startOfMismatch), snpAlleles).make());
                            numSinceMismatch = -1;
                            stopOfMismatch = -1;
                            startOfMismatch = -1;
                            refPosStartOfMismatch = -1;
                        }
                        refPos++;
                        alignmentPos++;
                    }
                    break;
                case N:
                case H:
                case P:
                default:
                    throw new ReviewedStingException( "Unsupported cigar operator created during SW alignment: " + ce.getOperator() );
            }
        }
        return vcs;
    }
}