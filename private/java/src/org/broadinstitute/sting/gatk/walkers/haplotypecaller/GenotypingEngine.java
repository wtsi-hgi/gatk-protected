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
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

public class GenotypingEngine {

    // Smith-Waterman parameters originally copied from IndelRealigner
    private final double SW_MATCH = 5.0;
    private final double SW_MISMATCH = -8.0;
    private final double SW_GAP; // command-line argument in HaplotypeCaller walker
    private final double SW_GAP_EXTEND; // command-line argument in HaplotypeCaller walker
    private final boolean DEBUG;

    private final static List<Allele> noCall = new ArrayList<Allele>(); // used to noCall all genotypes until the exact model is applied

    public GenotypingEngine( final boolean DEBUG, final double gop, final double gcp ) {
        this.DEBUG = DEBUG;
        SW_GAP = -1.0 * gop;
        SW_GAP_EXTEND = -1.0 * gcp;
        noCall.add(Allele.NO_CALL);
    }

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

    public ArrayList<VariantContext> callEventsFromHaplotypes( final VariantCallContext call, final ArrayList<Haplotype> haplotypes, final byte[] ref, final GenomeLoc refLoc,
                                                               final GenomeLoc activeRegionWindow, final GenomeLocParser genomeLocParser ) {

        final ArrayList<VariantContext> returnVCs = new ArrayList<VariantContext>();

        // Run Smith-Waterman on each called haplotype to figure out what events need to be written out in a VCF file
        final TreeSet<Integer> startPosKeySet = new TreeSet<Integer>();
        final ArrayList<HashMap<Integer,VariantContext>> eventDictionary = new ArrayList<HashMap<Integer, VariantContext>>();
        int count = 0;
        for( final Haplotype h : haplotypes ) {
            final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( ref, h.getBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
            if( DEBUG ) {
                System.out.println( h.toString() );
                System.out.println( "> Cigar = " + swConsensus.getCigar() );
            }
            // Walk along the alignment and turn any difference from the reference into an event
            final HashMap<Integer,VariantContext> eventMap = generateVCsFromAlignment( swConsensus, ref, h.getBases(), refLoc, "HC" + count++ );
            eventDictionary.add( eventMap );
            startPosKeySet.addAll( eventMap.keySet() );
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
                for( final HashMap<Integer,VariantContext> eventMap : eventDictionary ) {
                    final VariantContext vc = eventMap.get(loc);
                    if( vc != null && !containsVCWithMatchingAlleles(eventsAtThisLoc, vc) ) {
                        eventsAtThisLoc.add(vc);
                    }
                }
                
                // Create the allele mapping object which maps the original haplotype alleles to the alleles present in just this event
                final ArrayList<ArrayList<Integer>> alleleMapping = new ArrayList<ArrayList<Integer>>();
                final ArrayList<Integer> refList = new ArrayList<Integer>();
                for( int jjj = 0; jjj < eventDictionary.size(); jjj++ ) {
                    if( eventDictionary.get(jjj).get(loc) == null ) {
                        refList.add(jjj);
                    }
                }
                alleleMapping.add(refList);
                for( int iii = 0; iii < eventsAtThisLoc.size(); iii++ ) {
                    final ArrayList<Integer> list = new ArrayList<Integer>();
                    for( int jjj = 0; jjj < eventDictionary.size(); jjj++ ) {
                        if( eventDictionary.get(jjj).get(loc) != null && eventDictionary.get(jjj).get(loc).hasSameAllelesAs(eventsAtThisLoc.get(iii)) ) {
                            list.add(jjj);
                        }
                    }
                    alleleMapping.add(list);
                }
                
                final VariantContext mergedVC = VariantContextUtils.simpleMerge(genomeLocParser, eventsAtThisLoc, priorityList, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, VariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);

                if( DEBUG ) {
                    System.out.println("Genotyping event at " + loc + " with alleles = " + mergedVC.getAlleles());
                    System.out.println("Event/haplotype allele mapping = " + alleleMapping);
                }

                // Grab the genotype likelihoods from the appropriate places in the haplotype likelihood matrix -- calculation performed independently per sample
                final GenotypesContext genotypes = GenotypesContext.create(haplotypes.get(0).getSampleKeySet().size());
                for( final String sample : haplotypes.get(0).getSampleKeySet() ) { // BUGBUG: assume all haplotypes saw the same samples
                    final int numHaplotypes = alleleMapping.size();
                    final double[] genotypeLikelihoods = new double[numHaplotypes * (numHaplotypes+1) / 2];
                    final double[][] haplotypeLikelihoodMatrix = LikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods(haplotypes, sample, alleleMapping);
                    int glIndex = 0;
                    for( int iii = 0; iii < numHaplotypes; iii++ ) {
                        for( int jjj = 0; jjj <= iii; jjj++ ) {
                            genotypeLikelihoods[glIndex++] = haplotypeLikelihoodMatrix[iii][jjj]; // for example: AA,AB,BB,AC,BC,CC
                        }
                    }
                    final HashMap<String, Object> attributes = new HashMap<String, Object>();
                    attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods(genotypeLikelihoods));

                    // using the haplotype mapping object translate the haplotype allele into the event allele
                    final ArrayList<Allele> eventAllelesForSample = new ArrayList<Allele>();
                    final List<Allele> haplotypeAllelesForSample = call.getGenotype(sample).getAlleles();
                    for( final Allele a : haplotypeAllelesForSample ) {
                        final int haplotypeIndex = call.getAlleles().indexOf(a);
                        for( int iii = 0; iii < alleleMapping.size(); iii++ ) {
                            final ArrayList<Integer> haplotypeIndicesList = alleleMapping.get(iii);
                            if( haplotypeIndicesList.contains(haplotypeIndex) ) {
                                eventAllelesForSample.add(mergedVC.getAlleles().get(iii));
                                break;
                            }
                        }
                    }
                    genotypes.add(new Genotype(sample, eventAllelesForSample, Genotype.NO_LOG10_PERROR, null, attributes, false));
                }
                returnVCs.add(new VariantContextBuilder(mergedVC).log10PError(call.getLog10PError()).genotypesNoValidation(genotypes).make());
            }
        }

        return returnVCs;
    } 
    
    private boolean containsVCWithMatchingAlleles( final List<VariantContext> list, final VariantContext vcToTest ) {
        for( final VariantContext vc : list ) {
            if( vc.hasSameAllelesAs(vcToTest) ) {
                return true;
            }
        }
        return false;
    }
    
    private HashMap<Integer,VariantContext> generateVCsFromAlignment( final SWPairwiseAlignment swConsensus, final byte[] ref, final byte[] read, final GenomeLoc refLoc, final String sourceName ) {
        final HashMap<Integer,VariantContext> vcs = new HashMap<Integer,VariantContext>();

        int refPos = swConsensus.getAlignmentStart2wrt1();
        if( refPos < 0 ) { return null; } // Protection against SW failures
        int readPos = 0;
        final int lookAhead = 5;

        for( final CigarElement ce : swConsensus.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            //final ArrayList<Allele> alleles = new ArrayList<Allele>();
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
                        final int start = refLoc.getStart() + refPos - 1;
                        vcs.put(start, new VariantContextBuilder(sourceName, refLoc.getContig(), start, start, alleles).referenceBaseForIndel(ref[refPos-1]).make());
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
                    final int start = refLoc.getStart() + refPos - 1;
                    vcs.put(start, new VariantContextBuilder(sourceName, refLoc.getContig(), start, start + elementLength, alleles).referenceBaseForIndel(ref[refPos-1]).make());
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
                            final int start = refLoc.getStart() + refPosStartOfMismatch;
                            vcs.put(start, new VariantContextBuilder(sourceName, refLoc.getContig(), start, start + (stopOfMismatch - startOfMismatch), alleles).make());
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
                    throw new ReviewedStingException( "Unsupported cigar operator created during SW alignment: " + ce.getOperator() );
            }
        }

        if( DEBUG && vcs.size() == 0 ) {
            System.out.println("> Reference!");
        }

        return vcs;
    }
}