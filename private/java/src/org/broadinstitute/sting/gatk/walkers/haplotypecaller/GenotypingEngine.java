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
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

public class GenotypingEngine {

    private final boolean DEBUG;
    private final int MNP_LOOK_AHEAD;
    private final boolean OUTPUT_FULL_HAPLOTYPE_SEQUENCE;
    private final static List<Allele> noCall = new ArrayList<Allele>(); // used to noCall all genotypes until the exact model is applied

    public GenotypingEngine( final boolean DEBUG, final int MNP_LOOK_AHEAD, final boolean OUTPUT_FULL_HAPLOTYPE_SEQUENCE ) {
        this.DEBUG = DEBUG;
        this.MNP_LOOK_AHEAD = MNP_LOOK_AHEAD;
        this.OUTPUT_FULL_HAPLOTYPE_SEQUENCE = OUTPUT_FULL_HAPLOTYPE_SEQUENCE;
        noCall.add(Allele.NO_CALL);
    }

    // This function is the streamlined approach, currently not being used
    @Requires({"refLoc.containsP(activeRegionWindow)", "haplotypes.size() > 0"})
    public List<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>> assignGenotypeLikelihoodsAndCallHaplotypeEvents( final UnifiedGenotyperEngine UG_engine, final ArrayList<Haplotype> haplotypes, final byte[] ref, final GenomeLoc refLoc,
                                                                                                                             final GenomeLoc activeRegionWindow, final GenomeLocParser genomeLocParser ) {
        // Prepare the list of haplotype indices to genotype
        final ArrayList<Allele> allelesToGenotype = new ArrayList<Allele>();

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
            //if( DEBUG ) { System.out.println(sample + " --> " + Arrays.toString(genotypeLikelihoods)); }
            final HashMap<String, Object> attributes = new HashMap<String, Object>();
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods(genotypeLikelihoods));
            genotypes.add(new Genotype(sample, noCall, Genotype.NO_LOG10_PERROR, null, attributes, false));
        }
        final VariantCallContext call = UG_engine.calculateGenotypes(new VariantContextBuilder().loc(activeRegionWindow).alleles(allelesToGenotype).genotypes(genotypes).make(), UG_engine.getUAC().GLmodel);
        if( call == null ) { return Collections.emptyList(); } // exact model says that the call confidence is below the specified confidence threshold so nothing to do here

        // Prepare the list of haplotypes that need to be run through Smith-Waterman for output to VCF
        final ArrayList<Haplotype> haplotypesToRemove = new ArrayList<Haplotype>();
        for( final Haplotype h : haplotypes ) {
            if( call.getAllele(h.getBases()) == null ) { // exact model removed this allele from the list so no need to run SW and output to VCF
                haplotypesToRemove.add(h);
            }
        }
        haplotypes.removeAll(haplotypesToRemove);

        if( OUTPUT_FULL_HAPLOTYPE_SEQUENCE ) {
            final List<Pair<VariantContext, HashMap<Allele, ArrayList<Haplotype>>>> returnVCs = new ArrayList<Pair<VariantContext, HashMap<Allele, ArrayList<Haplotype>>>>();
            // set up the default 1-to-1 haplotype mapping object
            final HashMap<Allele,ArrayList<Haplotype>> haplotypeMapping = new HashMap<Allele,ArrayList<Haplotype>>();
            for( final Haplotype h : haplotypes ) {
                final ArrayList<Haplotype> list = new ArrayList<Haplotype>();
                list.add(h);
                haplotypeMapping.put(call.getAllele(h.getBases()), list);
            }
            returnVCs.add( new Pair<VariantContext, HashMap<Allele, ArrayList<Haplotype>>>(call,haplotypeMapping) );
            return returnVCs;
        }

        final ArrayList<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>> returnCalls = new ArrayList<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>>();

        // Using the cigar from each called haplotype figure out what events need to be written out in a VCF file
        final TreeSet<Integer> startPosKeySet = new TreeSet<Integer>();
        int count = 0;
        if( DEBUG ) { System.out.println("=== Best Haplotypes ==="); }
        for( final Haplotype h : haplotypes ) {
            if( DEBUG ) {
                System.out.println( h.toString() );
                System.out.println( "> Cigar = " + h.getCigar() );
            }
            // Walk along the alignment and turn any difference from the reference into an event
            h.setEventMap( generateVCsFromAlignment( h.getAlignmentStartHapwrtRef(), h.getCigar(), ref, h.getBases(), refLoc, "HC" + count++, MNP_LOOK_AHEAD ) );
            startPosKeySet.addAll(h.getEventMap().keySet());
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
                for( final Haplotype h : haplotypes ) {
                    final HashMap<Integer,VariantContext> eventMap = h.getEventMap();
                    final VariantContext vc = eventMap.get(loc);
                    if( vc != null && !containsVCWithMatchingAlleles(eventsAtThisLoc, vc) ) {
                        eventsAtThisLoc.add(vc);
                    }
                }
                
                // Create the allele mapping object which maps the original haplotype alleles to the alleles present in just this event
                final ArrayList<ArrayList<Haplotype>> alleleMapper = createAlleleMapper( loc, eventsAtThisLoc, haplotypes );

                // Merge the event to find a common reference representation
                final VariantContext mergedVC = VariantContextUtils.simpleMerge(genomeLocParser, eventsAtThisLoc, priorityList, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, VariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);

                final HashMap<Allele, ArrayList<Haplotype>> alleleHashMap = new HashMap<Allele, ArrayList<Haplotype>>();
                int aCount = 0;
                for( final Allele a : mergedVC.getAlleles() ) {
                    alleleHashMap.put(a, alleleMapper.get(aCount++)); // BUGBUG: needs to be cleaned up and merged with alleleMapper
                }

                if( DEBUG ) {
                    System.out.println("Genotyping event at " + loc + " with alleles = " + mergedVC.getAlleles());
                    //System.out.println("Event/haplotype allele mapping = " + alleleMapper);
                }

                // Grab the genotype likelihoods from the appropriate places in the haplotype likelihood matrix -- calculation performed independently per sample
                final GenotypesContext myGenotypes = GenotypesContext.create(haplotypes.get(0).getSampleKeySet().size());
                for( final String sample : haplotypes.get(0).getSampleKeySet() ) { // BUGBUG: assume all haplotypes saw the same samples
                    final int myNumHaplotypes = alleleMapper.size();
                    final double[] genotypeLikelihoods = new double[myNumHaplotypes * (myNumHaplotypes+1) / 2];
                    final double[][] haplotypeLikelihoodMatrix = LikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods(haplotypes, sample, alleleMapper);
                    int glIndex = 0;
                    for( int iii = 0; iii < myNumHaplotypes; iii++ ) {
                        for( int jjj = 0; jjj <= iii; jjj++ ) {
                            genotypeLikelihoods[glIndex++] = haplotypeLikelihoodMatrix[iii][jjj]; // for example: AA,AB,BB,AC,BC,CC
                        }
                    }
                    final HashMap<String, Object> attributes = new HashMap<String, Object>();
                    attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods(genotypeLikelihoods));

                    // using the allele mapping object translate the haplotype allele into the event allele
                    myGenotypes.add(new Genotype(sample, findEventAllelesInSample(mergedVC.getAlleles(), call.getAlleles(), call.getGenotype(sample).getAlleles(), alleleMapper, haplotypes),
                            Genotype.NO_LOG10_PERROR, null, attributes, loc != startPosKeySet.first()));
                }
                returnCalls.add( new Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>(
                                 new VariantContextBuilder(mergedVC).log10PError(call.getLog10PError()).genotypes(myGenotypes).make(), alleleHashMap) );
            }
        }
        return returnCalls;
    }

    @Requires({"refLoc.containsP(activeRegionWindow)", "haplotypes.size() > 0"})
    public List<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>> assignGenotypeLikelihoodsAndCallIndependentEvents( final UnifiedGenotyperEngine UG_engine, final ArrayList<Haplotype> haplotypes, final byte[] ref, final GenomeLoc refLoc,
                                                                                                                               final GenomeLoc activeRegionWindow, final GenomeLocParser genomeLocParser, final ArrayList<VariantContext> activeAllelesToGenotype ) {

        final ArrayList<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>> returnCalls = new ArrayList<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>>();

        // Using the cigar from each called haplotype figure out what events need to be written out in a VCF file
        final TreeSet<Integer> startPosKeySet = new TreeSet<Integer>();
        int count = 0;
        if( DEBUG ) { System.out.println("=== Best Haplotypes ==="); }
        for( final Haplotype h : haplotypes ) {
            // Walk along the alignment and turn any difference from the reference into an event
            h.setEventMap( generateVCsFromAlignment( h.getAlignmentStartHapwrtRef(), h.getCigar(), ref, h.getBases(), refLoc, "HC" + count++, MNP_LOOK_AHEAD ) );
            if( activeAllelesToGenotype.isEmpty() ) { startPosKeySet.addAll(h.getEventMap().keySet()); }
            if( DEBUG ) {
                System.out.println( h.toString() );
                System.out.println( "> Cigar = " + h.getCigar() );
                System.out.println( ">> Events = " + h.getEventMap().values());
            }
        }
        if( !activeAllelesToGenotype.isEmpty() ) { // we are in GGA mode!
            for( final VariantContext compVC : activeAllelesToGenotype ) {
                startPosKeySet.add( compVC.getStart() );
            }
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
                if( activeAllelesToGenotype.isEmpty() ) {
                    for( final Haplotype h : haplotypes ) {
                        final HashMap<Integer,VariantContext> eventMap = h.getEventMap();
                        final VariantContext vc = eventMap.get(loc);
                        if( vc != null && !containsVCWithMatchingAlleles(eventsAtThisLoc, vc) ) {
                            eventsAtThisLoc.add(vc);
                        }
                    }
                } else { // we are in GGA mode!
                    for( final VariantContext compVC : activeAllelesToGenotype ) {
                        if( compVC.getStart() == loc ) {
                            priorityList.clear();
                            int alleleCount = 0;
                            for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                                HashSet<Allele> alleleSet = new HashSet<Allele>(2);
                                alleleSet.add(compVC.getReference());
                                alleleSet.add(compAltAllele);
                                priorityList.add("Allele" + alleleCount);
                                eventsAtThisLoc.add(new VariantContextBuilder(compVC.subContextFromSamples(compVC.getSampleNames(), alleleSet)).source("Allele"+alleleCount).make());
                                alleleCount++;
                            }
                        }
                    }
                }

                //if( eventsAtThisLoc.isEmpty() ) { continue; }

                // Create the allele mapping object which maps the original haplotype alleles to the alleles present in just this event
                final ArrayList<ArrayList<Haplotype>> alleleMapper = createAlleleMapper( loc, eventsAtThisLoc, haplotypes );

                // Merge the event to find a common reference representation
                final VariantContext mergedVC = VariantContextUtils.simpleMerge(genomeLocParser, eventsAtThisLoc, priorityList, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, VariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);

                final HashMap<Allele, ArrayList<Haplotype>> alleleHashMap = new HashMap<Allele, ArrayList<Haplotype>>();
                int aCount = 0;
                for( final Allele a : mergedVC.getAlleles() ) {
                    alleleHashMap.put(a, alleleMapper.get(aCount++)); // BUGBUG: needs to be cleaned up and merged with alleleMapper
                }

                if( DEBUG ) {
                    System.out.println("Genotyping event at " + loc + " with alleles = " + mergedVC.getAlleles());
                    //System.out.println("Event/haplotype allele mapping = " + alleleMapper);
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
                    genotypes.add(new Genotype(sample, noCall, Genotype.NO_LOG10_PERROR, null, attributes, false));
                }
                final VariantCallContext call = UG_engine.calculateGenotypes(new VariantContextBuilder(mergedVC).genotypes(genotypes).make(), UG_engine.getUAC().GLmodel);

                if( call != null ) {
                    returnCalls.add( new Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>(call, alleleHashMap) );
                }
            }
        }
        return returnCalls;
    }

    @Requires({"haplotypes.size() >= eventsAtThisLoc.size() + 1"})
    @Ensures({"result.size() == eventsAtThisLoc.size() + 1"})
    protected static ArrayList<ArrayList<Haplotype>> createAlleleMapper( final int loc, final ArrayList<VariantContext> eventsAtThisLoc, final ArrayList<Haplotype> haplotypes ) {
        final ArrayList<ArrayList<Haplotype>> alleleMapper = new ArrayList<ArrayList<Haplotype>>();
        final ArrayList<Haplotype> refList = new ArrayList<Haplotype>();
        for( final Haplotype h : haplotypes ) {
            if( h.getEventMap().get(loc) == null ) {
                refList.add(h);
            }
        }
        alleleMapper.add(refList);
        for( final VariantContext vcAtThisLoc : eventsAtThisLoc ) {
            final ArrayList<Haplotype> list = new ArrayList<Haplotype>();
            for( final Haplotype h : haplotypes ) {
                if( h.getEventMap().get(loc) != null && h.getEventMap().get(loc).hasSameAllelesAs(vcAtThisLoc) ) {
                    list.add(h);
                }
            }
            alleleMapper.add(list);
        }
        return alleleMapper;
    }

    @Ensures({"result.size() == haplotypeAllelesForSample.size()"})
    protected static List<Allele> findEventAllelesInSample( final List<Allele> eventAlleles, final List<Allele> haplotypeAlleles, final List<Allele> haplotypeAllelesForSample, final ArrayList<ArrayList<Haplotype>> alleleMapper, final ArrayList<Haplotype> haplotypes ) {
        if( haplotypeAllelesForSample.contains(Allele.NO_CALL) ) { return noCall; }
        final ArrayList<Allele> eventAllelesForSample = new ArrayList<Allele>();
        for( final Allele a : haplotypeAllelesForSample ) {
            final Haplotype haplotype = haplotypes.get(haplotypeAlleles.indexOf(a));
            for( int iii = 0; iii < alleleMapper.size(); iii++ ) {
                final ArrayList<Haplotype> mappedHaplotypes = alleleMapper.get(iii);
                if( mappedHaplotypes.contains(haplotype) ) {
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

    protected static HashMap<Integer,VariantContext> generateVCsFromAlignment( final int alignmentStartHapwrtRef, final Cigar cigar, final byte[] ref, final byte[] alignment, final GenomeLoc refLoc, final String sourceNameToAdd, final int MNP_LOOK_AHEAD ) {
        final HashMap<Integer,VariantContext> vcs = new HashMap<Integer,VariantContext>();

        int refPos = alignmentStartHapwrtRef;
        if( refPos < 0 ) { return null; } // Protection against SW failures
        int alignmentPos = 0;

        for( final CigarElement ce : cigar.getCigarElements() ) {
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
                        if( numSinceMismatch > MNP_LOOK_AHEAD || (iii == elementLength - 1 && stopOfMismatch != -1) ) {
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