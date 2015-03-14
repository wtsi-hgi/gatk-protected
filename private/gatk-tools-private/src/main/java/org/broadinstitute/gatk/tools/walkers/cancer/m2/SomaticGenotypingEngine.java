/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.cancer.m2;

import com.google.java.contract.Ensures;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.*;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.genotyper.SampleList;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.*;

public class SomaticGenotypingEngine extends HaplotypeCallerGenotypingEngine {

    public static final String TUMOR_LOD = "TLOD";
    public static final String NORMAL_LOD = "NLOD";
    public static final String HAPLOTYPE_COUNT = "HCNT";

    protected M2ArgumentCollection MTAC;

    private final static Logger logger = Logger.getLogger(SomaticGenotypingEngine.class);

    // TODO: understand and remove if possible dependency on HCAC
    public SomaticGenotypingEngine(final HaplotypeCallerArgumentCollection configuration, final SampleList samples, final GenomeLocParser genomeLocParser, final AFCalculatorProvider afCalculatorProvider, final boolean doPhysicalPhasing, final M2ArgumentCollection MTAC) {
        super(configuration, samples, genomeLocParser, afCalculatorProvider, doPhysicalPhasing);
        this.MTAC = MTAC;
    }

    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * The list of samples we're working with is obtained from the readLikelihoods
     *
     * @param haplotypes                             Haplotypes to assign likelihoods to
     * @param readLikelihoods                       Map from reads->(haplotypes,likelihoods)
     * @param perSampleFilteredReadList              Map from sample to reads that were filtered after assembly and before calculating per-read likelihoods.
     * @param ref                                    Reference bytes at active region
     * @param refLoc                                 Corresponding active region genome location
     * @param activeRegionWindow                     Active window
     * @param genomeLocParser                        GenomeLocParser
     * @param activeAllelesToGenotype                Alleles to genotype
     * @param emitReferenceConfidence whether we should add a &lt;NON_REF&gt; alternative allele to the result variation contexts.
     *
     * @return                                       A CalledHaplotypes object containing a list of VC's with genotyped events and called haplotypes
     *
     */
//    @Requires({"refLoc.containsP(activeRegionWindow)", "haplotypes.size() > 0"})
    @Ensures("result != null")
    // TODO - can this be refactored? this is hard to follow!
    public HaplotypeCallerGenotypingEngine.CalledHaplotypes callMutations (
                                                       final List<Haplotype> haplotypes,
                                                       //final Map<String, PerReadAlleleLikelihoodMap> haplotypeReadMap,
                                                       final ReadLikelihoods<Haplotype> readLikelihoods,
                                                       final Map<String, List<GATKSAMRecord>> perSampleFilteredReadList,
                                                       final byte[] ref,
                                                       final GenomeLoc refLoc,
                                                       final GenomeLoc activeRegionWindow,
                                                       final GenomeLocParser genomeLocParser,
                                                       final RefMetaDataTracker tracker,
                                                       final List<VariantContext> activeAllelesToGenotype,
                                                       final boolean emitReferenceConfidence,
                                                       final String tumorSampleName,
                                                       final String matchedNormalSampleName,
                                                       final RodBinding<VariantContext> dbsnpRod,
                                                       final List<RodBinding<VariantContext>> cosmicRod,
                                                       final String DEBUG_READ_NAME

    ) {

        // sanity check input arguments
        if (haplotypes == null || haplotypes.isEmpty()) throw new IllegalArgumentException("haplotypes input should be non-empty and non-null, got "+haplotypes);
        if (readLikelihoods == null || readLikelihoods.sampleCount() == 0) throw new IllegalArgumentException("readLikelihoods input should be non-empty and non-null, got "+readLikelihoods);
        if (ref == null || ref.length == 0 ) throw new IllegalArgumentException("ref bytes input should be non-empty and non-null, got "+ref);
        if (refLoc == null || refLoc.size() != ref.length) throw new IllegalArgumentException(" refLoc must be non-null and length must match ref bytes, got "+refLoc);
        if (activeRegionWindow == null ) throw new IllegalArgumentException("activeRegionWindow must be non-null, got "+activeRegionWindow);
        if (activeAllelesToGenotype == null ) throw new IllegalArgumentException("activeAllelesToGenotype must be non-null, got "+activeAllelesToGenotype);
        if (genomeLocParser == null ) throw new IllegalArgumentException("genomeLocParser must be non-null, got "+genomeLocParser);


        // Somatic Tumor/Normal Sample Handling
        verifySamplePresence(tumorSampleName, readLikelihoods.samples());
        final boolean hasNormal = (matchedNormalSampleName != null);

        // update the haplotypes so we're ready to call, getting the ordered list of positions on the reference
        // that carry events among the haplotypes
        final TreeSet<Integer> startPosKeySet = decomposeHaplotypesIntoVariantContexts(haplotypes, readLikelihoods, ref, refLoc, activeAllelesToGenotype);

        // Walk along each position in the key set and create each event to be outputted
        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();

        for( final int loc : startPosKeySet ) {
            if( loc >= activeRegionWindow.getStart() && loc <= activeRegionWindow.getStop() ) { // genotyping an event inside this active region
                final List<VariantContext> eventsAtThisLoc = getVCsAtThisLocation(haplotypes, loc, activeAllelesToGenotype);

                if( eventsAtThisLoc.isEmpty() ) { continue; }

                // Create the event mapping object which maps the original haplotype events to the events present at just this locus
                final Map<Event, List<Haplotype>> eventMapper = createEventMapper(loc, eventsAtThisLoc, haplotypes);

                // Sanity check the priority list for mistakes
                final List<String> priorityList = makePriorityList(eventsAtThisLoc);

                // Merge the event to find a common reference representation

                VariantContext mergedVC = GATKVariantContextUtils.simpleMerge(eventsAtThisLoc, priorityList,
                        GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                        GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);

                if( mergedVC == null ) { continue; }

//                final VariantContextBuilder vcb = new VariantContextBuilder(mergedVC);

                final GenotypeLikelihoodsCalculationModel.Model calculationModel = mergedVC.isSNP()
                        ? GenotypeLikelihoodsCalculationModel.Model.SNP : GenotypeLikelihoodsCalculationModel.Model.INDEL;

                if (emitReferenceConfidence)
                    mergedVC = addNonRefSymbolicAllele(mergedVC);

                final Map<VariantContext, Allele> mergeMap = new LinkedHashMap<>();
                mergeMap.put(null, mergedVC.getReference()); // the reference event (null) --> the reference allele
                for(int iii = 0; iii < eventsAtThisLoc.size(); iii++) {
                    mergeMap.put(eventsAtThisLoc.get(iii), mergedVC.getAlternateAllele(iii)); // BUGBUG: This is assuming that the order of alleles is the same as the priority list given to simpleMerge function
                }

                final Map<Allele, List<Haplotype>> alleleMapper = createAlleleMapper(mergeMap, eventMapper);

                if( configuration.DEBUG && logger != null ) {
                    if (logger != null) logger.info("Genotyping event at " + loc + " with alleles = " + mergedVC.getAlleles());
                }

                ReadLikelihoods<Allele> readAlleleLikelihoods = readLikelihoods.marginalize(alleleMapper, genomeLocParser.createPaddedGenomeLoc(genomeLocParser.createGenomeLoc(mergedVC), ALLELE_EXTENSION));

                if (!mergedVC.isBiallelic()) {
                    logger.info("[UNSUPPORTED] Detected non-Biallelic VC" + mergedVC.toString());
                    continue;
                }

                // TODO: once tests are passing, refactor to use the new data structure (not the deprecated one)
                // handle overlapping fragments
                // TODO: CONFIRM WITH GSA IF IT IS OK TO REMOVE READS FROM THE PRALM (should be... they do it in filterPoorlyModeledReads!)
                PerReadAlleleLikelihoodMap tumorPRALM = readAlleleLikelihoods.toPerReadAlleleLikelihoodMap(readAlleleLikelihoods.sampleIndex(tumorSampleName));
                filterPRALMForOverlappingReads(tumorPRALM, mergedVC.getReference(), loc, false);
                M2.logReadInfo(DEBUG_READ_NAME, tumorPRALM.getLikelihoodReadMap().keySet(), "Present after filtering for overlapping reads");
                // extend to multiple samples
                double f = estimateAlleleFraction(mergedVC, tumorPRALM);
                double[] tumorGLs = getVariableGenotypeLikelihoods(mergedVC, tumorPRALM, f);

                PerReadAlleleLikelihoodMap normalPRALM = null;
                double[] normalGLs = null;
                if (hasNormal) {
                    normalPRALM = readAlleleLikelihoods.toPerReadAlleleLikelihoodMap(readAlleleLikelihoods.sampleIndex(matchedNormalSampleName));
                    filterPRALMForOverlappingReads(normalPRALM, mergedVC.getReference(), loc, true);
                    M2.logReadInfo(DEBUG_READ_NAME, normalPRALM.getLikelihoodReadMap().keySet(), "Present after filtering for overlapping reads");

                    normalGLs = getVariableGenotypeLikelihoods(mergedVC, normalPRALM, 0.5d);
                }



                double INIT_NORMAL_LOD_THRESHOLD = -Double.MAX_VALUE;
                double NORMAL_LOD_THRESHOLD = -Double.MAX_VALUE;

                int REF = 0, HET = 1;
                double tumorLod = tumorGLs[HET] - tumorGLs[REF];
                double normalLod = 0;
                if (hasNormal) {
                    GenomeLoc eventGenomeLoc = genomeLocParser.createGenomeLoc(activeRegionWindow.getContig(), loc);
                    Collection<VariantContext> cosmicVC = tracker.getValues(cosmicRod, eventGenomeLoc);
                    Collection<VariantContext> dbsnpVC = tracker.getValues(dbsnpRod, eventGenomeLoc);

                    // remove the effect of cosmic from dbSNP
                    boolean germlineAtRisk = (!dbsnpVC.isEmpty() && cosmicVC.isEmpty());

                    // TODO: expose this hardcoded threshold
                    INIT_NORMAL_LOD_THRESHOLD = 0.5;
                    NORMAL_LOD_THRESHOLD = (germlineAtRisk)?MTAC.NORMAL_DBSNP_LOD_THRESHOLD:MTAC.NORMAL_LOD_THRESHOLD;
                    normalLod = normalGLs[REF] - normalGLs[HET];
                }


                VariantContext call = null;
                if (tumorLod >= MTAC.INITIAL_TUMOR_LOD_THRESHOLD && normalLod >= INIT_NORMAL_LOD_THRESHOLD) {
                    VariantContextBuilder callVcb = new VariantContextBuilder(mergedVC);

                    if (normalLod < NORMAL_LOD_THRESHOLD) {
                        callVcb.filter("germline_risk");
                    }

                    // FIXME: can simply get first alternate since above we only deal with Bi-allelic sites...
                    int haplotypeCount = alleleMapper.get(mergedVC.getAlternateAllele(0)).size();
                    callVcb.attribute(HAPLOTYPE_COUNT, haplotypeCount);
                    callVcb.attribute(TUMOR_LOD, tumorLod);
                    callVcb.attribute(NORMAL_LOD, normalLod);

                    if (normalLod < NORMAL_LOD_THRESHOLD) {
                        callVcb.filter("germline_risk");
                    }

                    GenotypeBuilder tumorGenotype =
                            new GenotypeBuilder(tumorSampleName, mergedVC.getAlleles());

                    tumorGenotype.attribute("AF", f);

                    // how should we set the genotype properly here?
                    List<Allele> refAlleles = new ArrayList<>();
                    refAlleles.add(mergedVC.getReference());
                    refAlleles.add(mergedVC.getReference());


                    List<Genotype> genotypes = new ArrayList<>();
                    genotypes.add(tumorGenotype.make());

                    // if we are calling with a normal, add that sample in
                    if (hasNormal) {
                        int[] normalCounts = getRefAltCount(mergedVC, normalPRALM);
                        double normalF = (double) normalCounts[1] / ((double) normalCounts[0] + (double) normalCounts[1]);

                        GenotypeBuilder normalGenotype =
                                new GenotypeBuilder(matchedNormalSampleName, refAlleles).AD(normalCounts);
                        normalGenotype.attribute("AF", normalF);
                        genotypes.add(normalGenotype.make());
                    }

                    call = new VariantContextBuilder(callVcb).genotypes(genotypes).make();

                }

                // how should we be making use of _perSampleFilteredReadList_?
                if( call != null ) {
                    readAlleleLikelihoods = prepareReadAlleleLikelihoodsForAnnotation(readLikelihoods, perSampleFilteredReadList,
                            genomeLocParser, emitReferenceConfidence, alleleMapper, readAlleleLikelihoods, call);

                    VariantContext annotatedCall = annotationEngine.annotateContextForActiveRegion(tracker, readAlleleLikelihoods, call);

                    if( call.getAlleles().size() != mergedVC.getAlleles().size() )
                        annotatedCall = GATKVariantContextUtils.reverseTrimAlleles(annotatedCall);

                    // maintain the set of all called haplotypes
                    for ( final Allele calledAllele : call.getAlleles() ) {
                        final List<Haplotype> haplotypeList = alleleMapper.get(calledAllele);
                        if (haplotypeList == null) continue;
                        calledHaplotypes.addAll(haplotypeList);
                    }

                    returnCalls.add( annotatedCall );
                }

            }
        }

        // TODO: understand effect of enabling this for somatic calling...
        //final List<VariantContext> phasedCalls = doPhysicalPhasing ? phaseCalls(returnCalls, calledHaplotypes) : returnCalls;
        // return new CalledHaplotypes(phasedCalls, calledHaplotypes);
        return new CalledHaplotypes(returnCalls, calledHaplotypes);
    }

    private void verifySamplePresence(String sampleName, List<String> samples) {
        if (!samples.contains(sampleName)) {
            throw new IllegalArgumentException("Unable to find sample name "+sampleName+"in sample list of " + StringUtil.join(",", samples));
        }
    }

    private double[] getVariableGenotypeLikelihoods(VariantContext mergedVC, PerReadAlleleLikelihoodMap tumorPRALM, double f) {
        double[] genotypeLikelihoods = new double[2];
        int AA = 0, AB = 1;
        for(Map.Entry<GATKSAMRecord,Map<Allele, Double>> e : tumorPRALM.getLikelihoodReadMap().entrySet()) {
            Map<Allele, Double> m = e.getValue();
            Double refLL = m.get(mergedVC.getReference());

            // FIXME: what if it's not bi-allelic?  either support it or declare it in the contract
            Double altLL = m.get(mergedVC.getAlternateAllele(0));

            genotypeLikelihoods[AB] += Math.log10(Math.pow(10, refLL) * (1 - f) + Math.pow(10, altLL) * f);
            genotypeLikelihoods[AA] += Math.log10(Math.pow(10, refLL));
        }
        return genotypeLikelihoods;
    }

    // FIXME: calculate using the uncertainty rather than this cheap approach
    // Biallelic check above should ensure only 2 alleles
    private double estimateAlleleFraction(VariantContext vc, PerReadAlleleLikelihoodMap map) {
        int[] counts = getRefAltCount(vc, map);
        int refCount = counts[0];
        int altCount = counts[1];

//        logger.info("Counted " + refCount + " ref and " + altCount + " alt " );
        return (double) altCount / ((double) refCount + (double) altCount);
    }

    // TODO: ensure there are only two alleles in the VC
    private int[] getRefAltCount(VariantContext mergedVC, PerReadAlleleLikelihoodMap afMap) {
        int counts[] = new int[2];
        int REF = 0;
        int ALT = 1;

        for(Map.Entry<GATKSAMRecord,Map<Allele, Double>> e : afMap.getLikelihoodReadMap().entrySet()) {
            Map<Allele, Double> m = e.getValue();
            Double rl = m.get(mergedVC.getReference());
            Double al = m.get(mergedVC.getAlternateAllele(0));

            if (arePairHMMLikelihoodsInformative(rl, al)) {
                if (rl > al) {
                    counts[REF]++;
                } else {
                    counts[ALT]++;
                    logM2Debug("Using " + e.getKey().toString() + " towards alternate allele count");
                }
            }

//            if (al >= rl) logger.info("Alt found in " + e.getKey().getReadName());
        }
        return counts;
    }


    private void logM2Debug(String s) {
        if (MTAC.M2_DEBUG) {
            logger.info(s);
        }
    }

    // would have used org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap.getMostLikelyAllele but we have this case where
    // there is a read that doesn't overlap the variant site, and thus supports both alleles equally.
    private boolean arePairHMMLikelihoodsInformative(double l1, double l2) {
        // TODO: should this be parameterized, or simply encoded
        double EPSILON = 0.1;
        return (Math.abs(l1 - l2) >= EPSILON);
    }

    private void filterPRALMForOverlappingReads(PerReadAlleleLikelihoodMap pralm, Allele ref, int location, boolean retainMismatches) {

        Map<GATKSAMRecord, Map<Allele, Double>> m = pralm.getLikelihoodReadMap();


        // iterate through the reads, if the name has been seen before we have overlapping (potentially) fragments, so handle them
        Map<String, GATKSAMRecord> nameToRead = new HashMap<>();
        Set<GATKSAMRecord> readsToKeep = new HashSet<>();

        for(GATKSAMRecord rec : m.keySet()) {
            // if we haven't seen it... just record the name and add it to the list of reads to keep
            GATKSAMRecord existing = nameToRead.get(rec.getReadName());
            if (existing == null) {
                nameToRead.put(rec.getReadName(), rec);
                readsToKeep.add(rec);
            } else {
                logM2Debug("Found a paired read for " + rec.getReadName());

                // NOTE: Can we use FragmentUtils to do all of this processing (to find overlapping pairs?)
                // seems like maybe, but it has some requirements about the order of the reads supplied which may be painful to meet
                // TODO: CHECK IF THE READS BOTH OVERLAP THE POSITION!!!!
                if ( ReadUtils.isInsideRead(existing, location) && ReadUtils.isInsideRead(rec, location) ) {

                    MostLikelyAllele existingMLA = pralm.getMostLikelyAllele(pralm.getLikelihoodReadMap().get(existing));
                    Allele existingAllele = existingMLA.getMostLikelyAllele();

                    MostLikelyAllele recMLA = pralm.getMostLikelyAllele(pralm.getLikelihoodReadMap().get(rec));
                    Allele recAllele = recMLA.getMostLikelyAllele();

                    // if the reads disagree at this position...
                    if (!existingAllele.equals(recAllele)) {
                        //... and we're not retaining mismatches, throw them both out
                        if (!retainMismatches) {
                            logM2Debug("Discarding read-pair due to disagreement" + rec.getReadName() + " and allele " + existingAllele);
                            readsToKeep.remove(existing);

                            //... and we are retaining mismatches, keep the mismatching one
                        } else {
                            if (existingAllele.equals(ref)) {
                                logM2Debug("Discarding read to keep mismatching " + rec.getReadName() + " and allele " + existingAllele);
                                readsToKeep.remove(existing);
                                readsToKeep.add(rec);
                            }
                        }
                        // Otherwise, keep the element with the higher quality score
                    } else {
                        logM2Debug("Discarding lower quality read of overlapping pair " + rec.getReadName() + " and allele " + existingAllele);
                        if (existingMLA.getLog10LikelihoodOfMostLikely() < recMLA.getLog10LikelihoodOfMostLikely()) {
                            readsToKeep.remove(existing);
                            readsToKeep.add(rec);
                        }
                    }
                } else {
                    // although these are overlapping fragments, they don't overlap at the position in question
                    // so keep the read
                    readsToKeep.add(rec);
                }
            }

        }

        // perhaps moved into PRALM
        final Iterator<Map.Entry<GATKSAMRecord, Map<Allele, Double>>> it = m.entrySet().iterator();
        while ( it.hasNext() ) {
            final Map.Entry<GATKSAMRecord, Map<Allele, Double>> record = it.next();
            if(!readsToKeep.contains(record.getKey())) {
                it.remove();
                logM2Debug("Dropping read " + record.getKey() + " due to overlapping read fragment rules");
            }
        }
    }

    // Move to utility class so we can use one shared with HaplotypeCallerGenotypingEngine
    private VariantContext addNonRefSymbolicAllele(final VariantContext mergedVC) {
        final VariantContextBuilder vcb = new VariantContextBuilder(mergedVC);
        final List<Allele> originalList = mergedVC.getAlleles();
        final List<Allele> alleleList = new ArrayList<>(originalList.size() + 1);
        alleleList.addAll(mergedVC.getAlleles());
        alleleList.add(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        vcb.alleles(alleleList);
        return vcb.make();
    }

}
