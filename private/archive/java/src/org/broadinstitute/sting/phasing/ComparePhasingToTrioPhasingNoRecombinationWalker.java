/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.GenotypePhasingEvaluator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.*;

import java.io.PrintStream;
import java.util.*;

import static org.broadinstitute.sting.utils.variant.GATKVCFUtils.getVCFHeadersFromRods;

/**
 * Walks along all variant ROD loci and compares the phasing between RBP and trio phasing.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})

@ReadFilters({MappingQualityZeroFilter.class})
// Filter out all reads with zero mapping quality

public class ComparePhasingToTrioPhasingNoRecombinationWalker extends RodWalker<CompareResult, CompareToTrioPhasingStats> {
    @Input(fullName="trio", doc="trio VCF", required=true)
    public RodBinding<VariantContext> trio;

    @Input(fullName="phasing", doc="read-backed phasing VCF", required=true)
    public RodBinding<VariantContext> phasing;

    private final static int NUM_IN_TRIO = 3;

    private final static int DIPLOID = 2;

    @Output
    protected PrintStream out;

    @Argument(fullName = "trioAugmentedPhasing", shortName = "trioAugmentedPhasing", doc = "File to which trio-phased variants should be written", required = false)
    protected VCFWriter writer = null;
    
    @Argument(fullName = "diffTrioAndPhasingTracks", shortName = "diffTrioAndPhasingTracks", doc = "File to which comparisons of phasing information in 'trio' and 'phasing' tracks should be written", required = false)
    protected PrintStream diffTrioAndPhasingTracks = null;

    @Argument(fullName = "phasingSample", shortName = "phasingSample", doc = "Name of child sample", required = false)
    protected String phasingSample = null;

    private CompareTrioAndPhasingTracks diffTrioAndPhasingCounts = null;

    private enum TrioStatus {
        PRESENT, MISSING, TRIPLE_HET
    }

    private GenomeLoc prevLoc = null;
    private VariantContext prevTrioVc = null;
    private TrioStatus prevTrioStatus = TrioStatus.MISSING;

    private Genotype prevPhasingGt = null;


    public void initialize() {
        initializeVcfWriter();

        // Will compare the phasing ALREADY present in the trio track [without regards to what this trio phasing mechanism (without recombination) would do]:
        if (diffTrioAndPhasingTracks != null)
            diffTrioAndPhasingCounts = new CompareTrioAndPhasingTracks();
    }

    private void initializeVcfWriter() {
        if (writer == null)
            return;

        // setup the header fields:
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(getToolkit(), Arrays.asList(phasing.getName()));
        Set<String> samples = new TreeSet<String>(rodNameToHeader.get(phasing.getName()).getGenotypeSamples());
        writer.writeHeader(new VCFHeader(hInfo, samples));
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public CompareToTrioPhasingStats reduceInit() {
        return new CompareToTrioPhasingStats();
    }

    /**
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public CompareResult map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        GenomeLoc curLoc = ref.getLocus();
        VariantContext phasingVc = tracker.getFirstValue(phasing, curLoc);

        CompareToTrioPhasingStats stats = new CompareToTrioPhasingStats();
        CompareResult result = new CompareResult(phasingVc, stats);

        if (phasingVc == null || phasingVc.isFiltered())
            return result;

        if (phasingSample == null) {
            GenotypeCollection phasingSampleToGt = phasingVc.getGenotypes();
            if (phasingSampleToGt.size() != 1)
                throw new UserException("Must provide EXACTLY one sample in " + phasing.getName() + " track!");
            phasingSample = phasingSampleToGt.get(0).getSampleName();
        }

        Genotype curPhasingGt = phasingVc.getGenotype(phasingSample);
        if (curPhasingGt == null || !curPhasingGt.isHet()) // can ignore this missing/irrelevant genotype
            return result;

        VariantContext curTrioVc = tracker.getFirstValue(trio, curLoc);
        boolean useTrioVc = (curTrioVc != null && !curTrioVc.isFiltered());

        Genotype sampleCurGtInTrio = null;
        if (useTrioVc) {
            sampleCurGtInTrio = curTrioVc.getGenotype(phasingSample);

            if (curTrioVc.getNSamples() > NUM_IN_TRIO || sampleCurGtInTrio == null)
                throw new UserException("Must provide " + trio.getName() + " track data for sample: " + phasingSample);

            if (!curPhasingGt.sameGenotype(sampleCurGtInTrio)) {
                logger.warn("Locus " + curLoc + " breaks phase, since " + phasing.getName() + " and " + trio.getName() + " tracks have different genotypes for " + phasingSample + "!");
                prevLoc = null;
                return result;
            }
        }

        // Now, we have a [trio-consistent] het genotype that may be phased or not [and we want to know if it could be phased based on trio information]:
        TrioStatus currentTrioStatus = TrioStatus.MISSING;
        if (useTrioVc)
            currentTrioStatus = determineTrioStatus(curTrioVc);

        if (prevLoc != null && curLoc.onSameContig(prevLoc)) {
            String trioPhaseStatus;
            stats.comparedSites++;
            String addToOutput = "";

            if (prevTrioStatus == TrioStatus.TRIPLE_HET || currentTrioStatus == TrioStatus.TRIPLE_HET) {
                trioPhaseStatus = "Het3";
            }
            else if (prevTrioStatus == TrioStatus.MISSING || currentTrioStatus == TrioStatus.MISSING) {
                trioPhaseStatus = "Missing";
            }
            else {
                if (prevTrioStatus != TrioStatus.PRESENT || currentTrioStatus != TrioStatus.PRESENT)
                    throw new ReviewedStingException("LOGICAL error: prevTrioStatus != TrioStatus.PRESENT || currentTrioStatus != TrioStatus.PRESENT");

                trioPhaseStatus = "trio_phased";
                stats.trioPhaseableSites++;

                if (writer != null) { // Phase the genotypes using the trio information:
                    String parent1 = null;
                    String parent2 = null;
                    for (final Genotype trio : curTrioVc.getGenotypes()) {
                        String trioSample = trio.getSampleName();
                        if (trio.getPloidy() != DIPLOID)
                            throw new UserException("Each sample in trio must be diploid!");
                        if (trioSample.equals(phasingSample))
                            continue;

                        if (parent1 == null)
                            parent1 = trioSample;
                        else if (parent2 == null)
                            parent2 = trioSample;
                        else
                            throw new ReviewedStingException("Cannot be more than 2 parents in TRIO!");
                    }
                    if (parent1 == null || parent2 == null)
                        throw new ReviewedStingException("Must have 2 parents in TRIO!");

                    Genotype samplePrevGtInTrio = prevTrioVc.getGenotype(phasingSample);

                    Genotype parent1PrevGt = prevTrioVc.getGenotype(parent1);
                    Genotype parent1CurGt = curTrioVc.getGenotype(parent1);

                    Genotype parent2PrevGt = prevTrioVc.getGenotype(parent2);
                    Genotype parent2CurGt = curTrioVc.getGenotype(parent2);

                    int prevHomIndex, prevOtherIndex;
                    Allele prevHomAllele;
                    Set<Allele> prevOtherAlleles;
                    if (parent1PrevGt.isHom()) {
                        prevHomIndex = 1;
                        prevOtherIndex = 2;
                        prevHomAllele = parent1PrevGt.getAllele(0);
                        prevOtherAlleles = new TreeSet<Allele>(parent2PrevGt.getAlleles());
                    }
                    else if (parent2PrevGt.isHom()) {
                        prevHomIndex = 2;
                        prevOtherIndex = 1;
                        prevHomAllele = parent2PrevGt.getAllele(0);
                        prevOtherAlleles = new TreeSet<Allele>(parent1PrevGt.getAlleles());
                    }
                    else
                        throw new ReviewedStingException("LOGICAL ERROR: at least one parent is hom!");

                    int curHomIndex, curOtherIndex;
                    Allele curHomAllele;
                    Set<Allele> curOtherAlleles;
                    if (parent1CurGt.isHom()) {
                        curHomIndex = 1;
                        curOtherIndex = 2;
                        curHomAllele = parent1CurGt.getAllele(0);
                        curOtherAlleles = new TreeSet<Allele>(parent2CurGt.getAlleles());
                    }
                    else if (parent2CurGt.isHom()) {
                        curHomIndex = 2;
                        curOtherIndex = 1;
                        curHomAllele = parent2CurGt.getAllele(0);
                        curOtherAlleles = new TreeSet<Allele>(parent1CurGt.getAlleles());
                    }
                    else
                        throw new ReviewedStingException("LOGICAL ERROR: at least one parent is hom!");

                    boolean phased = true;

                    Map<Allele, Integer> prevAlleleToParent = new TreeMap<Allele, Integer>();
                    for (Allele prevAllele : samplePrevGtInTrio.getAlleles()) {
                        if (prevAllele.equals(prevHomAllele))
                            prevAlleleToParent.put(prevAllele, prevHomIndex);
                        else if (prevOtherAlleles.contains(prevAllele))
                            prevAlleleToParent.put(prevAllele, prevOtherIndex);
                        else {
                            logger.warn("CANNOT trio phase, due to inconsistent inheritance of alleles at: " + GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), prevTrioVc));
                            phased = false;
                            break;
                        }
                    }

                    Map<Integer, Allele> parentToCurAllele = new HashMap<Integer, Allele>();
                    for (Allele curAllele : sampleCurGtInTrio.getAlleles()) {
                        if (curAllele.equals(curHomAllele))
                            parentToCurAllele.put(curHomIndex, curAllele);
                        else if (curOtherAlleles.contains(curAllele))
                            parentToCurAllele.put(curOtherIndex, curAllele);
                        else {
                            logger.warn("CANNOT trio phase, due to inconsistent inheritance of alleles at: " + curLoc);
                            phased = false;
                            break;
                        }
                    }

                    if (phased) {
                        List<Allele> phasedCurAlleles = new LinkedList<Allele>();
                        for (Allele prevAllele : prevPhasingGt.getAlleles()) {
                            Integer prevIndex = prevAlleleToParent.get(prevAllele);
                            if (prevIndex == null)
                                throw new ReviewedStingException("LOGICAL error: expecting to find prev allele in trio parents");
                            Allele curAllele = parentToCurAllele.get(prevIndex);
                            if (curAllele == null)
                                throw new ReviewedStingException("LOGICAL error: expecting to find cur allele in trio parents");
                            phasedCurAlleles.add(curAllele);
                        }

                        boolean useTrioPhase = true;
                        Genotype phasedGt = new Genotype(phasingSample, phasedCurAlleles, curPhasingGt.getNegLog10PError(), curPhasingGt.getFilters(), curPhasingGt.getAttributes(), phased);

                        if (curPhasingGt.isPhased()) {
                            stats.bothCanPhase++;
                            useTrioPhase = false;

                            boolean ignorePhase = false;
                            if (!phasedGt.sameGenotype(curPhasingGt, ignorePhase)) {
                                String contradictMessage = "Phase from " + phasing.getName() + " track at " + curLoc + " contradicts the trio-based phasing.";
                                stats.contradictoryPhaseSites++;
                                addToOutput += "\tcontradictory";

                                if (phasingVc.hasAttribute(ReadBackedPhasing.PHASING_INCONSISTENT_KEY)) {
                                    stats.contradictoryPhaseSitesWithPhaseInconsistency++;
                                    addToOutput += "\tphaseInconsistent";                                    
                                    useTrioPhase = true;
                                    contradictMessage += " Ignoring " + phasing.getName() + " phase due to phase-inconsistency.";
                                }
                                else {
                                    contradictMessage += " Maintaining phase from " + phasing.getName() + ".";
                                }
                                logger.warn(contradictMessage);
                            }
                        }

                        if (useTrioPhase) { // trio phasing adds PREVIOUSLY UNKNOWN phase information:
                            GenotypeCollection genotypes = phasingVc.getGenotypes();
                            genotypes.add(phasedGt);

                            phasingVc = VariantContext.modifyGenotypes(phasingVc, genotypes);
                            result.phasedVc = phasingVc;
                        }
                    }
                }
            }
            out.println(prevLoc + "\t" + curLoc + "\t" + trioPhaseStatus + "\t" + curPhasingGt.isPhased() + addToOutput);
            
            if (diffTrioAndPhasingTracks != null && prevTrioStatus != TrioStatus.MISSING && currentTrioStatus != TrioStatus.MISSING && sampleCurGtInTrio.isPhased() && curPhasingGt.isPhased()) {
                AllelePair prevTrioAll = new AllelePair(prevTrioVc.getGenotype(phasingSample));
                AllelePair curTrioAll = new AllelePair(sampleCurGtInTrio);
                
                AllelePair prevPhasingAll = new AllelePair(prevPhasingGt);
                AllelePair curPhasingAll = new AllelePair(curPhasingGt);
                
                boolean topsMatch = (GenotypePhasingEvaluator.topMatchesTop(prevTrioAll, prevPhasingAll) && GenotypePhasingEvaluator.topMatchesTop(curTrioAll, curPhasingAll));
                boolean bottomsMatch = (GenotypePhasingEvaluator.bottomMatchesBottom(prevTrioAll, prevPhasingAll) && GenotypePhasingEvaluator.bottomMatchesBottom(curTrioAll, curPhasingAll));

                boolean topMatchesBottom = (GenotypePhasingEvaluator.topMatchesBottom(prevTrioAll, prevPhasingAll) && GenotypePhasingEvaluator.topMatchesBottom(curTrioAll, curPhasingAll));
                boolean bottomMatchesTop = (GenotypePhasingEvaluator.bottomMatchesTop(prevTrioAll, prevPhasingAll) && GenotypePhasingEvaluator.bottomMatchesTop(curTrioAll, curPhasingAll));

                boolean phasesAgree = ((topsMatch && bottomsMatch) || (topMatchesBottom && bottomMatchesTop));

                diffTrioAndPhasingTracks.println(prevLoc + "\t" + curLoc + "\t" + trioPhaseStatus + "\t" + phasesAgree);
                diffTrioAndPhasingCounts.addComparison(trioPhaseStatus, phasesAgree);
            }
        }

        prevLoc = curLoc;
        prevTrioVc = curTrioVc;
        prevTrioStatus = currentTrioStatus;
        prevPhasingGt = curPhasingGt;

        return result;
    }

    private static TrioStatus determineTrioStatus(VariantContext trioVc) {
        if (trioVc.getNSamples() != NUM_IN_TRIO)
            return TrioStatus.MISSING;

        for (int i = 0; i < NUM_IN_TRIO; i++) {
            Genotype gtI = trioVc.getGenotype(i);
            if (gtI.isNoCall() || gtI.isFiltered())
                return TrioStatus.MISSING;

            if (!gtI.isHet())
                return TrioStatus.PRESENT;
        }

        return TrioStatus.TRIPLE_HET;
    }

    public CompareToTrioPhasingStats reduce(CompareResult addIn, CompareToTrioPhasingStats runningCount) {
        if (addIn == null)
            addIn = new CompareResult();

        if (writer != null && addIn.phasedVc != null)
            WriteVCF.writeVCF(addIn.phasedVc, writer, logger);

        return runningCount.addIn(addIn.stats);
    }

    /**
     * @param result the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(CompareToTrioPhasingStats result) {
        System.out.println("Compared " + result.comparedSites + " sites.");
        System.out.println("Trio can phase " + result.trioPhaseableSites + " sites.");
        System.out.println("Trio and " + phasing.getName() + " track can both phase " + result.bothCanPhase + " sites.");
        System.out.println("Contradiction between phase inferred from " + trio.getName() + " and phase present in " + phasing.getName() + " tracks at " + result.contradictoryPhaseSites + " sites.");
        System.out.println("Of those, " + phasing.getName() + " track is phase-inconsistent at " + result.contradictoryPhaseSitesWithPhaseInconsistency + " sites.");

        if (diffTrioAndPhasingCounts != null) {
            System.out.println("");
            diffTrioAndPhasingCounts.printSummary(System.out);
        }
    }
}

class CompareToTrioPhasingStats {
    public int comparedSites;
    public int trioPhaseableSites;
    public int contradictoryPhaseSites;
    public int contradictoryPhaseSitesWithPhaseInconsistency;
    public int bothCanPhase;

    public CompareToTrioPhasingStats() {
        this.comparedSites = 0;
        this.trioPhaseableSites = 0;
        this.contradictoryPhaseSites = 0;
        this.contradictoryPhaseSitesWithPhaseInconsistency = 0;
        this.bothCanPhase = 0;
    }

    public CompareToTrioPhasingStats addIn(CompareToTrioPhasingStats other) {
        this.comparedSites += other.comparedSites;
        this.trioPhaseableSites += other.trioPhaseableSites;
        this.contradictoryPhaseSites += other.contradictoryPhaseSites;
        this.contradictoryPhaseSitesWithPhaseInconsistency += other.contradictoryPhaseSitesWithPhaseInconsistency;
        this.bothCanPhase += other.bothCanPhase;

        return this;
    }
}

class CompareResult {
    public VariantContext phasedVc;
    public CompareToTrioPhasingStats stats;

    public CompareResult() {
        this.phasedVc = null;
        this.stats = new CompareToTrioPhasingStats();
    }

    public CompareResult(VariantContext phasedVc, CompareToTrioPhasingStats stats) {
        this.phasedVc = phasedVc;
        this.stats = stats;
    }
}

class CompareTrioAndPhasingTracks {
    private Map<String, AgreeDisagreeCounts> trioStatusToAgreement;

    public CompareTrioAndPhasingTracks() {
        this.trioStatusToAgreement = new HashMap<String, AgreeDisagreeCounts>();
    }

    public void addComparison(String trioStatus, boolean agree) {
        AgreeDisagreeCounts counts = trioStatusToAgreement.get(trioStatus);
        if (counts == null) {
            counts = new AgreeDisagreeCounts();
            trioStatusToAgreement.put(trioStatus, counts);
        }

        if (agree)
            counts.incrementAgree();
        else
            counts.incrementDisagree();
    }

    public void printSummary(PrintStream out) {
        out.println("--------------------------------------------");
        out.println("Summary of trio vs. phasing tracks' phasing:");
        out.println("--------------------------------------------");        

        int globalAgree = 0;
        int globalDisagree = 0;
        for (AgreeDisagreeCounts counts : trioStatusToAgreement.values()) {
            globalAgree += counts.agree;
            globalDisagree += counts.disagree;
        }
        int globalTotal = globalAgree + globalDisagree;

        out.println("Concordant phase:\t" + percentString(globalAgree, globalTotal));
        out.println("Discordant phase:\t" + percentString(globalDisagree, globalTotal));

        for (Map.Entry<String, AgreeDisagreeCounts> statusCounts : trioStatusToAgreement.entrySet()) {
            String status = statusCounts.getKey();
            AgreeDisagreeCounts counts = statusCounts.getValue();

            out.println("");
            out.println("'" + status + "'" + " Concordant phase:\t" + percentString(counts.agree, counts.total()));
            out.println("'" + status + "'" + " Discordant phase:\t" + percentString(counts.disagree, counts.total()));
        }
        out.println("--------------------------------------------");
        out.println("");
    }

    private static String percentString(int numerator, int denominator) {
        int NUM_DECIMAL_PLACES = 1;
        String percent = new Formatter().format("%." + NUM_DECIMAL_PLACES + "f", MathUtils.percentage(numerator, denominator)).toString();

        StringBuilder sb = new StringBuilder();
        sb.append(numerator).append(" (").append(percent).append("%)");

        return sb.toString();
    }
}

class AgreeDisagreeCounts {
    protected int agree;
    protected int disagree;

    public AgreeDisagreeCounts() {
        this.agree = 0;
        this.disagree = 0;
    }

    public void incrementAgree() {
        agree++;
    }

    public void incrementDisagree() {
        disagree++;
    }

    public int total() {
        return agree + disagree;
    }
}