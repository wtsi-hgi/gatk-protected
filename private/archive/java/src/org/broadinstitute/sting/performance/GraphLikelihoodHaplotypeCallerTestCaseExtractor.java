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
package org.broadinstitute.sting.gatk.walkers.haplotypecaller;


import com.google.java.contract.Requires;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.na12878kb.assess.AssessNA12878;
import org.broadinstitute.sting.gatk.walkers.na12878kb.assess.AssessmentType;
import org.broadinstitute.sting.gatk.walkers.na12878kb.assess.Assessor;
import org.broadinstitute.sting.gatk.walkers.na12878kb.assess.BadSitesWriter;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.*;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.InvalidRecordHandler;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.GATKSamRecordFactory;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.*;

import java.io.File;
import java.util.*;

/**
 * Walker to generate test case with FP FN errors between the pair-hmm and the graph-based likelihoods. Not meant to
 * be outside testing (no sense whatsoever to do that).
 */
public class GraphLikelihoodHaplotypeCallerTestCaseExtractor extends HaplotypeCaller {

    @Argument(fullName="activeReadsOutput", shortName="readsout",doc="Write output to this BAM filename instead of STDOUT", required = false)
    protected File activeReadsOutputFile;

    protected SAMFileWriter activeReadsOutput;

    @Argument(fullName="activeRegionsOutput", shortName="arout", doc="Write output with the information about active regions that present problems", required = false)
    protected File activeRegionsOutputFile;

    protected VariantContextWriter activeRegionsOutput;

    @Argument(fullName="activeRegionsInput", shortName="arin", doc="VCF describe the test active regions",required=false)
    protected List<RodBinding<VariantContext>> activeRegionsInput;

    @Argument(shortName="dontDiscardLoglessErrors", fullName="dontDiscardLoglessErrors",
            doc="whether we only should output data on errors by the graph likelihood based method that are not produced using the regular pair-hmm method.",required= false)
    public boolean dontDiscardLoglessErrors = false;


    @Argument(shortName="mode", doc="mode of execution", required = true)
    public Mode mode;

    public enum Mode {
        GENERATION, TESTING
    }

    @ArgumentCollection
    private NA12878DBArgumentCollection kbArgs = new NA12878DBArgumentCollection();

    protected NA12878KnowledgeBase kbDb;
    protected SiteIterator<MongoVariantContext> consensusSiteIterator;


    public static final String ACTIVE_REGION_ERRORS_KEY = "Errors";
    public static final String ACTIVE_REGION_START_KEY = "Start";
    public static final String ACTIVE_REGION_END_KEY = "End";
    public static final String ACTIVE_REGION_EXTENDED_START_KEY = "ExtStart";
    public static final String ACTIVE_REGION_EXTENDED_END_KEY = "ExtEnd";
    public static final String ACTIVE_REGION_ERROR_CLASSES_KEY = "ErrorClasses";

    @Argument(fullName="excludeCallset", shortName = "excludeCallset", doc="Don't count calls that come from only these excluded callsets", required=false)
    public Set<String> excludeCallset = Collections.singleton("CEUTrio_best_practices");

    /**
     * An output VCF file containing the bad sites (FN/FP) that were found in the input callset w.r.t. the current NA12878 knowledge base
     */
    @Argument(fullName = "badSites", shortName = "badSites", doc="VCF file containing information on FP/FNs in the input callset", required=false)
    protected File badSitesFile;
    protected VariantContextWriter badSites = null;

    @Argument(fullName = "badSitesInput", shortName = "badSitesInput", doc="input badsites", required=false)
    protected List<RodBinding<VariantContext>> badSitesInput;

    @Argument(fullName="maxToWrite", shortName = "maxToWrite", doc="Max. number of bad sites to write out", required=false)
    public int maxToWrite = 10000;

    @Argument(fullName="minDepthForLowCoverage", shortName = "minDepthForLowCoverage", doc="A false negative will be flagged as due to low coverage if the (optional) BAM is provided and the coverage overlapping the site is less than this value", required=false)
    public int minDepthForLowCoverage = 5;

    @Argument(fullName="requireReviewed", shortName = "requireReviewed", doc="If true, we will only use reviewed sites for the analysis", required=false)
    public boolean onlyReviewed = false;

    @Argument(fullName="typesToInclude", shortName = "typesToInclude", doc="Should we analyze SNPs, INDELs, or both?", required=false)
    public AssessNA12878.TypesToInclude typesToInclude = AssessNA12878.TypesToInclude.BOTH;

    @Argument(fullName="AssessmentsToExclude", shortName = "AssessmentsToExclude", doc="If provided, we will prevent any of these states from being written out to the badSites VCF.", required=false)
    public Set<AssessmentType> AssessmentsToExclude = EnumSet.noneOf(AssessmentType.class);

    private SAMFileReader bamReader = null;
    private MyBadSitesWriter badSitesWriter;



    private Assessor assessor;



    @Override
    public void initialize() {
        super.initialize();
        if (mode == Mode.TESTING)
            testingModeInitialize();
        else
            generationModeInitialize();
    }


    private List<TestCaseData> testDataList;

    private void testingModeInitialize()  {
        if (activeRegionsInput == null || activeRegionsInput.size() == 0)
            throw new StingException("in testing mode you need to provide an active regions input vcf");
        if (badSitesInput == null || badSitesInput.size() == 0)
            throw new StingException("in testing mode you need to provide an active bad site input vcf");

        badSitesOutputInitialize();
        initializeKBConnection();
        try {
            testDataList = GraphLikelihoodVsLoglessAccuracyIntegrationTest.getTestCaseDataList();
        } catch (Exception e) {
            logger.warn("Cannot resolve the testing class and test-case-data list");
            testDataList = new LinkedList<>();
            throw new RuntimeException(e);
        }
    }


    private void badSitesOutputInitialize() {
        badSites = VariantContextWriterFactory.sortOnTheFly(VariantContextWriterFactory.create(badSitesFile, getToolkit().getSAMFileHeader().getSequenceDictionary(), EnumSet.of(Options.INDEX_ON_THE_FLY,
                Options.ALLOW_MISSING_FIELDS_IN_HEADER)), 10000);
        badSitesWriter = new MyBadSitesWriter(maxToWrite, AssessmentsToExclude, badSites);
        badSitesWriter.initialize(GATKVCFUtils.getHeaderFields(getToolkit()));
    }


    public void generationModeInitialize() {
        if (activeReadsOutputFile == null) throw new StingException("in generation mode the active-reads-output must be specified");
        if (activeRegionsOutputFile == null) throw new StingException("in generation mode the active-region-output must be specified");
        if (badSitesFile == null) throw new StingException("in generation mode the bad-sites-output must be specifed");
        activeReadsOutput =  new SAMFileWriterFactory().makeBAMWriter(getToolkit().getSAMFileHeader(),false,activeReadsOutputFile);

        badSitesOutputInitialize();
        initializeKBConnection();

        activeRegionsOutput = VariantContextWriterFactory.create(activeRegionsOutputFile,getToolkit().getSAMFileHeader().getSequenceDictionary(), EnumSet.of(Options.INDEX_ON_THE_FLY,
                Options.DO_NOT_WRITE_GENOTYPES,
                Options.ALLOW_MISSING_FIELDS_IN_HEADER));
        activeRegionsOutput.writeHeader(activeRegionsHeader());
        if (getToolkit().getSAMFileHeaders().size() != 1)
            throw new UserException.BadInput("this walker need and supports only one input BAM file");
        bamReader = new SAMFileReader(getToolkit().getReadsDataSource().getSAMFile(getToolkit().getReadsDataSource().getReaderIDs().iterator().next()));
        bamReader.setSAMRecordFactory(new GATKSamRecordFactory());
    }

    private void initializeKBConnection() {
        if (kbArgs.dbToUse == NA12878DBArgumentCollection.DBType.DEFAULT)
            kbArgs.dbToUse = NA12878DBArgumentCollection.DBType.PRODUCTION;
        kbDb = new NA12878KnowledgeBase(getToolkit().getGenomeLocParser(), kbArgs);

        assessor = new Assessor("NA12878",typesToInclude, excludeCallset, badSitesWriter, bamReader, minDepthForLowCoverage, -1,0,false);
        consensusSiteIterator = new MySiteIterator(kbDb);
    }


    protected class MyBadSitesWriter extends BadSitesWriter {

        protected ActiveRegion currentActiveRegion;

        public void setCurrentActiveRegion(final ActiveRegion ar) {
            currentActiveRegion = ar;
        }

        @Override
        public void initialize(final Set<VCFHeaderLine> lines) {

            final Set<VCFHeaderLine> augmentedLines = new HashSet<>(lines);
            augmentedLines.add(new VCFInfoHeaderLine("AR",1,VCFHeaderLineType.Character,"active region where the bad-site occurs"));
            super.initialize(augmentedLines);
        }

        public void resetActiveRegion() {
            currentActiveRegion = null;
        }

        protected boolean outputBadSites = true;
        public void switchOffOutput() {
            outputBadSites = false;
        }
        public void switchOnOutput() {
            outputBadSites = true;

        }

        protected List<VariantContext> badSiteBuffer = new LinkedList<>();

        /**
         * Create a new BadSitesWriter
         *
         * @param maxToWrite           the maximum number of records to write
         * @param assessmentsToExclude don't include assessments in this set, even if they would normally be emitted
         * @param writer               the underlying VCWriter we'll use to emit bad sites.  Can be null if captureBadSites is false
         */
        public MyBadSitesWriter(final int maxToWrite, final Set<AssessmentType> assessmentsToExclude, final VariantContextWriter writer) {
            super(maxToWrite, assessmentsToExclude, writer);
        }

        @Override
        protected void outputSite(final VariantContext vc) {
            if (outputBadSites) {
                if (currentActiveRegion != null) {
                    final VariantContextBuilder vcb = new VariantContextBuilder(vc);
                    vcb.attribute("AR",currentActiveRegion.getLocation().getContig() + "_" + currentActiveRegion.getLocation().getStart());
                    super.outputSite(vcb.make());
                } else {
                    super.outputSite(vc);
                }
            }
            final AssessmentType at = AssessmentType.valueOf(vc.getAttributeAsString("WHY","NOT_RELEVANT"));
            switch (at) {
                case TRUE_POSITIVE:
                case TRUE_NEGATIVE:
                case CORRECTLY_FILTERED:
                case CORRECTLY_UNCALLED:
                case CALLED_IN_DB_UNKNOWN_STATUS:
                case CALLED_NOT_IN_DB_AT_ALL:
                case NOT_RELEVANT:
                case REASONABLE_FILTERS_WOULD_FILTER_FP_SITE:
                    break;
                default:
                    badSiteBuffer.add(vc);
            }
        }

        public void outputBadSites(final List<VariantContext> vcs) {
            Collections.sort(vcs,vcfSorter);
            for (final VariantContext vc : vcs)
                outputSite(vc);
        }

        protected ArrayList<VariantContext> flushBadSiteBuffer() {
            final ArrayList<VariantContext> result = new ArrayList<>(badSiteBuffer);
            badSiteBuffer.clear();
            return result;
        }

    }

    private VCFHeader activeRegionsHeader() {
        final VCFHeader result = new VCFHeader();
        result.addMetaDataLine(new VCFInfoHeaderLine(ACTIVE_REGION_ERROR_CLASSES_KEY, VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.Character,"classes of lk or genotyping errors encounterd in that active region"));
        result.addMetaDataLine(new VCFInfoHeaderLine(ACTIVE_REGION_ERRORS_KEY, VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.Character, "errors found in the active region as compared with knowledge based and optionally alternative HC method."));
        result.addMetaDataLine(new VCFInfoHeaderLine(ACTIVE_REGION_START_KEY, 1, VCFHeaderLineType.Integer,"Active region start"));
        result.addMetaDataLine(new VCFInfoHeaderLine(ACTIVE_REGION_END_KEY, 1, VCFHeaderLineType.Integer,"Active region end"));
        result.addMetaDataLine(new VCFInfoHeaderLine(ACTIVE_REGION_EXTENDED_START_KEY, 1, VCFHeaderLineType.Integer, "Active region extended start"));
        result.addMetaDataLine(new VCFInfoHeaderLine(ACTIVE_REGION_EXTENDED_END_KEY, 1, VCFHeaderLineType.Integer,"Active region extended end"));
        return result;
    }

    private static final List<VariantContext> NO_CALLS = Collections.emptyList();


    private SortedSet<GATKSAMRecord> outputReads = new TreeSet<GATKSAMRecord>(new Comparator<GATKSAMRecord>() {
        @Override
        public int compare(final GATKSAMRecord o1, final GATKSAMRecord o2) {
            if (o1 == o2) return 0;
            final int r1 = o1.getReferenceIndex();
            final int r2 = o2.getReferenceIndex();
            if (r1 == r2) {
                  final int e1 = o1.getAlignmentEnd();
                  final int e2 = o2.getAlignmentEnd();
                  if (e1 == e2) {
                      return o1.getReadName().compareTo(o2.getReadName());
                  }
                  else return e1 < e2 ? 1 : -1;
            } else {
                return r1 < r2 ? 1 : -1;
            }
        }
    });


    /**
     * VCF Record sorting comparator
     */
    private Comparator<VariantContext> vcfSorter = new Comparator<VariantContext>() {

        @Override
        public int compare(final VariantContext o1, final VariantContext o2) {
            final int r1 = getToolkit().getGenomeLocParser().getContigIndex(o1.getChr());
            final int r2 = getToolkit().getGenomeLocParser().getContigIndex(o2.getChr());
            if (r1 == r2) {
                final int s1 = o1.getStart();
                final int s2 = o2.getStart();
                if (s1 == s2) return 0;
                else if (s1 <= 0) return 1;
                else if (s2 <= 0) return -1;
                else return s1 < s2 ? -1 : 1;

            } else {
                return r1 < r2 ? -1 : 1;
            }
        }
    };


    /**
     * Read sorter comparator
     */
    private static Comparator<GATKSAMRecord> readsSorter = new Comparator<GATKSAMRecord>() {

        @Override
        public int compare(final GATKSAMRecord o1, final GATKSAMRecord o2) {
            final int r1 = o1.getReferenceIndex();
            final int r2 = o2.getReferenceIndex();
            if (r1 == r2) {
                final int s1 = o1.getAlignmentStart();
                final int s2 = o2.getAlignmentStart();
                if (s1 == s2) return o1.getReadName().compareTo(o2.getReadName());
                else if (s1 <= 0) return 1;
                else if (s2 <= 0) return -1;
                else return s1 < s2 ? -1 : 1;
            } else {
                return r1 < r2 ? -1 : 1;
            }
        }
    };

    /**
     * Map method for this walker.
     *
     * <p>
     *     It delegates on {@link #mapTesting} and {@link #mapGeneration} to implement each mode modes.
     * </p>
     * @param originalActiveRegion
     * @param metaDataTracker
     * @return never null.
     */
    @Override
    public List<VariantContext> map(final ActiveRegion originalActiveRegion, final RefMetaDataTracker metaDataTracker) {
        if (activeRegionsOutputFile == null) {
            return mapTesting(originalActiveRegion, metaDataTracker);
        } else {
            return mapGeneration(originalActiveRegion, metaDataTracker);
        }
    }


    /**
     * Implements the map function for the TESTING mode.
     *
     * TODO there quite a bit of code shared with mapGeneration so we chould try to dry it a bit more. Yet this is
     * TODO just code for testing purposes thus not priority.
     */
    private List<VariantContext> mapTesting(final ActiveRegion originalActiveRegion, final RefMetaDataTracker metaDataTracker) {
        final List<GATKSAMRecord> reads = copyReads(originalActiveRegion);

        // Ignore sites between active regions as that is in any case a problem with the active region traversal.
        final List<MongoVariantContext> missedSites = consensusSiteIterator.getSitesBefore(originalActiveRegion.getLocation().getStartLocation());


        // We remove output reads from the buffer.
        removePassedReadsFromOutputBuffer(originalActiveRegion);
        badSitesWriter.setCurrentActiveRegion(originalActiveRegion);

        final VariantContext generatedActiveRegion = findGeneratedActiveRegion(originalActiveRegion, metaDataTracker);

        final List<VariantContext> generatedBadSites = generatedActiveRegion == null ? Collections.EMPTY_LIST : metaDataTracker.getValues(badSitesInput);

        if (generatedBadSites.size() == 0 && generatedActiveRegion != null)
            throw new StingException("perhaps you used different data, as there is no reported bad-sites for a repoted active region " + originalActiveRegion.getLocation());

        //Second we run the method on this active region using graph-likelihoods and loglessPairHMM
        final boolean originalUseGraphLikelihoods = likelihoodEngineImplementation == LikelihoodCalculationEngine.Implementation.GraphBased;
        likelihoodEngineImplementation = LikelihoodCalculationEngine.Implementation.GraphBased;
        final List<VariantContext> graphResults = super.map(originalActiveRegion, metaDataTracker);
        final List<VariantContext> loglessResults = new LinkedList<>();
        if (!originalUseGraphLikelihoods || !dontDiscardLoglessErrors) {
            likelihoodEngineImplementation = LikelihoodCalculationEngine.Implementation.PairHMM;
            loglessResults.addAll(super.map(originalActiveRegion, metaDataTracker));
        }

        final List<MongoVariantContext> assessmentSites = consensusSiteIterator.getSitesBefore(originalActiveRegion.getLocation().getStopLocation().incPos());
        final List<VariantContext> relevantBadSites;

        // We prevent the assessor to output badSites:
        badSitesWriter.switchOffOutput();
        if (dontDiscardLoglessErrors) {
            assessor.assessSite(graphResults, assessmentSites, onlyReviewed);
            relevantBadSites = badSitesWriter.flushBadSiteBuffer();
        } else {
            assessor.assessSite(loglessResults, assessmentSites, onlyReviewed);
            final List<VariantContext> loglessBadSites = badSitesWriter.flushBadSiteBuffer();
            assessor.assessSite(graphResults,assessmentSites,onlyReviewed);
            final List<VariantContext> graphBadSites = badSitesWriter.flushBadSiteBuffer();
            relevantBadSites = badSiteSetDiff(graphBadSites,loglessBadSites);
        }
        // Now we allow the output of badSites (those in relevantBadSites).
        badSitesWriter.switchOnOutput();
        //reportBadSitesAndActiveRegions(originalActiveRegion, relevantBadSites, reads);
        final List<Pair<VariantContext,VariantContext>> pairs = badSiteMatch(relevantBadSites, generatedBadSites);
        if (pairs.size() != 0) testDataList.add(new TestCaseData(originalActiveRegion,generatedActiveRegion,pairs));
        return originalUseGraphLikelihoods ? graphResults : loglessResults;
    }

    // Deep-clone the active-region read list.
    private List<GATKSAMRecord> copyReads(final ActiveRegion ar) {
       final List<GATKSAMRecord> result = new ArrayList<>(ar.getReads().size());
       for (final GATKSAMRecord read : ar.getReads())
           try {
               result.add((GATKSAMRecord)read.clone());
           } catch (final CloneNotSupportedException e) {
               throw new StingException(e.getMessage(),e);
           }
        return result;
    }

    /**
     * Search for a matching active region in the input active region ROD
     * @param originalActiveRegion target active region
     * @param metaDataTracker
     * @return null if none was found.
     */
    private VariantContext findGeneratedActiveRegion(final ActiveRegion originalActiveRegion, final RefMetaDataTracker metaDataTracker) {
        final List<VariantContext> generatedActiveRegions = metaDataTracker.getValues(activeRegionsInput,originalActiveRegion.getLocation().getStartLocation());
        VariantContext generatedActiveRegion;
        if (generatedActiveRegions.size() == 0)
            generatedActiveRegion = null;
        else {
            generatedActiveRegion = null;
            for (final VariantContext candidate : generatedActiveRegions) {
                final int candidateStart = candidate.getStart();
                final int candidateStop = candidate.getAttributeAsInt(ACTIVE_REGION_END_KEY,-1);
                if (candidateStop == -1)
                    throw new StingException("Missing " + ACTIVE_REGION_END_KEY + " tag in input active regions ");
                if (originalActiveRegion.getLocation().getStart() == candidateStart && originalActiveRegion.getLocation().getStop() == candidateStop) {
                    generatedActiveRegion = candidate;
                    break;
                }
            }
        }
        return generatedActiveRegion;
    }

    /**
     * Returns matched pairs of detected bad-site and previously generated ones. Pairs will have null at first or second
     * if there is no match for either bad-site type.
     */
    private List<Pair<VariantContext, VariantContext>> badSiteMatch(final List<VariantContext> directBadSiteVCs,
                                                                    final List<VariantContext> generatedBadSiteVCs) {
        final Map<BadSite,VariantContext> directBadSites = new LinkedHashMap<>(directBadSiteVCs.size());
        final Map<BadSite,VariantContext> generatedBadSites = new LinkedHashMap<>(generatedBadSiteVCs.size());
        final List<Pair<VariantContext,VariantContext>> result = new ArrayList<>(directBadSiteVCs.size() + generatedBadSiteVCs.size());
        for (final VariantContext vc : directBadSiteVCs)
            directBadSites.put(new BadSite(vc),vc);
        for (final VariantContext vc : generatedBadSiteVCs)
            generatedBadSites.put(new BadSite(vc),vc);
        for (final BadSite bs : directBadSites.keySet())
            result.add(new Pair<>(directBadSites.get(bs),generatedBadSites.remove(bs)));
        for (final BadSite bs : generatedBadSites.keySet())
            result.add(new Pair<>((VariantContext)null,generatedBadSites.get(bs)));
        return result;
    }

    /**
     * Returns the subset of variant-context bad-sites detected by graph-based method but not using logless pair-hmm.
     */
    private List<VariantContext> badSiteSetDiff(final List<VariantContext> graphBadSiteVCs, final List<VariantContext> loglessBadSiteVCs) {
        final List<BadSite> graphBadSites = new ArrayList<>(graphBadSiteVCs.size());
        for (final VariantContext vc : graphBadSiteVCs)
            graphBadSites.add(new BadSite(vc));
        final Set<BadSite> loglessBadSites = new HashSet<>(loglessBadSiteVCs.size());
        for (final VariantContext vc : loglessBadSiteVCs)
            loglessBadSites.add(new BadSite(vc));
        final List<VariantContext> result = new ArrayList<>(graphBadSiteVCs.size());
        for (final BadSite bs : graphBadSites)
            if (!loglessBadSites.contains(bs))
                result.add(bs.toVariantContext());
        return result;
    }

    /**
     * Outputs bad-site and active-region records.
     */
    private void reportBadSitesAndActiveRegions(final ActiveRegion originalActiveRegion, final List<VariantContext> relevantBadSites, final List<GATKSAMRecord> reads) {

        if (relevantBadSites.size() == 0) return;

        badSitesWriter.outputBadSites(relevantBadSites);
        badSitesWriter.flushBadSiteBuffer();
        Collections.sort(reads, readsSorter);

        if (activeReadsOutput != null)
            for (final GATKSAMRecord read : reads) {
                if (outputReads.contains(read))
                    continue;

                activeReadsOutput.addAlignment(read);
                outputReads.add(read);
            }
        if (activeRegionsOutput != null)
            activeRegionsOutput.add(buildActiveRegionVariantContext(originalActiveRegion, relevantBadSites));
    }

    /**
     * Map method to run in the generation mode.
     *
     * @param originalActiveRegion
     * @param metaDataTracker
     * @return
     */
    private List<VariantContext> mapGeneration(final ActiveRegion originalActiveRegion, final RefMetaDataTracker metaDataTracker) {
        final List<GATKSAMRecord> reads =copyReads(originalActiveRegion);
        badSitesWriter.resetActiveRegion();
        badSitesWriter.flushBadSiteBuffer();
        // First we deal with the sites in KB that have missed by both methods:
        final List<MongoVariantContext> missedSites = consensusSiteIterator.getSitesBefore(originalActiveRegion.getLocation().getStartLocation());

        final List<VariantContext> relevantBadSites = new LinkedList<>();
        removePassedReadsFromOutputBuffer(originalActiveRegion);
        badSitesWriter.setCurrentActiveRegion(originalActiveRegion);

        //Second we run the method on this active region using graph-likelihoods and loglessPairHMM
        final boolean originalUseGraphLikelihoods = likelihoodEngineImplementation == LikelihoodCalculationEngine.Implementation.GraphBased;
        likelihoodEngineImplementation = LikelihoodCalculationEngine.Implementation.GraphBased;
        final List<VariantContext> graphResults = super.map(originalActiveRegion, metaDataTracker);
        final List<VariantContext> loglessResults = new LinkedList<>();
        if (!originalUseGraphLikelihoods || !dontDiscardLoglessErrors) {
          likelihoodEngineImplementation = LikelihoodCalculationEngine.Implementation.PairHMM;
          loglessResults.addAll(super.map(originalActiveRegion, metaDataTracker));
        }
        final List<MongoVariantContext> assessmentSites = consensusSiteIterator.getSitesBefore(originalActiveRegion.getLocation().getStopLocation().incPos());
        badSitesWriter.switchOffOutput();
        if (dontDiscardLoglessErrors) {
            assessor.assessSite(graphResults, assessmentSites, onlyReviewed);
            relevantBadSites.addAll(badSitesWriter.flushBadSiteBuffer());
        } else {
            assessor.assessSite(loglessResults, assessmentSites, onlyReviewed);
            final List<VariantContext> loglessBadSites = badSitesWriter.flushBadSiteBuffer();
            assessor.assessSite(graphResults,assessmentSites,onlyReviewed);
            final List<VariantContext> graphBadSites = badSitesWriter.flushBadSiteBuffer();
            relevantBadSites.addAll(badSiteSetDiff(graphBadSites,loglessBadSites));
        }

        badSitesWriter.switchOnOutput();
        reportBadSitesAndActiveRegions(originalActiveRegion, relevantBadSites, reads);
        badSitesWriter.switchOffOutput();
        return originalUseGraphLikelihoods ? graphResults : loglessResults;
    }

    /**
     * Construct the variant context representing an active region.
     * @param originalActiveRegion the active region to represent into a variant-context.
     * @param badSitesInActiveRegion list of bad-site records generated in the active region.
     * @return never null.
     */
    private VariantContext buildActiveRegionVariantContext(final ActiveRegion originalActiveRegion, final List<VariantContext> badSitesInActiveRegion) {
        final ActiveRegionVariantContextBuilder builder = new ActiveRegionVariantContextBuilder();
        final GenomeLoc originalLoc = originalActiveRegion.getLocation();
        final GenomeLoc extendedLoc = originalActiveRegion.getExtendedLoc();
        final byte[] referenceBytes = originalActiveRegion.getReference(referenceReader,0,originalActiveRegion.getLocation().getStartLocation());
        builder.alleles(Collections.singletonList(Allele.create(referenceBytes[0], true)));
        builder.loc(originalLoc);
        builder.extendedLoc(extendedLoc);
        builder.badSites(badSitesInActiveRegion);
        return builder.make();
    }

    /**
     * Give the current active region, remove reads from the buffer that do not have any overlap with it.
     * <p>
     *     These are always found before the active region.
     * </p>
     *
     * @param originalActiveRegion
     */
    private void removePassedReadsFromOutputBuffer(final ActiveRegion originalActiveRegion) {
        // Remove output reads in buffer before the beginning of the extended active region.
        final GenomeLoc extendedLoc = originalActiveRegion.getExtendedLoc();
        final int extendedLocReferenceIndex = extendedLoc.getContigIndex();
        final int extendedLocStart = extendedLoc.getStart();
        final Iterator<GATKSAMRecord> outputReadsIterator = outputReads.iterator();
        GATKSAMRecord firstPassedRead = null;
        while (outputReadsIterator.hasNext()) {
            final GATKSAMRecord candidate = outputReadsIterator.next();
            if (candidate.getReferenceIndex() < extendedLocReferenceIndex)
                firstPassedRead = candidate;
            else if (candidate.getAlignmentEnd() < extendedLocStart)
                firstPassedRead = candidate;
            if (firstPassedRead != null) break;
        }
        if (firstPassedRead != null)
            outputReads.tailSet(firstPassedRead).clear();
    }

    @Override
    public boolean isReduceByInterval() {
        return false;
    }

    /**
     * Closes output files.
     *
     * @param i ignored.
     */
    @Override
    public void onTraversalDone(final Integer i) {

        super.onTraversalDone(i);
        if (activeReadsOutput != null) {
            activeReadsOutput.close();
        }
        if (activeRegionsOutput != null)
            activeRegionsOutput.close();
        if (badSites != null)
            badSites.close();
        if (testDataList != null)
            logger.info("Test Data entry generated: " + testDataList.size());
    }

    public static class TestCaseData {

        //public final ActiveRegion activeRegion;
        public final List<Pair<VariantContext,VariantContext>> badSiteMatches;
        public TestCaseData(final ActiveRegion oar, final VariantContext gar, final List<Pair<VariantContext, VariantContext>> matchedBadSites) {
            //activeRegion = oar;
            badSiteMatches = matchedBadSites;
        }
    }

    protected class MySiteIterator implements SiteIterator<MongoVariantContext> {

        protected int count = 0;
        protected SiteIterator<MongoVariantContext> actualIterator;
        protected NA12878KnowledgeBase kbDb;

        protected GenomeLoc lastSuccessfullLocation;

        @Requires("sm != null")
        protected MySiteIterator(final NA12878KnowledgeBase kbDb) {
            this.kbDb = kbDb;
            actualIterator = kbDb.getConsensusSites(new SiteManager(getToolkit().getGenomeLocParser(),
                    getToolkit().getIntervals(),
                    getToolkit().getMasterSequenceDictionary()));
        }

        @Override
        public void setErrorHandler(final InvalidRecordHandler<MongoVariantContext> errorHandler) {
            actualIterator.setErrorHandler(errorHandler);
        }

        @Override
        public List<MongoVariantContext> getSitesBefore(final GenomeLoc loc) {
            List<MongoVariantContext> result;
            try {
                result = actualIterator.getSitesBefore(loc);
            } catch (RuntimeException ex) {
                List<GenomeLoc> newLocs = new LinkedList<>();
                for (final GenomeLoc l : getToolkit().getIntervals()) {
                    if (lastSuccessfullLocation.isPast(l))
                        continue;
                    else if (lastSuccessfullLocation.overlapsP(l))
                        newLocs.add(getToolkit().getGenomeLocParser().createGenomeLoc(l.getContig(),lastSuccessfullLocation.getStart(),l.getStop()));
                    else
                        newLocs.add(l);
                }
                final GenomeLocSortedSet newIntervals = new GenomeLocSortedSet(getToolkit().getGenomeLocParser(),newLocs);
                actualIterator = kbDb.getConsensusSites(new SiteManager(getToolkit().getGenomeLocParser(),newIntervals,
                        getToolkit().getMasterSequenceDictionary()));
                result = actualIterator.getSitesBefore(loc);
            }
            lastSuccessfullLocation = loc;
            return result;
        }

        @Override
        public List<MongoVariantContext> getSitesAtLocation(final GenomeLoc loc) {
            return actualIterator.getSitesAtLocation(loc);
        }


        @Override
        public List<MongoVariantContext> getNextEquivalents() {
            return actualIterator.getNextEquivalents();
        }

        @Override
        public List<MongoVariantContext> toList() {
            return actualIterator.toList();
        }

        @Override
        public void close() {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public Iterator<MongoVariantContext> iterator() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public boolean hasNext() {
            return false;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public MongoVariantContext next() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public void remove() {
            //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    /**
     * Redefines the equals and hashCode methods of variant context for bad-site matching.
     */
    private class BadSite extends VariantContext {

        private final int hashCode;

        protected BadSite(final VariantContext other) {
            super(other);
            hashCode = calculateHashCode();
        }

        @Override
        public boolean equals(final Object o) {
            if (o == null)
                return false;
            else if (o == this)
                return true;
            else if (o.hashCode() != hashCode())
                return false;
            else if (!(o instanceof VariantContext))
                return false;
            else {
                final VariantContext vc = (VariantContext) o;
                return (vc.isBiallelic() == this.isBiallelic()) &&
                        vc.getChr().equals(getChr()) &&
                        vc.getStart() == getStart() &&
                        vc.getEnd() == getEnd() &&
                        vc.getReference().getBaseString().equals(getReference().getBaseString()) &&
                        vc.getAlternateAllele(0).getBaseString().equals(getAlternateAllele(0).getBaseString()) &&
                        vc.getGenotype(0).getGenotypeString().equals(getGenotype(0).getGenotypeString());
            }
        }

        @Override
        public int hashCode() {
            return hashCode;
        }

        // Calculates a ok hash code.
        private int calculateHashCode() {
            final String uniqueString = "" + isBiallelic() + "-" + getChr() + "-" + getStart() + "-" +
                    getEnd() + "-" + getReference().getBaseString() +
                    "-" + getAlternateAllele(0).getBaseString() + "-" + getGenotype(0).getGenotypeString();
            return uniqueString.hashCode();
        }

        /**
         * Returns a {@link VariantContext} casted representation of this BadSite.
         * @return never {@code null}.
         */
        public VariantContext toVariantContext() {
            return this;
        }




    }
}
