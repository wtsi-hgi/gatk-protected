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

import htsjdk.samtools.SAMFileWriter;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.engine.arguments.DbsnpArgumentCollection;
import org.broadinstitute.gatk.engine.filters.BadMateFilter;
import org.broadinstitute.gatk.engine.io.DirectOutputTracker;
import org.broadinstitute.gatk.engine.io.stubs.SAMFileWriterStub;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypingOutputMode;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.FixedAFCalculatorProvider;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.*;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegion;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegionReadState;
import org.broadinstitute.gatk.utils.activeregion.ActivityProfileState;
import org.broadinstitute.gatk.utils.clipping.ReadClipper;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContextUtils;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.DownsamplingUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.fragments.FragmentCollection;
import org.broadinstitute.gatk.utils.fragments.FragmentUtils;
import org.broadinstitute.gatk.utils.genotyper.*;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.haplotypeBAMWriter.HaplotypeBAMWriter;
import org.broadinstitute.gatk.utils.pairhmm.PairHMM;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.sam.*;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

import static java.lang.Math.pow;

@PartitionBy(PartitionType.LOCUS)
@BAQMode(ApplicationTime = ReadTransformer.ApplicationTime.FORBIDDEN)
@ActiveRegionTraversalParameters(extension=100, maxRegion=300)
public class M2 extends ActiveRegionWalker<List<VariantContext>, Integer> implements AnnotatorCompatible, NanoSchedulable {
    public static final String BAM_TAG_TUMOR = "tumor";
    public static final String BAM_TAG_NORMAL = "normal";

    protected Set<SAMReaderID> tumorSAMReaderIDs = new HashSet<>();
    protected Set<SAMReaderID> normalSAMReaderIDs = new HashSet<>();
    protected String tumorSampleName;
    protected String normalSampleName;

    protected SampleList samplesList;

    // fasta reference reader to supplement the edges of the reference sequence
    protected CachingIndexedFastaSequenceFile referenceReader;

    // the assembly engine
    protected LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    protected ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;

    // the genotyping engine
    protected HaplotypeCallerGenotypingEngine genotypingEngine = null;


    private byte MIN_TAIL_QUALITY;
    private double log10GlobalReadMismappingRate;



    @ArgumentCollection
    protected M2ArgumentCollection MTAC = new M2ArgumentCollection();

    @Argument(fullName = "debug_read_name", required = false, doc="trace this read name through the calling process")
    public String DEBUG_READ_NAME = null;


    /***************************************/
    // Reference Metadata inputs
    /***************************************/
    @Input(fullName="cosmic", shortName = "cosmic", doc="VCF file of COSMIC sites", required=false)
    public List<RodBinding<VariantContext>> cosmicRod = Collections.emptyList();

    @Input(fullName="normal_panel", shortName = "PON", doc="VCF file of sites observed in normal", required=false)
    public List<RodBinding<VariantContext>> normalPanelRod = Collections.emptyList();

    @Override
    public void initialize() {
        super.initialize();

        samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSAMFileSamples(getToolkit().getSAMFileHeader())));

        // MUTECT: check that we have at least one tumor bam
        for(SAMReaderID id : getToolkit().getReadsDataSource().getReaderIDs()) {
            if (id.getTags().getPositionalTags().size() == 0) {
                throw new RuntimeException("BAMs must be tagged as either 'tumor' or 'normal'");
            }

            // only supports single-sample BAMs (ie first read group is representative)
            String bamSampleName = getToolkit().getReadsDataSource().getHeader(id).getReadGroups().get(0).getSample();

            for(String tag : id.getTags().getPositionalTags()) {
                if (BAM_TAG_TUMOR.equalsIgnoreCase(tag)) {
                    tumorSAMReaderIDs.add(id);
                    if (tumorSampleName == null) {
                        tumorSampleName = bamSampleName;
                    } else {
                        if (!tumorSampleName.equals(bamSampleName)) {
                            throw new UserException.BadInput("Found more than one tumor sample name in read data");
                        }
                    }
                } else if (BAM_TAG_NORMAL.equalsIgnoreCase(tag)) {
                    normalSAMReaderIDs.add(id);
                    if (normalSampleName == null) {
                        normalSampleName = bamSampleName;
                    } else {
                        if (!normalSampleName.equals(bamSampleName)) {
                            throw new UserException.BadInput("Found more than one normal sample name in read data");
                        }
                    }
                } else {
                    throw new RuntimeException("Unknown BAM tag '" + tag + "' must be either 'tumor' or 'normal'");
                }
            }
        }

        final VariantAnnotatorEngine annotationEngine = initializeVCFOutput();

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile, e);
        }

        // create and setup the assembler
        assemblyEngine = new ReadThreadingAssembler(maxNumHaplotypesInPopulation, kmerSizes, dontIncreaseKmerSizesForCycles, allowNonUniqueKmersInRef, numPruningSamples);

        assemblyEngine.setErrorCorrectKmers(errorCorrectKmers);
        assemblyEngine.setPruneFactor(MIN_PRUNE_FACTOR);
        assemblyEngine.setDebug(SCAC.DEBUG);
        assemblyEngine.setDebugGraphTransformations(debugGraphTransformations);
        assemblyEngine.setAllowCyclesInKmerGraphToGeneratePaths(allowCyclesInKmerGraphToGeneratePaths);
        assemblyEngine.setRecoverDanglingBranches(!doNotRecoverDanglingBranches);
        assemblyEngine.setMinBaseQualityToUseInAssembly(MIN_BASE_QUALTY_SCORE);

        MIN_TAIL_QUALITY = (byte)(MIN_BASE_QUALTY_SCORE - 1);

        if ( graphWriter != null ) assemblyEngine.setGraphWriter(graphWriter);

        // setup the likelihood calculation engine
        if ( phredScaledGlobalReadMismappingRate < 0 ) phredScaledGlobalReadMismappingRate = -1;

        // configure the global mismapping rate
        if ( phredScaledGlobalReadMismappingRate < 0 ) {
            log10GlobalReadMismappingRate = - Double.MAX_VALUE;
        } else {
            log10GlobalReadMismappingRate = QualityUtils.qualToErrorProbLog10(phredScaledGlobalReadMismappingRate);
            logger.info("Using global mismapping rate of " + phredScaledGlobalReadMismappingRate + " => " + log10GlobalReadMismappingRate + " in log10 likelihood units");
        }

        //static member function - set number of threads
        PairHMM.setNumberOfThreads(getToolkit().getTotalNumberOfThreads());
        // create our likelihood calculation engine
        likelihoodCalculationEngine = createLikelihoodCalculationEngine();

        final MergeVariantsAcrossHaplotypes variantMerger = new MergeVariantsAcrossHaplotypes();

        final GenomeAnalysisEngine toolkit = getToolkit();
        final GenomeLocParser genomeLocParser = toolkit.getGenomeLocParser();

        genotypingEngine = new SomaticGenotypingEngine( SCAC, samplesList, genomeLocParser, FixedAFCalculatorProvider.createThreadSafeProvider(getToolkit(), SCAC, logger), !doNotRunPhysicalPhasing, MTAC);

        genotypingEngine.setCrossHaplotypeEventMerger(variantMerger);
        genotypingEngine.setAnnotationEngine(annotationEngine);


        if ( bamWriter != null ) {
            // we currently do not support multi-threaded BAM writing, so exception out
            if ( getToolkit().getTotalNumberOfThreads() > 1 )
                throw new UserException.BadArgumentValue("bamout", "Currently cannot emit a BAM file from the HaplotypeCaller in multi-threaded mode.");
            haplotypeBAMWriter = HaplotypeBAMWriter.create(bamWriterType, bamWriter, getToolkit().getSAMFileHeader());
        }

        // why isn't this a constructor (instead of initialize)?  Since the method is package-friendly
        trimmer.initialize(getToolkit().getGenomeLocParser(), SCAC.DEBUG,
                SCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES, false);

        // KCIBUL: what's the right way to set this sensible default for somatic mutation calling from here?
        trimmer.snpPadding = 50;



    }

    private VariantAnnotatorEngine initializeVCFOutput() {
        // initialize the output VCF header
        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationClassesToUse), annotationsToUse, annotationsToExclude, this, getToolkit());

        Set<VCFHeaderLine> headerInfo = new HashSet<>();

        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        // MuTect
        headerInfo.add(new VCFInfoHeaderLine(SomaticGenotypingEngine.NORMAL_LOD, 1, VCFHeaderLineType.String, "Normal LOD score"));
        headerInfo.add(new VCFInfoHeaderLine(SomaticGenotypingEngine.TUMOR_LOD, 1, VCFHeaderLineType.String, "Tumor LOD score"));
        headerInfo.add(new VCFInfoHeaderLine("PON", 1, VCFHeaderLineType.String, "Count from Panel of Normals"));

        headerInfo.add(new VCFInfoHeaderLine(SomaticGenotypingEngine.HAPLOTYPE_COUNT, 1, VCFHeaderLineType.String, "Number of haplotypes that support this variant"));
        headerInfo.add(new VCFInfoHeaderLine("ECNT", 1, VCFHeaderLineType.String, "Number of events in this haplotype"));
        headerInfo.add(new VCFInfoHeaderLine("MIN_ED", 1, VCFHeaderLineType.Integer, "Minimum distance between events in this active region"));
        headerInfo.add(new VCFInfoHeaderLine("MAX_ED", 1, VCFHeaderLineType.Integer, "Maximum distance between events in this active region"));
        headerInfo.add(new VCFFormatHeaderLine("AF", 1, VCFHeaderLineType.Float, "Allele fraction of the event in the tumor"));

        headerInfo.add(new VCFFilterHeaderLine("panel_of_normals", "Seen in at least 2 samples in the panel of normals"));
        headerInfo.add(new VCFFilterHeaderLine("alt_allele_in_normal", "Evidence seen in the normal sample"));
        headerInfo.add(new VCFFilterHeaderLine("multi_event_alt_allele_in_normal", "TBD"));
        headerInfo.add(new VCFFilterHeaderLine("homologous_mapping_event", "TBD"));
        headerInfo.add(new VCFFilterHeaderLine("clustered_events", "TBD"));

        headerInfo.add(new VCFFilterHeaderLine("t_lod_fstar", "TBD"));
        headerInfo.add(new VCFFilterHeaderLine("germline_risk", "TBD"));

        vcfWriter.writeHeader(new VCFHeader(headerInfo, SampleListUtils.asList(samplesList)));

        return annotationEngine;
    }

    @Override
    public ActivityProfileState isActive(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if( context == null || context.getBasePileup().isEmpty() )
            // if we don't have any data, just abort early
            return new ActivityProfileState(ref.getLocus(), 0.0);

        final Map<String, AlignmentContext> splitContexts = AlignmentContextUtils.splitContextBySampleName(context);
        AlignmentContext tumorContext = splitContexts.get(tumorSampleName);
        AlignmentContext normalContext = splitContexts.get(normalSampleName);

        // if there are no tumor reads... there is no activity!
        if (tumorContext == null) {
            return new ActivityProfileState(ref.getLocus(), 0);
        }

        // KCIBUL -- this method was inlined and modified from ReferenceConfidenceModel
        ReadBackedPileup tumorPileup = tumorContext.getBasePileup().getMappingFilteredPileup(20);
        final double[] tumorGLs = calcGenotypeLikelihoodsOfRefVsAny(tumorPileup, ref.getBase(), MIN_BASE_QUALTY_SCORE);
        final double tumorLod = tumorGLs[1] - tumorGLs[0];

        // NOTE: do I want to convert to a probability (or just keep this as a LOD score)

        // also at this point, we should throw out noisy sites (hence the nonRefInNormalCheck) but this is non-optimal
        double prob = 0;
        if (tumorLod > MTAC.INITIAL_TUMOR_LOD_THRESHOLD) {

            // TODO: should we even do this performance optimization?
            // in any case, we have to handle the case where there is no normal (and thus no normal context) which is
            // different than having a normal but having no reads (where we should not enter the active region)
            if (normalSampleName != null && normalContext != null) {
                int nonRefInNormal = getCountOfNonRefEvents(normalContext.getBasePileup(), ref.getBase(), MIN_BASE_QUALTY_SCORE);

                final double[] normalGLs = calcGenotypeLikelihoodsOfRefVsAny(normalContext.getBasePileup(), ref.getBase(), MIN_BASE_QUALTY_SCORE, 0.5f);
                final double normalLod = normalGLs[0] - normalGLs[1];

                // TODO: parameterize these
                if (normalLod > 1.0 && nonRefInNormal < 4) {
                    prob = 1;
                    logger.debug("At " + ref.getLocus().toString() + " tlod: " + tumorLod + " nlod: " + normalLod + " with normal non-ref of " + nonRefInNormal);
                }
            } else {
                prob = 1;
                logger.debug("At " + ref.getLocus().toString() + " tlod: " + tumorLod + " and no-normal calling");
            }

        }

        return new ActivityProfileState( ref.getLocus(), prob, ActivityProfileState.Type.NONE, null);
    }

    private final static List<VariantContext> NO_CALLS = Collections.emptyList();
    @Override
    public List<VariantContext> map( final ActiveRegion originalActiveRegion, final RefMetaDataTracker metaDataTracker ) {
        if ( justDetermineActiveRegions )
            // we're benchmarking ART and/or the active region determination code in the HC, just leave without doing any work
            return NO_CALLS;

        if( !originalActiveRegion.isActive() )
            // Not active so nothing to do!
            return referenceModelForNoVariation(originalActiveRegion, true);

        // No reads here so nothing to do!
        if( originalActiveRegion.size() == 0 ) { return referenceModelForNoVariation(originalActiveRegion, true); }

        logReadInfo(DEBUG_READ_NAME, originalActiveRegion.getReads(), "Present in original active region");

        // create the assembly using just high quality reads (Q20 or higher).  We want to use lower
        // quality reads in the PairHMM (and especially in the normal) later, so we can't use a ReadFilter
        ActiveRegion assemblyActiveRegion = new ActiveRegion(originalActiveRegion.getLocation(), originalActiveRegion.getSupportingStates(),originalActiveRegion.isActive(), getToolkit().getGenomeLocParser(), originalActiveRegion.getExtension());
        for (GATKSAMRecord rec : originalActiveRegion.getReads()) {
            if (rec.getMappingQuality() >= 20 ) {
                assemblyActiveRegion.add(rec);
            }
        }

        logReadInfo(DEBUG_READ_NAME, assemblyActiveRegion.getReads(), "Present in assembly active region");

        // run the local assembler, getting back a collection of information on how we should proceed
        final List<VariantContext> givenAlleles = new ArrayList<>();
        final AssemblyResultSet untrimmedAssemblyResult = assembleReads(assemblyActiveRegion, givenAlleles);


        final TreeSet<VariantContext> allVariationEvents = untrimmedAssemblyResult.getVariationEvents();
        // TODO - line bellow might be unecessary : it might be that assemblyResult will always have those alleles anyway
        // TODO - so check and remove if that is the case:
        allVariationEvents.addAll(givenAlleles);

        final ActiveRegionTrimmer.Result trimmingResult = trimmer.trim(originalActiveRegion,allVariationEvents);


        // Stop the trimming madness!!!
        if (!trimmingResult.isVariationPresent())
            return referenceModelForNoVariation(originalActiveRegion,false);

        logReadInfo(DEBUG_READ_NAME, trimmingResult.getCallableRegion().getReads(), "Present in trimming result");

        final AssemblyResultSet assemblyResult =
                trimmingResult.needsTrimming() ? untrimmedAssemblyResult.trimTo(trimmingResult.getCallableRegion()) : untrimmedAssemblyResult;

//        final AssemblyResultSet assemblyResult = untrimmedAssemblyResult;

        // after talking to Ryan -- they grab the reads out of the assembly (and trim then) to pass into the PairHMM
        // because at one point they were trying error correcting of the reads based on the haplotypes.. but that is not
        // working out, so it's safe for us just to take the reads
//
        final ActiveRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();
        logReadInfo(DEBUG_READ_NAME, regionForGenotyping.getReads(), "Present in region for genotyping");

//
//        final ActiveRegion regionForGenotyping = trimmingResult.getCallableRegion();

//        final ActiveRegion regionForGenotyping = originalActiveRegion;

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        final Collection<GATKSAMRecord> filteredReads = filterNonPassingReads( regionForGenotyping );

        logReadInfo(DEBUG_READ_NAME, regionForGenotyping.getReads(), "Present in region for genotyping after filtering reads");

        final Map<String, List<GATKSAMRecord>> perSampleFilteredReadList = splitReadsBySample( filteredReads );

        // abort early if something is out of the acceptable range
        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if( ! assemblyResult.isVariationPresent() )
            return referenceModelForNoVariation(originalActiveRegion, false);

        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if( regionForGenotyping.size() == 0 ) {
            // no reads remain after filtering so nothing else to do!
            return referenceModelForNoVariation(originalActiveRegion, false);
        }

        // evaluate each sample's reads against all haplotypes

        final List<Haplotype> haplotypes = assemblyResult.getHaplotypeList();
        final Map<String,List<GATKSAMRecord>> reads = splitReadsBySample( regionForGenotyping.getReads() );

        // modify MAPQ scores in normal to be high so that we don't do any base quality score capping
        for(GATKSAMRecord rec : regionForGenotyping.getReads()) {
            if (isReadFromNormal(rec)) {
                rec.setMappingQuality(60);
            }
        }

        logger.debug("Computing read likelihoods with " + regionForGenotyping.getReads().size() + " reads against " + haplotypes.size() + " haplotypes across region " + assemblyResult.getRegionForGenotyping().toString());


        // Calculate the likelihoods: CPU intensive part.
        final ReadLikelihoods<Haplotype> readLikelihoods =
                likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult,samplesList,reads);

        // Realign reads to their best haplotype.
        // KCIBUL: this is new stuff -- review it!
        final Map<GATKSAMRecord,GATKSAMRecord> readRealignments = realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc());
        readLikelihoods.changeReads(readRealignments);

        for (GATKSAMRecord rec : readRealignments.keySet()) {
            logReadInfo(DEBUG_READ_NAME, rec, "Present after computing read likelihoods");
        }

        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = ((SomaticGenotypingEngine)genotypingEngine).callMutations(
                haplotypes,
                readLikelihoods,
                perSampleFilteredReadList,
                assemblyResult.getFullReferenceWithPadding(),
                assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getLocation(),
                getToolkit().getGenomeLocParser(),
                metaDataTracker,
                givenAlleles, false ,
                tumorSampleName,
                normalSampleName,
                dbsnp.dbsnp,
                cosmicRod,
                DEBUG_READ_NAME
                );

        if ( bamWriter != null ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            if (disableOptimizations)
                calledHaplotypeSet.add(assemblyResult.getReferenceHaplotype());
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                    haplotypes,
                    assemblyResult.getPaddedReferenceLoc(),
                    haplotypes,
                    calledHaplotypeSet,
                    readLikelihoods);
        }

        if( SCAC.DEBUG ) { logger.info("----------------------------------------------------------------------------------"); }


        List<VariantContext> annotatedCalls = new ArrayList<>();
        int eventCount = calledHaplotypes.getCalls().size();
        Integer minEventDistance = null;
        Integer maxEventDistance = null;
        Integer lastPosition = null;
        for (VariantContext vc : calledHaplotypes.getCalls()) {
            if (lastPosition == null) {
                lastPosition = vc.getStart();
            } else {
                int dist = Math.abs(vc.getStart() - lastPosition);
                if (maxEventDistance == null || dist > maxEventDistance) {
                    maxEventDistance = dist;
                }
                if (minEventDistance == null || dist < minEventDistance) {
                    minEventDistance = dist;
                }
            }
        }

        // can we do this with the Annotation classes instead?
        for (VariantContext originalVC : calledHaplotypes.getCalls()) {

            // TODO: clean up filtering logic... there must be a better way...
            boolean wasFiltered = originalVC.getFilters().size() > 0;
            VariantContextBuilder vcb = new VariantContextBuilder(originalVC);

            vcb.attribute("ECNT", eventCount);
            vcb.attribute("MIN_ED", minEventDistance);
            vcb.attribute("MAX_ED", maxEventDistance);

            Collection<VariantContext> panelOfNormalsVC = metaDataTracker.getValues(normalPanelRod,
                    getToolkit().getGenomeLocParser().createGenomeLoc(originalVC.getChr(), originalVC.getStart()));
            VariantContext ponVc = panelOfNormalsVC.isEmpty()?null:panelOfNormalsVC.iterator().next();

            if (ponVc != null) {
                vcb.filter("panel_of_normals");
                wasFiltered = true;
            }

            // FIXME: how do we sum qscores here?
            // FIXME: parameterize thresholds
            // && sum of alt likelihood scores > 20

            // TODO: make the change to have only a single normal sample (but multiple tumors is ok...)
            int normalAltCounts = 0;
            double normalF = 0;
            int normalAltQualityScoreSum = 0;
            if (hasNormal()) {
                Genotype normalGenotype = originalVC.getGenotype(normalSampleName);

                // NOTE: how do we get the non-ref depth here?
                normalAltCounts = normalGenotype.getAD()[1];
                normalF = (Double) normalGenotype.getExtendedAttribute("AF");

                Object qss = normalGenotype.getExtendedAttribute("QSS");
                if (qss == null) {
                    logger.error("Null qss at " + originalVC.getStart());
                }
                normalAltQualityScoreSum = (Integer) ((Object[]) qss)[1];
            }

            if ( (normalAltCounts >= MTAC.MAX_ALT_ALLELES_IN_NORMAL_COUNT || normalF >= MTAC.MAX_ALT_ALLELE_IN_NORMAL_FRACTION ) && normalAltQualityScoreSum > MTAC.MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM) {
                vcb.filter("alt_allele_in_normal");
                wasFiltered = true;
            } else {

                // NOTE: does normal alt counts presume the normal had all these events in CIS?
                if ( eventCount > 1 && normalAltCounts >= 1) {
                    vcb.filter("multi_event_alt_allele_in_normal");
                    wasFiltered = true;
                } else if (eventCount >= 3) {
                    vcb.filter("homologous_mapping_event");
                    wasFiltered = true;
                }

            }

            // NOTE: what if there is a 3bp indel followed by a snp... we are comparing starts
            // so it would be thrown but it's really an adjacent event
            if ( eventCount >= 2 && maxEventDistance >= 3) {
                vcb.filter("clustered_events");
                wasFiltered = true;
            }

            double tumorLod = (Double) originalVC.getAttributeAsDouble(SomaticGenotypingEngine.TUMOR_LOD, -1);
            if (tumorLod < MTAC.TUMOR_LOD_THRESHOLD) {
                vcb.filter("t_lod_fstar");
                wasFiltered = true;
            }

            if (!wasFiltered) vcb.passFilters();
            annotatedCalls.add(vcb.make());
        }






        return annotatedCalls;
    }


    private final static byte REF_MODEL_DELETION_QUAL = (byte) 30;
    /**
     * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
     *
     * @param pileup the read backed pileup containing the data we want to evaluate
     * @param refBase the reference base at this pileup position
     * @param minBaseQual the min base quality for a read in the pileup at the pileup position to be included in the calculation
     * @return genotype likelihoods of [AA,AB]
     */
    protected double[] calcGenotypeLikelihoodsOfRefVsAny(final ReadBackedPileup pileup, final byte refBase, final byte minBaseQual, final double f) {
        final double[] genotypeLikelihoods = new double[2];
        int AA = 0, AB=1;
        for( final PileupElement p : pileup ) {
            final byte qual = (p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual());
            if( p.isDeletion() || qual > minBaseQual ) {

                // TODO: why not use base qualities here?
                //double pobs = QualityUtils.qualToErrorProbLog10(qual);
                double pobs = 1.0d - pow(10, (30 / -10.0));
                if( isNonRef(refBase, p)) {
                    genotypeLikelihoods[AB] += Math.log10(f*pobs + (1-f)*pobs/3.0d);
                    genotypeLikelihoods[AA] += Math.log10((1-pobs)/3);
                } else {
                    genotypeLikelihoods[AB] += Math.log10(f*(1-pobs)/3.0d + (1-f)*pobs);
                    genotypeLikelihoods[AA] += Math.log10(pobs);
                }
            }
        }

        return genotypeLikelihoods;
    }

    private boolean hasNormal() {
        return (normalSampleName != null);
    }

    protected int getCountOfNonRefEvents(final ReadBackedPileup pileup, final byte refBase, final byte minBaseQual) {
        int i=0;
        for( final PileupElement p : pileup ) {
            final byte qual = (p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual());
            if( p.isDeletion() || qual > minBaseQual ) {
                if( isNonRef(refBase, p)) {
                    i++;
                }
            }
        }
        return i;
    }

    protected double[] calcGenotypeLikelihoodsOfRefVsAny(final ReadBackedPileup pileup, final byte refBase, final byte minBaseQual) {
        double f = calculateF(pileup, refBase, minBaseQual);
        return calcGenotypeLikelihoodsOfRefVsAny(pileup, refBase, minBaseQual, f);
    }

    private double calculateF(final ReadBackedPileup pileup, final byte refBase, final byte minBaseQual) {
        int refCount = 0, altCount = 0;
        for( final PileupElement p : pileup ) {
            final byte qual = (p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual());

            // only consider deletions AND sites of sufficient quality
            if( p.isDeletion() || qual > minBaseQual ) {
                if( isNonRef(refBase, p)) {
                    altCount++;
                } else {
                    refCount++;
                }
            }
        }
        double f = (double) altCount / ((double) refCount + (double) altCount);
        return f;
    }

    private boolean isNonRef(byte refBase, PileupElement p) {
        return p.getBase() != refBase || p.isDeletion() || p.isBeforeDeletionStart() || p.isAfterDeletionEnd() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip();
    }

    int MIN_READ_LENGTH = 30; // private in superclass

    protected Set<GATKSAMRecord> filterNonPassingReads( final ActiveRegion activeRegion) {
        final Set<GATKSAMRecord> readsToRemove = new LinkedHashSet<>();
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {

            // KCIBUL: only perform read quality filtering on tumor reads...
            if (isReadFromNormal(rec)) {

                if( rec.getReadLength() < MIN_READ_LENGTH ) {
                    readsToRemove.add(rec);
                }

            } else {



                if( rec.getReadLength() < MIN_READ_LENGTH ||
                    rec.getMappingQuality() < 20 ||
                    BadMateFilter.hasBadMate(rec) ||

                    (keepRG != null && !rec.getReadGroup().getId().equals(keepRG)) ) {
                    readsToRemove.add(rec);
                }
            }
        }
        activeRegion.removeAll( readsToRemove );
        return readsToRemove;
    }

    private static GATKSAMRecord findReadByName(Collection<GATKSAMRecord> reads, String name) {
        for(GATKSAMRecord read : reads) {
            if (name.equals(read.getReadName())) return read;
        }
        return null;
    }

    /**
     * Instantiates the appropriate likelihood calculation engine.
     *
     * @return never {@code null}.
     */
    private ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine() {
        return new PairHMMLikelihoodCalculationEngine( (byte)gcpHMM, pairHMM, pairHMMSub, alwaysLoadVectorLoglessPairHMMLib, log10GlobalReadMismappingRate, noFpga, pcrErrorModel );
    }

    /**
     * FROM HC
     *
     * Create an ref model result (ref model or no calls depending on mode) for an active region without any variation
     * (not is active, or assembled to just ref)
     *
     * @param region the region to return a no-variation result
     * @param needsToBeFinalized should the region be finalized before computing the ref model (should be false if already done)
     * @return a list of variant contexts (can be empty) to emit for this ref region
     */
    protected List<VariantContext> referenceModelForNoVariation(final ActiveRegion region, final boolean needsToBeFinalized) {
            return NO_CALLS;
    }

    protected Map<String, List<GATKSAMRecord>> splitReadsBySample( final Collection<GATKSAMRecord> reads ) {
        return HaplotypeCaller.splitReadsBySample(samplesList, reads);
    }

    // enable deletions in the pileup
    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable non primary and extended reads in the active region
    @Override
    public EnumSet<ActiveRegionReadState> desiredReadStates() {
//        if ( includeUnmappedReads )
//            throw new UserException.BadArgumentValue("includeUnmappedReads", "is not yet functional");
//        else
            return EnumSet.of(
                    ActiveRegionReadState.PRIMARY,
                    ActiveRegionReadState.NONPRIMARY,
                    ActiveRegionReadState.EXTENDED);
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(List<VariantContext> callsInRegion, Integer numCalledRegions) {
        for( final VariantContext call : callsInRegion ) {
            vcfWriter.add( call );
        }
        return (callsInRegion.isEmpty() ? 0 : 1) + numCalledRegions;
    }

    @Override
    public void onTraversalDone(Integer result) {
//        if ( SCAC.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) ((GVCFWriter)vcfWriter).close(false); // GROSS -- engine forces us to close our own VCF writer since we wrapped it
//        referenceConfidenceModel.close();
        //TODO remove the need to call close here for debugging, the likelihood output stream should be managed
        //TODO (open & close) at the walker, not the engine.
        likelihoodCalculationEngine.close();
        logger.info("Ran local assembly on " + result + " active regions");
    }


    // The following are not used but are required by the AnnotatorCompatible interface
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }

    /**
     * rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    public RodBinding<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     *  as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field).
     *  Records that are filtered in the comp track will be ignored.
     *  Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Advanced
    @Input(fullName="comp", shortName = "comp", doc="comparison VCF file", required=false)
    public List<RodBinding<VariantContext>> comps = Collections.emptyList();
    public List<RodBinding<VariantContext>> getCompRodBindings() { return comps; }



    /**
     * Which annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available annotations.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
//    protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{"ClippingRankSumTest", "DepthPerSampleHC"}));
    protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{"DepthPerAlleleBySample", "BaseQualitySumPerAlleleBySample"}));

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the -A or -G arguments,
     * so annotations will be excluded even if they are explicitly included with the other options.
     */
    @Advanced
    @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", required=false)
    protected List<String> annotationsToExclude = new ArrayList<>(Arrays.asList(new String[]{"SpanningDeletions", "TandemRepeatAnnotator"}));

    /**
     * Which groups of annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available groups.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
//    protected String[] annotationClassesToUse = { "Standard" };
    protected String[] annotationClassesToUse = { };

    /// HC-Related Parameters

    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    @ArgumentCollection
    protected HaplotypeCallerArgumentCollection SCAC = new HaplotypeCallerArgumentCollection();

    /**
     * Active region trimmer reference.
     */
    @ArgumentCollection
    protected ActiveRegionTrimmer trimmer = new ActiveRegionTrimmer();

    @Hidden
    @Argument(fullName="keepRG", shortName="keepRG", doc="Only use read from this read group when making calls (but use all reads to build the assembly)", required = false)
    protected String keepRG = null;



    /// Assembler parameters

    /**
     * Assembly graph can be quite complex, and could imply a very large number of possible haplotypes.  Each haplotype
     * considered requires N PairHMM evaluations if there are N reads across all samples.  In order to control the
     * run of the haplotype caller we only take maxNumHaplotypesInPopulation paths from the graph, in order of their
     * weights, no matter how many paths are possible to generate from the graph.  Putting this number too low
     * will result in dropping true variation because paths that include the real variant are not even considered.
     */
    @Advanced
    @Argument(fullName="maxNumHaplotypesInPopulation", shortName="maxNumHaplotypesInPopulation", doc="Maximum number of haplotypes to consider for your population. This number will probably need to be increased when calling organisms with high heterozygosity.", required = false)
    protected int maxNumHaplotypesInPopulation = 128;

    @Hidden
    @Argument(fullName="errorCorrectKmers", shortName="errorCorrectKmers", doc = "Use an exploratory algorithm to error correct the kmers used during assembly.  May cause fundamental problems with the assembly graph itself", required=false)
    protected boolean errorCorrectKmers = false;

    /**
     * The minimum confidence needed for a given base for it to be used in variant calling.
     */
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public byte MIN_BASE_QUALTY_SCORE = 10;

    /**
     * Users should be aware that this argument can really affect the results of the variant calling and should exercise caution.
     * Using a prune factor of 1 (or below) will prevent any pruning from the graph which is generally not ideal; it can make the
     * calling much slower and even less accurate (because it can prevent effective merging of "tails" in the graph).  Higher values
     * tend to make the calling much faster, but also lowers the sensitivity of the results (because it ultimately requires higher
     * depth to produce calls).
     */
    @Advanced
    @Argument(fullName="minPruning", shortName="minPruning", doc = "The minimum allowed pruning factor in assembly graph. Paths with < X supporting kmers are pruned from the graph", required = false)
    protected int MIN_PRUNE_FACTOR = 2;

    @Hidden
    @Argument(fullName="allowCyclesInKmerGraphToGeneratePaths", shortName="allowCyclesInKmerGraphToGeneratePaths", doc="If specified, we will allow cycles in the kmer graphs to generate paths with multiple copies of the path sequenece rather than just the shortest paths", required = false)
    protected boolean allowCyclesInKmerGraphToGeneratePaths = false;

    @Hidden
    @Argument(fullName="debugGraphTransformations", shortName="debugGraphTransformations", doc="If specified, we will write DOT formatted graph files out of the assembler for only this graph size", required = false)
    protected boolean debugGraphTransformations = false;
    /**
     * By default, the read threading assembler will attempt to recover dangling heads and tails. See the `minDanglingBranchLength` argument documentation for more details.
     */
    @Hidden
    @Argument(fullName="doNotRecoverDanglingBranches", shortName="doNotRecoverDanglingBranches", doc="Disable dangling head and tail recovery", required = false)
    protected boolean doNotRecoverDanglingBranches = false;






    // -----------------------------------------------------------------------------------------------
    // arguments to control internal behavior of the read threading assembler
    // -----------------------------------------------------------------------------------------------

    @Advanced
    @Argument(fullName="kmerSize", shortName="kmerSize", doc="Kmer size to use in the read threading assembler", required = false)
    protected List<Integer> kmerSizes = Arrays.asList(10, 25);

    @Advanced
    @Argument(fullName="dontIncreaseKmerSizesForCycles", shortName="dontIncreaseKmerSizesForCycles", doc="Should we disable the iterating over kmer sizes when graph cycles are detected?", required = false)
    protected boolean dontIncreaseKmerSizesForCycles = false;

    @Advanced
    @Argument(fullName="allowNonUniqueKmersInRef", shortName="allowNonUniqueKmersInRef", doc="Should we allow graphs which have non-unique kmers in the reference?", required = false)
    protected boolean allowNonUniqueKmersInRef = false;

    @Advanced
    @Argument(fullName="numPruningSamples", shortName="numPruningSamples", doc="The number of samples that must pass the minPuning factor in order for the path to be kept", required = false)
    protected int numPruningSamples = 1;




    @Output(fullName="graphOutput", shortName="graph", doc="File to which debug assembly graph information should be written", required = false, defaultToStdout = false)
    protected PrintStream graphWriter = null;

    /**
     * If set, certain "early exit" optimizations in HaplotypeCaller, which aim to save compute and time by skipping
     * calculations if an ActiveRegion is determined to contain no variants, will be disabled. This is most likely to be useful if
     * you're using the -bamout argument to examine the placement of reads following reassembly and are interested in seeing the mapping of
     * reads in regions with no variations. Setting the -forceActive and -dontTrimActiveRegions flags may also be necessary.
     */
    @Advanced
    @Argument(fullName = "disableOptimizations", shortName="disableOptimizations", doc="Don't skip calculations in ActiveRegions with no variants",
            required = false)
    private boolean disableOptimizations = false;

    /**
     * The assembled haplotypes will be written as BAM to this file if requested.  Really for debugging purposes only.
     * Note that the output here does not include uninformative reads so that not every input read is emitted to the bam.
     *
     * Turning on this mode may result in serious performance cost for the HC.  It's really only appropriate to
     * use in specific areas where you want to better understand why the HC is making specific calls.
     *
     * The reads are written out containing a HC tag (integer) that encodes which haplotype each read best matches
     * according to the haplotype caller's likelihood calculation.  The use of this tag is primarily intended
     * to allow good coloring of reads in IGV.  Simply go to Color Alignments By > Tag and enter HC to more
     * easily see which reads go with these haplotype.
     *
     * Note that the haplotypes (called or all, depending on mode) are emitted as single reads covering the entire
     * active region, coming from read HC and a special read group.
     *
     * Note that only reads that are actually informative about the haplotypes are emitted.  By informative we mean
     * that there's a meaningful difference in the likelihood of the read coming from one haplotype compared to
     * its next best haplotype.
     *
     * The best way to visualize the output of this mode is with IGV.  Tell IGV to color the alignments by tag,
     * and give it the HC tag, so you can see which reads support each haplotype.  Finally, you can tell IGV
     * to group by sample, which will separate the potential haplotypes from the reads.  All of this can be seen
     * in the following screenshot: https://www.dropbox.com/s/xvy7sbxpf13x5bp/haplotypecaller%20bamout%20for%20docs.png
     *
     */
    @Advanced
    @Output(fullName="bamOutput", shortName="bamout", doc="File to which assembled haplotypes should be written", required = false, defaultToStdout = false)
    protected GATKSAMFileWriter bamWriter = null;
    protected HaplotypeBAMWriter haplotypeBAMWriter;

    /**
     * The type of BAM output we want to see.
     */
    @Advanced
    @Argument(fullName="bamWriterType", shortName="bamWriterType", doc="How should haplotypes be written to the BAM?", required = false)
    public HaplotypeBAMWriter.Type bamWriterType = HaplotypeBAMWriter.Type.CALLED_HAPLOTYPES;


// PAIR-HMM-Related Goodness


    /**
     * The phredScaledGlobalReadMismappingRate reflects the average global mismapping rate of all reads, regardless of their
     * mapping quality.  This term effects the probability that a read originated from the reference haplotype, regardless of
     * its edit distance from the reference, in that the read could have originated from the reference haplotype but
     * from another location in the genome.  Suppose a read has many mismatches from the reference, say like 5, but
     * has a very high mapping quality of 60.  Without this parameter, the read would contribute 5 * Q30 evidence
     * in favor of its 5 mismatch haplotype compared to reference, potentially enough to make a call off that single
     * read for all of these events.  With this parameter set to Q30, though, the maximum evidence against the reference
     * that this (and any) read could contribute against reference is Q30.
     *
     * Set this term to any negative number to turn off the global mapping rate
     */
    @Advanced
    @Argument(fullName="phredScaledGlobalReadMismappingRate", shortName="globalMAPQ", doc="The global assumed mismapping rate for reads", required = false)
    protected int phredScaledGlobalReadMismappingRate = 45;

    @Hidden
    @Argument(fullName="noFpga", shortName="noFpga", doc="If provided, disables the use of the FPGA HMM implementation", required = false)
    protected boolean noFpga = false;

//    public PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.CONSERVATIVE;
//    public PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.AGGRESSIVE;
    public PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.HOSTILE;

    /**
     * The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime.
     */
    @Hidden
    @Argument(fullName = "pair_hmm_implementation", shortName = "pairHMM", doc = "The PairHMM implementation to use for genotype likelihood calculations", required = false)
    public PairHMM.HMM_IMPLEMENTATION pairHMM = PairHMM.HMM_IMPLEMENTATION.VECTOR_LOGLESS_CACHING;

    @Advanced
    @Argument(fullName="gcpHMM", shortName="gcpHMM", doc="Flat gap continuation penalty for use in the Pair HMM", required = false)
    protected int gcpHMM = 10;


    /**
     * This argument is intended for use in the test suite only. It gives developers the ability to select of the
     * hardware dependent vectorized implementation of the vectorized PairHMM library (pairHMM=VECTOR_LOGLESS_CACHING).
     * For normal usage, you should rely on the architecture auto-detection.
     */
    @Hidden
    @Advanced
    @Argument(fullName = "pair_hmm_sub_implementation", shortName = "pairHMMSub", doc = "The PairHMM machine-dependent sub-implementation to use for genotype likelihood calculations", required = false)
    public PairHMM.HMM_SUB_IMPLEMENTATION pairHMMSub = PairHMM.HMM_SUB_IMPLEMENTATION.ENABLE_ALL;

    /**
     * This argument is intended for use in the test suite only. It gives developers the ability to load different
     * hardware dependent sub-implementations (-pairHMMSub) of the vectorized PairHMM library (-pairHMM=VECTOR_LOGLESS_CACHING)
     * for each test. Without this option, the library is only loaded once (for the first test executed in the suite) even if
     * subsequent tests specify a different implementation.
     * Each test will output the corresponding library loading messages.
     */
    @Hidden
    @Advanced
    @Argument(fullName = "always_load_vector_logless_PairHMM_lib", shortName = "alwaysloadVectorHMM", doc = "Load the vector logless PairHMM library each time a GATK run is initiated in the test suite", required = false)
    public boolean alwaysLoadVectorLoglessPairHMMLib = false;


    // CODE TO RUN LOCAL ASSEMBLY -- SHOULD BE PACKAGED OUT OF GATK!!!

    // Parameters to control read error correction
    @Hidden
    @Argument(fullName="errorCorrectReads", shortName="errorCorrectReads", doc = "Use an exploratory algorithm to error correct the kmers used during assembly.  May cause fundamental problems with the assembly graph itself", required=false)
    protected boolean errorCorrectReads = false;

    @Hidden
    @Argument(fullName="kmerLengthForReadErrorCorrection", shortName="kmerLengthForReadErrorCorrection", doc = "Use an exploratory algorithm to error correct the kmers used during assembly.  May cause fundamental problems with the assembly graph itself", required=false)
    protected int kmerLengthForReadErrorCorrection = 25;

    @Hidden
    @Argument(fullName="minObservationsForKmerToBeSolid", shortName="minObservationsForKmerToBeSolid", doc = "A k-mer must be seen at least these times for it considered to be solid", required=false)
    protected int minObservationsForKmerToBeSolid = 20;

    @Hidden
    @Argument(fullName="captureAssemblyFailureBAM", shortName="captureAssemblyFailureBAM", doc="If specified, we will write a BAM called assemblyFailure.bam capturing all of the reads that were in the active region when the assembler failed for any reason", required = false)
    protected boolean captureAssemblyFailureBAM = false;

    @Advanced
    @Argument(fullName="dontUseSoftClippedBases", shortName="dontUseSoftClippedBases", doc="If specified, we will not analyze soft clipped bases in the reads", required = false)
    protected boolean dontUseSoftClippedBases = false;


    @Hidden
    @Argument(fullName="justDetermineActiveRegions", shortName="justDetermineActiveRegions", doc = "If specified, the HC won't actually do any assembly or calling, it'll just run the upfront active region determination code.  Useful for benchmarking and scalability testing", required=false)
    protected boolean justDetermineActiveRegions = false;




    // reference base padding size
    private static final int REFERENCE_PADDING = 500;

    private static final byte MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION = 6;
    private final static int maxReadsInRegionPerSample = 1000; // TODO -- should be an argument
    private final static int minReadsPerAlignmentStart = 5; // TODO -- should be an argument



    /**
     * High-level function that runs the assembler on the active region reads,
     * returning a data structure with the resulting information needed
     * for further HC steps
     *
     * @param activeRegion the region we should assemble
     * @param giveAlleles additional alleles we might need to genotype (can be empty)
     * @return the AssemblyResult describing how to proceed with genotyping
     */
    protected AssemblyResultSet assembleReads(final ActiveRegion activeRegion, final List<VariantContext> giveAlleles) {
        // Create the reference haplotype which is the bases from the reference that make up the active region
        finalizeActiveRegion(activeRegion); // handle overlapping fragments, clip adapter and low qual tails

        final byte[] fullReferenceWithPadding = activeRegion.getActiveRegionReference(referenceReader, REFERENCE_PADDING);
        final GenomeLoc paddedReferenceLoc = getPaddedLoc(activeRegion);
        final Haplotype referenceHaplotype = createReferenceHaplotype(activeRegion, paddedReferenceLoc);

        // Create ReadErrorCorrector object if requested - will be used within assembly engine.
        ReadErrorCorrector readErrorCorrector = null;
        if (errorCorrectReads)
            readErrorCorrector = new ReadErrorCorrector(kmerLengthForReadErrorCorrection, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION, minObservationsForKmerToBeSolid, SCAC.DEBUG, fullReferenceWithPadding);

        try {
            final AssemblyResultSet assemblyResultSet = assemblyEngine.runLocalAssembly( activeRegion, referenceHaplotype, fullReferenceWithPadding, paddedReferenceLoc, giveAlleles,readErrorCorrector );
            assemblyResultSet.debugDump(logger);
            return assemblyResultSet;

        } catch ( final Exception e ) {
            // Capture any exception that might be thrown, and write out the assembly failure BAM if requested
            if ( captureAssemblyFailureBAM ) {
                final SAMFileWriter writer = SAMFileWriterStub.createSAMFileWriter("assemblyFailure.bam", getToolkit());
                new DirectOutputTracker().addOutput((SAMFileWriterStub) writer);
                for ( final GATKSAMRecord read : activeRegion.getReads() ) {
                    writer.addAlignment(read);
                }
                writer.close();
            }
            throw e;
        }
    }

    private void finalizeActiveRegion( final ActiveRegion activeRegion ) {
        if (activeRegion.isFinalized()) return;

        if( SCAC.DEBUG ) { logger.info("Assembling " + activeRegion.getLocation() + " with " + activeRegion.size() + " reads:    (with overlap region = " + activeRegion.getExtendedLoc() + ")"); }

        // Loop through the reads hard clipping the adaptor and low quality tails
        final List<GATKSAMRecord> readsToUse = new ArrayList<>(activeRegion.getReads().size());
        for( final GATKSAMRecord myRead : activeRegion.getReads() ) {
            GATKSAMRecord clippedRead;
            if (errorCorrectReads)
                clippedRead = ReadClipper.hardClipLowQualEnds( myRead, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION );
            else  // default case: clip low qual ends of reads
                clippedRead= ReadClipper.hardClipLowQualEnds( myRead, MIN_TAIL_QUALITY );

            if ( dontUseSoftClippedBases || ! ReadUtils.hasWellDefinedFragmentSize(clippedRead) ) {
                // remove soft clips if we cannot reliably clip off adapter sequence or if the user doesn't want to use soft clips at all
                clippedRead = ReadClipper.hardClipSoftClippedBases(clippedRead);
            } else {
                // revert soft clips so that we see the alignment start and end assuming the soft clips are all matches
                // TODO -- WARNING -- still possibility that unclipping the soft clips will introduce bases that aren't
                // TODO -- truly in the extended region, as the unclipped bases might actually include a deletion
                // TODO -- w.r.t. the reference.  What really needs to happen is that kmers that occur before the
                // TODO -- reference haplotype start must be removed
                clippedRead = ReadClipper.revertSoftClippedBases(clippedRead);
            }

            clippedRead = ( clippedRead.getReadUnmappedFlag() ? clippedRead : ReadClipper.hardClipAdaptorSequence( clippedRead ) );
            if( !clippedRead.isEmpty() && clippedRead.getCigar().getReadLength() > 0 ) {
                clippedRead = ReadClipper.hardClipToRegion(clippedRead, activeRegion.getExtendedLoc().getStart(), activeRegion.getExtendedLoc().getStop());
                if( activeRegion.readOverlapsRegion(clippedRead) && clippedRead.getReadLength() > 0 ) {
                    //logger.info("Keeping read " + clippedRead + " start " + clippedRead.getAlignmentStart() + " end " + clippedRead.getAlignmentEnd());
                    readsToUse.add(clippedRead);
                }
            }
        }

        // TODO -- Performance optimization: we partition the reads by sample 4 times right now; let's unify that code.

        final List<GATKSAMRecord> downsampledReads = DownsamplingUtils.levelCoverageByPosition(ReadUtils.sortReadsByCoordinate(readsToUse), maxReadsInRegionPerSample, minReadsPerAlignmentStart);

        // handle overlapping read pairs from the same fragment
        // KC: commented out as we handle overlapping read pairs in a different way...
        //cleanOverlappingReadPairs(downsampledReads, normalSampleNames);

        activeRegion.clearReads();
        activeRegion.addAll(downsampledReads);
        activeRegion.setFinalized(true);
    }

    private GenomeLoc getPaddedLoc( final ActiveRegion activeRegion ) {
        final int padLeft = Math.max(activeRegion.getExtendedLoc().getStart()-REFERENCE_PADDING, 1);
        final int padRight = Math.min(activeRegion.getExtendedLoc().getStop()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(activeRegion.getExtendedLoc().getContig()).getSequenceLength());
        return getToolkit().getGenomeLocParser().createGenomeLoc(activeRegion.getExtendedLoc().getContig(), padLeft, padRight);
    }

    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     * @param activeRegion the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the GenomeLoc which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    private Haplotype createReferenceHaplotype(final ActiveRegion activeRegion, final GenomeLoc paddedReferenceLoc) {
        return ReferenceConfidenceModel.createReferenceHaplotype(activeRegion, activeRegion.getActiveRegionReference(referenceReader), paddedReferenceLoc);
    }

    /**
     * Clean up reads/bases that overlap within read pairs
     *
     * @param reads the list of reads to consider
     */
    private void cleanOverlappingReadPairs(final List<GATKSAMRecord> reads, Set<String> normalSampleNames) {
        Map<String, List<GATKSAMRecord>> data = splitReadsBySample(reads);
        for ( String sampleName : data.keySet() ) {
            final boolean isTumor = !normalSampleNames.contains(sampleName);
            final List<GATKSAMRecord> perSampleReadList = data.get(sampleName);

            final FragmentCollection<GATKSAMRecord> fragmentCollection = FragmentUtils.create(perSampleReadList);
            for ( final List<GATKSAMRecord> overlappingPair : fragmentCollection.getOverlappingPairs() )

                // in MuTect -- right now we compare the
                FragmentUtils.adjustQualsOfOverlappingPairedFragments(overlappingPair);


        }
    }

    public static void logReadInfo(String readName, Collection<GATKSAMRecord> records, String message) {
        if (readName != null) {
            for (GATKSAMRecord rec : records) {
                logReadInfo(readName, rec, message);
            }

        }
    }

    public static void logReadInfo(String readName, GATKSAMRecord rec, String message) {
        if (readName != null && rec != null && readName.equals(rec.getReadName())) {
            logger.info("Found " + rec.toString() + " - " + message);
        }
    }

    /**
     * Returns a map with the original read as a key and the realigned read as the value.
     * <p>
     *     Missing keys or equivalent key and value pairs mean that the read was not realigned.
     * </p>
     * @return never {@code null}
     */
    // TODO: migrate from HC -> HCUtils Class and share it!
    private Map<GATKSAMRecord,GATKSAMRecord> realignReadsToTheirBestHaplotype(final ReadLikelihoods<Haplotype> originalReadLikelihoods, final Haplotype refHaplotype, final GenomeLoc paddedReferenceLoc) {

        final Collection<ReadLikelihoods<Haplotype>.BestAllele> bestAlleles = originalReadLikelihoods.bestAlleles();
        final Map<GATKSAMRecord,GATKSAMRecord> result = new HashMap<>(bestAlleles.size());

        for (final ReadLikelihoods<Haplotype>.BestAllele bestAllele : bestAlleles) {
            final GATKSAMRecord originalRead = bestAllele.read;
            final Haplotype bestHaplotype = bestAllele.allele;
            final boolean isInformative = bestAllele.isInformative();
            final GATKSAMRecord realignedRead = AlignmentUtils.createReadAlignedToRef(originalRead, bestHaplotype, refHaplotype, paddedReferenceLoc.getStart(), isInformative);
            result.put(originalRead,realignedRead);
        }
        return result;
    }

    private boolean isReadFromNormal(GATKSAMRecord rec) {
        return normalSampleName != null && normalSampleName.equals(rec.getReadGroup().getSample());

    }
    // KCIBUL: new stuff -- read up on this!!
    /**
     * As of GATK 3.3, HaplotypeCaller outputs physical (read-based) information (see version 3.3 release notes and documentation for details). This argument disables that behavior.
     */
    @Advanced
    @Argument(fullName="doNotRunPhysicalPhasing", shortName="doNotRunPhysicalPhasing", doc="Disable physical phasing", required = false)
    protected boolean doNotRunPhysicalPhasing = false;

    public static final String HAPLOTYPE_CALLER_PHASING_ID_KEY = "PID";
    public static final String HAPLOTYPE_CALLER_PHASING_GT_KEY = "PGT";
}


