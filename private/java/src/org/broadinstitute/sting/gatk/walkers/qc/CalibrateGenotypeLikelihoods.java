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

package org.broadinstitute.sting.gatk.walkers.qc;

import net.sf.samtools.SAMReadGroupRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.R.RScriptExecutor;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.io.Resource;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Computes raw GL calibration data for read groups in BAMs against a comp VCF track of genotypes
 *
 * <p>
 *  This walker generates plots that display the empirical accuracy of the likelihoods of the genotypes calculated in a
 *  given BAM file. To calculate the empirical error, this walker needs a truth ROD (e.g. OMNI or a gold standard callset)
 *  and it is <i>essential</i> that the truth callset includes all possible genotypes (AA, AB, BB).
 * </p>
 *
 *
 * <h3>Input</h3>
 * <p>
 *  This walker takes two inputs:
 *  <ul>
 *      <li>The BAM file to genotype (given truth alleles) and evaluate the accuracy</li>
 *      <li>A truth ROD with the gold standard genotypes</li>
 *  </ul>
 *  <b>Warning: The truth ROD MUST include all possible genotypes to build the error model appropriately (AA, AB, BB)</b>
 * </p>
 *
 * <h3>Output</h3>
 *  <p>
 *      Two intermediate tables:
 *      <ul>
 *          <li>Raw output of all sites with the three genotypes and their likelihoods</li>
 *          <li>A digested table stratified by read group and platforms (for plotting)</li>
 *      </ul>
 *
 *      Three plots:
 *      <ul>
 *          <li>The calibration curve of the three possible genotypes (AA, AB, BB)</li>
 *          <li>The calibration curve stratified by read group</li>
 *          <li>The calibration curve stratified by platform</li>
 *      </ul>
 *  <p>
 *
 * <h3>Example</h3>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -T CalibrateGenotypeLikelihoods \
 *   -R reference/human_g1k_v37.fasta \
 *   -alleles omni.vcf \
 *   -I testBAM.bam \
 *   -o calibration
 * </pre>
 *
 *
 * @author depristo
 * @since May, 2011
 */

@Requires(value={DataSource.REFERENCE})
@Allows(value={DataSource.READS, DataSource.REFERENCE})

// Ugly fix because RodWalkers don't have access to reads
@By(DataSource.REFERENCE)
@Reference(window=@Window(start=-200,stop=200))
public class CalibrateGenotypeLikelihoods extends RodWalker<CalibrateGenotypeLikelihoods.Data, CalibrateGenotypeLikelihoods.Data> implements TreeReducible<CalibrateGenotypeLikelihoods.Data> {

    @Input(fullName="alleles", shortName = "alleles", doc="The set of alleles at which to genotype when in GENOTYPE_MODE = GENOTYPE_GIVEN_ALLELES", required=false)
    public RodBinding<VariantContext> alleles;

    @Input(fullName="externalLikelihoods",shortName="el",doc="Any number of VCFs with external likelihoods for which to evaluate their calibration.",required=false)
    public List<RodBinding<VariantContext>> externalLikelihoods = Collections.emptyList();

    @Argument(fullName="minimum_base_quality_score", shortName="mbq", doc="Minimum base quality score for calling a genotype", required=false)
    private int mbq = -1;

    @Argument(fullName="maximum_deletion_fraction", shortName="deletions", doc="Maximum deletion fraction for calling a genotype", required=false)
    private double deletions = -1;
    
    @Argument(fullName="indels", shortName="indels", doc="Do indel evaluation", required=false)
    private boolean doIndels = false;

    @Argument(fullName="repeats", shortName="repeats", doc="Do indel evaluation", required=false)
    private boolean doRepeats = false;

    @Argument(fullName="skipFilteredRecords", shortName="skipFiltered", doc="Skip filtered records when evaluating external likelihoods", required=false)
    private boolean skipFilteredRecords = false;

    //@Argument(fullName="standard_min_confidence_threshold_for_calling", shortName="stand_call_conf", doc="the minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls", required=false)
    private double callConf = 0;

    @Output(doc="The name of the output files for both tables and pdf (name will be prepended to the appropriate extensions)")
    private File moltenDatasetFileName;

    PrintStream moltenDataset;
    Set<String> samples;
    String SCRIPT_FILE = "CalibrateGenotypeLikelihoods.R";

    /**
     * Trivial wrapper class.  Data is a collection of Datum.
     */
    public static class Data {
        final Collection<Datum> values;

        public Data() { this(new LinkedList<Datum>()); }
        public Data(Collection<Datum> data) { values = data; }

        final public static Data EMPTY_DATA = new Data(Collections.<Datum>emptyList());
    }

    /**
     * The raw datapoints we are tracking for a specific site for a specific sample.
     * read group id and sample name.  The PL object.
     * the ref and alt alleles. The type of the variant context.  And the genotype of the
     * comp. track at this site.
     */
    public static class Datum implements Comparable<Datum> {
        final String rgID, sample;
        final GenomeLoc loc;
        final GenotypeLikelihoods pl;
        final String ref, alt;
        final VariantContext.Type siteType;
        final GenotypeType genotypeType;
        final boolean isRepeat;

        @Override
        public int compareTo(Datum o) {
            int bySample = sample.compareTo(o.sample);
            int byRG = rgID.compareTo(o.rgID);
            return bySample != 0 ? bySample : byRG;
        }

        public Datum(final GenomeLoc loc, String ref, String alt, String sample, String rgID, GenotypeLikelihoods pl, VariantContext.Type siteType, GenotypeType genotypeType, boolean isrep) {
            this.loc = loc;
            this.ref = ref;
            this.alt = alt;
            this.sample = sample;
            this.rgID = rgID;
            this.pl = pl;
            this.siteType = siteType;
            this.genotypeType = genotypeType;
            this.isRepeat = isrep;
        }
    }

    private UnifiedGenotyperEngine snpEngine;
    private UnifiedGenotyperEngine indelEngine;

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        // because either can be specified, check that exactly one was
        if ( externalLikelihoods.isEmpty() && ! alleles.isBound()  ) {
            throw new UserException("No input VCF was given. Please provide either both external likelihoods and alleles VCFs, or alleles VCF and read data.");
        } else if ( !externalLikelihoods.isEmpty() && ! alleles.isBound() ) {
            throw new UserException("Alleles VCF not provided to compare against external likelihoods.");
        }

        try {
            moltenDataset = new PrintStream(moltenDatasetFileName);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(moltenDatasetFileName, e);
        }
        
        // Get the samples from the VCF file
        Set<String> vcfSamples = null;
        List<ReferenceOrderedDataSource> rods = getToolkit().getRodDataSources();
        for (ReferenceOrderedDataSource rod : rods) {
            VCFHeader header =  (VCFHeader) rod.getHeader();
            vcfSamples = new HashSet<String>();
            vcfSamples.addAll(header.getGenotypeSamples());
        }

        // Get the samples from the BAM file
        if( externalLikelihoods.isEmpty() ) {
            samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
            if (samples.isEmpty())
                throw new UserException.BadInput("Bam file has no samples in the ReadGroup tag");
            logger.info("Samples in your BAM file: " + samples);
        } else {
            samples = vcfSamples;
        }

        // Assert that the user provided a VCF (not some other tribble format)
        if (vcfSamples == null)
            throw new UserException.BadInput("Must provide a VCF. This walker is designed to work with VCFs only");

        // Assert that the VCFs have samples and genotypes
        if (vcfSamples.isEmpty())
            throw new UserException.MalformedVCFHeader("The VCF must have samples so we can find the correct genotype for the samples in the BAM file");

        // This is the list of samples represented in both the BAM and the VCF
        vcfSamples.retainAll(samples);

        // An empty list means the VCF doesn't have any samples from the BAM file
        if (vcfSamples.isEmpty())
            throw new UserException.BadInput("The VCF file does not contain any of the samples in your BAM file");

        // Uneven sets means there are some samples missing from the BAM in the VCF -- warn the user but don't throw exception
        if (vcfSamples.size() < samples.size()) {
            Set<String> samplesPresent = new HashSet<String>();
            samplesPresent.addAll(samples);
            samplesPresent.removeAll(vcfSamples);
            logger.warn("VCF file does not contain ALL samples in your BAM file.  Missing samples: " + Utils.join(", ", samplesPresent));
        }

        // Filling in SNP calling arguments for UG
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES;
        uac.GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
        if (mbq >= 0) uac.MIN_BASE_QUALTY_SCORE = mbq;
        if (deletions >= 0) uac.MAX_DELETION_FRACTION = deletions;
        uac.STANDARD_CONFIDENCE_FOR_CALLING = callConf;
        uac.CONTAMINATION_FRACTION = 0.0;
        uac.alleles = alleles;
        // Adding the INDEL calling arguments for UG
        if (doIndels)  {
            uac.GLmodel = GenotypeLikelihoodsCalculationModel.Model.INDEL;
            indelEngine = new UnifiedGenotyperEngine(getToolkit(), uac, Logger.getLogger(UnifiedGenotyperEngine.class), null, null, samples, 2 * samples.size() );
        }
        else {
            uac.GLmodel = GenotypeLikelihoodsCalculationModel.Model.SNP;
            snpEngine = new UnifiedGenotyperEngine(getToolkit(), uac, Logger.getLogger(UnifiedGenotyperEngine.class), null, null, samples, 2 * samples.size() );

        }

        if (doRepeats)
            SCRIPT_FILE = "CalibrateGenotypeLikelihoodsByRepeat.R";
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public Data map( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {
        if ( tracker == null || tracker.getNTracksWithBoundFeatures() == 0 )
            return Data.EMPTY_DATA;

        // Grabs a usable VariantContext from the Alleles ROD
        final VariantContext vcComp = getVCComp(tracker, ref, context);
        if( vcComp == null )
            return Data.EMPTY_DATA;

        return ( externalLikelihoods.isEmpty() ? calculateGenotypeDataFromAlignments(tracker, ref, context, vcComp) : calculateGenotypeDataFromExternalVC(tracker, ref, context, vcComp) );
    }

    /**
     * Convenience function that determines the genotype in the comp VC for sample
     *
     * @param tracker
     * @param ref
     * @param sample
     * @return
     */
    private static Genotype getGenotype( final RefMetaDataTracker tracker, final ReferenceContext ref, final String sample, final UnifiedGenotyperEngine engine ) {
        for ( final VariantContext vc : tracker.getValues(engine.getUAC().alleles, ref.getLocus()) ) {
            if ( vc.isNotFiltered() && vc.hasGenotype(sample) )
                return vc.getGenotype(sample);
            else
                return null;
        }

        return null;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Data reduceInit() {
        return new Data();
    }

    @Override
    public Data treeReduce( final Data sum1, final Data sum2 ) {
        sum2.values.addAll(sum1.values);
        return sum2;
    }

    @Override
    public Data reduce( final Data mapValue, final Data reduceSum ) {
        return treeReduce(mapValue, reduceSum);
    }

    @Override
    public void onTraversalDone( final Data data ) {
        // print the header
        final List<String> pGNames = Arrays.asList("QofAAGivenD", "QofABGivenD", "QofBBGivenD");
        final List<String> fields = Arrays.asList("sample", "rg", "loc", "ref", "alt", "siteType", "pls", "comp", "pGGivenDType", "pGGivenD", "pDGivenG");
        if (doRepeats)
            fields.add("isRepeat");
        moltenDataset.println(Utils.join("\t", fields));

        // determine the priors by counting all of the events we've seen in comp
        final double[] counts = new double[]{1, 1, 1};
        for ( final Datum d : data.values ) { counts[d.genotypeType.ordinal()-1]++; }
        double sum = MathUtils.sum(counts);
        logger.info(String.format("Types %s %s %s", GenotypeType.values()[1], GenotypeType.values()[2], GenotypeType.values()[3]));
        logger.info(String.format("Counts %.0f %.0f %.0f %.0f", counts[0], counts[1], counts[2], sum));
        double[] log10priors = new double[]{Math.log10(counts[0] / sum), Math.log10(counts[1] / sum), Math.log10(counts[2] / sum)};
        logger.info(String.format("Priors %.2f %.2f %.2f", log10priors[0], log10priors[1], log10priors[2]));

        // emit the molten data set
        for ( final Datum d : data.values ) {
            final double[] log10pDGivenG = d.pl.getAsVector();
            final double[] log10pGGivenD = log10pDGivenG.clone();
            for ( int i = 0; i < log10priors.length; i++ ) log10pGGivenD[i] += log10priors[i];
            final double[] pOfDGivenG = MathUtils.normalizeFromLog10(log10pDGivenG, false);
            final double[] pOfGGivenD = MathUtils.normalizeFromLog10(log10pGGivenD, false);
            for ( int i = 0; i < pGNames.size(); i++ ) {
                final int qDGivenG = QualityUtils.trueProbToQual(pOfDGivenG[i]);
                final int qGGivenD = QualityUtils.trueProbToQual(pOfGGivenD[i]);
                if ( qGGivenD > 1 ) { // tons of 1s, and not interesting
                    if (doRepeats)
                    moltenDataset.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%b%n",
                            d.sample, d.rgID, d.loc, d.ref, d.alt, d.siteType, d.pl.getAsString(), d.genotypeType.toString(),
                            pGNames.get(i), qGGivenD, qDGivenG, d.isRepeat);
                    else
                        moltenDataset.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d%n",
                                d.sample, d.rgID, d.loc, d.ref, d.alt, d.siteType, d.pl.getAsString(), d.genotypeType.toString(),
                                pGNames.get(i), qGGivenD, qDGivenG);

                }
            }
        }

        moltenDataset.close();
        RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(SCRIPT_FILE, CalibrateGenotypeLikelihoods.class));
        executor.addArgs(moltenDatasetFileName.getAbsolutePath());
        executor.addArgs( externalLikelihoods.isEmpty() ? 1 : 0 ); // only plot the "ALL" read group if running without external likelihood files
        if (doRepeats)
            executor.addArgs(1); // stratify by repeat class if in indel mode
        logger.info("Generating plots...");
        executor.exec();
    }

    private VariantContext getVCComp( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {
        if (doIndels)
            return UnifiedGenotyperEngine.getVCFromAllelesRod( tracker, ref, context.getLocation(), false, logger, indelEngine.getUAC().alleles );
        return UnifiedGenotyperEngine.getVCFromAllelesRod( tracker, ref, context.getLocation(), false, logger, snpEngine.getUAC().alleles );
    }

    private Data calculateGenotypeDataFromAlignments( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context, final VariantContext vcComp ) {
        final Data data = new Data();
        final Map <String,AlignmentContext> contextBySample = AlignmentContextUtils.splitContextBySampleName(context);
        for ( final Map.Entry<String,AlignmentContext> sAC : contextBySample.entrySet())  {
            String sample = sAC.getKey();
            AlignmentContext sampleAC = sAC.getValue();
            Genotype compGT = getGenotype(tracker, ref, sample, snpEngine != null ? snpEngine : indelEngine);
            if ( compGT == null || compGT.isNoCall() )
                continue;

            // now split by read group
            Map<SAMReadGroupRecord,AlignmentContext> byRG = AlignmentContextUtils.splitContextByReadGroup(sampleAC, getToolkit().getSAMFileHeader().getReadGroups());
            byRG.put(new SAMReadGroupRecord("ALL"), context);     // uncomment to include a synthetic RG for all RG for the sample
            for ( final Map.Entry<SAMReadGroupRecord, AlignmentContext> rgAC : byRG.entrySet() ) {
                VariantCallContext call;
                if ( (vcComp.isIndel() || vcComp.isMixed()) && doIndels ) {
                    //throw new UserException.BadInput("CalibrateGenotypeLikelihoods does not currently support indel GL calibration.  This capability needs to be tested and verified to be working with the new genotyping code for indels in UG");
                    call = indelEngine.calculateLikelihoodsAndGenotypes(tracker, ref, rgAC.getValue()).get(0);
                } else if (vcComp.isSNP() && !doIndels) {
                    call = snpEngine.calculateLikelihoodsAndGenotypes(tracker, ref, rgAC.getValue()).get(0);
                }
                else
                    break;

                if ( call == null )
                    throw new ReviewedStingException("Unexpected genotyping failure " + sample + " at " + ref.getLocus());

                final Genotype rgGT = call.getGenotype(sample);

                if ( rgGT != null && ! rgGT.isNoCall() && rgGT.hasLikelihoods() ) {
                    addValue(data, vcComp, ref, sample, rgAC.getKey().getReadGroupId(), rgGT, compGT);
                }
            }
        }

        return data;
    }

    private final void addValue( final Data data, final VariantContext vcComp, final ReferenceContext ref, final String sample,
                                 final String rgID, final Genotype calledGT, final Genotype compGT) {
        String refs, alts;
        boolean isRepeat = false;
        if (vcComp.isIndel()) {
            refs = vcComp.getReference().getDisplayString();
            alts = vcComp.getAlternateAllele(0).getDisplayString();
            isRepeat = (GATKVariantContextUtils.isTandemRepeat(vcComp, ref.getForwardBases()));
        } else {
            refs = vcComp.getReference().getBaseString();
            alts = vcComp.getAlternateAllele(0).getBaseString();
        }
        final GenomeLoc loc = getToolkit().getGenomeLocParser().createGenomeLoc(vcComp);
        final Datum d = new Datum(loc, refs, alts, sample, rgID, calledGT.getLikelihoods(),
                vcComp.getType(), compGT.getType(), isRepeat);
        data.values.add(d);
    }

    private Data calculateGenotypeDataFromExternalVC( final RefMetaDataTracker tracker,
                                                      final ReferenceContext ref,
                                                      final AlignmentContext context,
                                                      final VariantContext vcComp ) {

        // the tracker should contain a VCF with external likelihoods.
        final Data data = new Data();
        for( final RodBinding<VariantContext> rod : externalLikelihoods ) {
            final VariantContext extVC = tracker.getFirstValue(rod);
            if ( extVC == null ) {
                return Data.EMPTY_DATA;
            }

            if ( extVC.isFiltered() && skipFilteredRecords )
                return Data.EMPTY_DATA; // skip filtered eval records

            // make sure there is an alternate allele and that it matches exactly the extVC allele
            if ( vcComp.getAlternateAlleles() == null || vcComp.getAlternateAlleles().size() == 0 || !vcComp.hasSameAllelesAs(extVC) ) {
                return Data.EMPTY_DATA;
            }
    
            for ( final Genotype genotype : extVC.getGenotypes() ) {
                final String sample = genotype.getSampleName();
                final Genotype compGT = vcComp.hasGenotype(sample) ? vcComp.getGenotype(sample) : null;
                if ( compGT == null || genotype.isNoCall() || compGT.isNoCall() )
                    continue;

                addValue(data, vcComp, ref, sample, rod.getName() + "." + genotype.getSampleName(), genotype, compGT);
            }
        }

        return data;
    }
}
