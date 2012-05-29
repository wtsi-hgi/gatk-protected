/*
 * Copyright (c) 2010 The Broad Institute
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
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.io.Resource;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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
 * <h2>Input</h2>
 * <p>
 *  This walker takes two inputs:
 *  <ul>
 *      <li>The BAM file to genotype (given truth alleles) and evaluate the accuracy</li>
 *      <li>A truth ROD with the gold standard genotypes</li>
 *  </ul>
 *  <b>Warning: The truth ROD MUST include all possible genotypes to build the error model appropriately (AA, AB, BB)</b>
 * </p>
 *
 * <h2>Output</h2>
 *  <p>
 *      Two intermediate tables:
 *      <ul>
 *          <li>Raw ouput of all sites with the three genotypes and their likelihoods</li>
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
 * <h2>Example</h2>
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
    
    @Argument(fullName="noBanded", shortName="noBanded", doc="No Banded indel GL computation", required=false)
    private boolean noBandedIndelGLs = false;

    //@Argument(fullName="standard_min_confidence_threshold_for_calling", shortName="stand_call_conf", doc="the minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls", required=false)
    private double callConf = 0;

    @Output(doc="The name of the output files for both tables and pdf (name will be prepended to the appropriate extensions)", required=true)
    private File moltenDatasetFileName;

    PrintStream moltenDataset;
    Set<String> samples;
    final String SCRIPT_FILE = "CalibrateGenotypeLikelihoods.R";

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
        final Genotype.Type genotypeType;

        @Override
        public int compareTo(Datum o) {
            int bySample = sample.compareTo(o.sample);
            int byRG = rgID.compareTo(o.rgID);
            return bySample != 0 ? bySample : byRG;
        }

        public Datum(final GenomeLoc loc, String ref, String alt, String sample, String rgID, GenotypeLikelihoods pl, VariantContext.Type siteType, Genotype.Type genotypeType) {
            this.loc = loc;
            this.ref = ref;
            this.alt = alt;
            this.sample = sample;
            this.rgID = rgID;
            this.pl = pl;
            this.siteType = siteType;
            this.genotypeType = genotypeType;
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
        uac.alleles = alleles;
        // Adding the INDEL calling arguments for UG
        if (doIndels)  {
            uac.GLmodel = GenotypeLikelihoodsCalculationModel.Model.INDEL;
            uac.DONT_DO_BANDED_INDEL_COMPUTATION = noBandedIndelGLs;
            indelEngine = new UnifiedGenotyperEngine(getToolkit(), uac, Logger.getLogger(UnifiedGenotyperEngine.class), null, null, samples,2*samples.size() );
        }
        else {
            uac.GLmodel = GenotypeLikelihoodsCalculationModel.Model.SNP;
            snpEngine = new UnifiedGenotyperEngine(getToolkit(), uac, Logger.getLogger(UnifiedGenotyperEngine.class), null, null, samples,2*samples.size() );

        }
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
        final VariantContext vcComp = getVCComp(tracker,ref,context);
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
        List<String> pGNames = Arrays.asList("QofAAGivenD", "QofABGivenD", "QofBBGivenD");
        List<String> fields = Arrays.asList("sample", "rg", "loc", "ref", "alt", "siteType", "pls", "comp", "pGGivenDType", "pGGivenD", "pDGivenG");
        moltenDataset.println(Utils.join("\t", fields));

        // determine the priors by counting all of the events we've seen in comp
        final double[] counts = new double[]{1, 1, 1};
        for ( final Datum d : data.values ) { counts[d.genotypeType.ordinal()-1]++; }
        double sum = MathUtils.sum(counts);
        logger.info(String.format("Types %s %s %s", Genotype.Type.values()[1], Genotype.Type.values()[2], Genotype.Type.values()[3]));
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
                final int qDGivenG = QualityUtils.probToQual(pOfDGivenG[i], QualityUtils.ERROR_RATE_OF_MAX_QUAL_SCORE);
                final int qGGivenD = QualityUtils.probToQual(pOfGGivenD[i], QualityUtils.ERROR_RATE_OF_MAX_QUAL_SCORE);
                if ( qGGivenD > 1 ) { // tons of 1s, and not interesting
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
        Map <String,AlignmentContext> contextBySample = AlignmentContextUtils.splitContextBySampleName(context);
        for (Map.Entry<String,AlignmentContext> sAC : contextBySample.entrySet()) {
            String sample = sAC.getKey();
            AlignmentContext sampleAC = sAC.getValue();
            Genotype compGT = getGenotype(tracker, ref, sample, snpEngine != null ? snpEngine : indelEngine);
            if ( compGT == null || compGT.isNoCall() )
                continue;

            // now split by read group
            Map<SAMReadGroupRecord,AlignmentContext> byRG = AlignmentContextUtils.splitContextByReadGroup(sampleAC, getToolkit().getSAMFileHeader().getReadGroups());
            byRG.put(new SAMReadGroupRecord("ALL"), context);     // uncomment to include a synthetic RG for all RG for the sample
            for ( Map.Entry<SAMReadGroupRecord, AlignmentContext> rgAC : byRG.entrySet() ) {
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
                    throw new ReviewedStingException("Unexpected genotyping failure " + sample + " at " + ref.getLocus() + " call " + call);

                final Genotype rgGT = call.getGenotype(sample);

                if ( rgGT != null && ! rgGT.isNoCall() && rgGT.hasLikelihoods() ) {
                    addValue(data, vcComp, ref, sample, rgAC.getKey().getReadGroupId(), rgGT, compGT);
                }
            }
        }

        return data;
    }

    private final void addValue( final Data data, final VariantContext vcComp, final ReferenceContext ref, final String sample,
                                 final String rgID, final Genotype calledGT, final Genotype compGT ) {
        String refs, alts;
        if (vcComp.isIndel()) {
            refs = ((char)ref.getBase()) + vcComp.getReference().getDisplayString();
            alts = ((char)ref.getBase()) + vcComp.getAlternateAllele(0).getDisplayString();
        } else {
            refs = vcComp.getReference().getBaseString();
            alts = vcComp.getAlternateAllele(0).getBaseString();
        }
        final GenomeLoc loc = getToolkit().getGenomeLocParser().createGenomeLoc(vcComp);
        final Datum d = new Datum(loc, refs, alts, sample, rgID, calledGT.getLikelihoods(),
                vcComp.getType(), compGT.getType());
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
    
            // make sure there is an alternate allele
            if ( vcComp.getAlternateAlleles() == null || vcComp.getAlternateAlleles().size() == 0 ) {
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
