package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

/**
 * Implementation of the replication and validation framework with reference based error model
 * for pooled sequencing.
 *
 *  <p>
 *     [Long description of the walker]
 *  </p>
 *
 *
 * <h2>Input</h2>
 *  <p>
 *      The input should be a BAM file with pooled sequencing data where each pool is represented by
 *      samples with the same barcode. A reference sample name must be provided and it must be barcoded
 *      uniquely.
 *  </p>
 *
 * <h2>Output</h2>
 *  <p>
 *      A VCF file with the following annotations:
 *      <ul>
 *          <li>Likelihood that the site is not AC=0</li>
 *          <li>site AC</li>
 *          <li>each pool AC, 95% confidence interval for AC, likelihood</li>
 *      </ul>
 *  </p>
 *
 * <h2>Examples</h2>
 *  <pre>
 *    java
 *      -javaagent:/home/unix/carneiro/src/gatk/repval/lib/cofoja-1.0-20110609.jar \
 *      -jar GenomeAnalysisTK.jar
 *      -T PoolCallerWalker
 *      -R reference.fasta \
 *	    -I mySequences.bam \
 *	    -B:reference,VCF /humgen/gsa-hpprojects/dev/carneiro/repval/data/reference_sample.vcf \
 *	    -refsample NA12878
 *  </pre>
 *
 * @author Mauricio Carneiro
 * @since 5/4/11
 */

@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
@ReadFilters( {BadMateFilter.class, MappingQualityUnavailableFilter.class} )
@Reference(window=@Window(start=-200,stop=200))
@By(DataSource.REFERENCE)
@Downsample(by= DownsampleType.BY_SAMPLE, toCoverage=250)
public class OldPoolCaller extends LocusWalker<Integer, Long> implements TreeReducible<Long> {

    @Output(doc="File to which variants should be written", required=true)
    protected VCFWriter vcfWriter = null;

    @Input(fullName="reference_sample", shortName = "reference", doc="VCF file with the truth callset for the reference sample", required=true)
    RodBinding<VariantContext> referenceSampleRod;

    @Argument(shortName="refsample", fullName="reference_sample_name", doc="Reference sample name.", required=true)
    String referenceSampleName;

    @Argument(shortName="sp", fullName="samples_per_pool", doc="Number of samples in each pool (must be the same for all pools).", required=true)
    int nSamplesPerPool;

    @Argument(shortName="minqs", fullName="min_quality_score", doc="Min quality score to consider. Smaller numbers process faster. Default: Q1.", required=false)
    byte minQualityScore= 1;

    @Argument(shortName="maxqs", fullName="max_quality_score", doc="Max quality score to consider. Smaller numbers process faster. Default: Q40.", required=false)
    byte maxQualityScore= 40;

    @Argument(shortName="prior", fullName="site_quality_prior", doc="Phred-Scaled prior quality of the site. Default: Q20.", required=false)
    byte phredScaledPrior = 20;

     @Argument(shortName = "min_call_power", fullName = "min_power_threshold_for_calling", doc="The minimum confidence in the error model to make a call. Number should be between 0 (no power requirement) and 1 (maximum power required).", required = false)
    double minPower = 0.95;

    @Argument(shortName = "min_depth", fullName = "min_reference_depth", doc="The minimum depth required in the reference sample in order to make a call.", required = false)
    int minReferenceDepth = 100;

    @Argument(shortName="ef", fullName="exclude_filtered_reference_sites", doc="Don't include in the analysis sites where the reference sample VCF is filtered. Default: false.", required=false)
    boolean EXCLUDE_FILTERED_REFERENCE_SITES = false;

    @Hidden
    @Argument(shortName = "dl", doc="DEBUG ARGUMENT -- treats all reads as coming from the same lane", required=false)
    boolean DEBUG_IGNORE_LANES = false;

    // the unified argument collection
    @ArgumentCollection
     private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();


    private Set<String> samples;
    int nSamples;
    int maxAlleleCount;

    // enable deletions in the pileup
    public boolean includeReadsWithDeletionAtLoci() { return true; }


    /**
     * Returns the true bases for the reference sample in this locus. Homozygous loci will return one base
     * but heterozygous will return two bases (hence why it returns a collection).
     *
     * @param referenceSampleContext the variant context from the reference sample ROD track
     * @param ref the reference sequence context
     * @return the true bases for the reference sample.
     */
    private Collection<Allele> getTrueAlleles(VariantContext referenceSampleContext, ReferenceContext ref) {

        ArrayList<Allele> trueReferenceAlleles = new ArrayList<Allele>();

        // Site is not a variant, take from the reference
        if (referenceSampleContext == null) {
            trueReferenceAlleles.add(Allele.create(ref.getBase(),true));
        }
        // Site has a VCF entry -- is variant
        else {
            Genotype referenceGenotype = referenceSampleContext.getGenotype(referenceSampleName);
            List<Allele> referenceAlleles = referenceGenotype.getAlleles();
            for (Allele allele : referenceAlleles) {
                 if (!trueReferenceAlleles.contains(allele))
                    trueReferenceAlleles.add(allele);
            }
        }
        return trueReferenceAlleles;
    }


    public void initialize() {

        // Set the max allele count (defines the size of the error model array)
        maxAlleleCount = 2*nSamplesPerPool;

        // Initialize the VCF
        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.add(new VCFInfoHeaderLine("AC", 1, VCFHeaderLineType.Integer, "Allele count in the site, number of alternate alleles across all pools"));
        headerLines.add(new VCFInfoHeaderLine("AF", 1, VCFHeaderLineType.Float, "Allele frequency in the site. Proportion of the alternate alleles across all pools"));
        headerLines.add(new VCFInfoHeaderLine("AN", 1, VCFHeaderLineType.Integer, "Total number of alleles in the site. Total number of chromosomes represented on this site across all pools"));
        headerLines.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Total depth in the site. Sum of the depth of all pools"));
        headerLines.add(new VCFInfoHeaderLine("MQ", 1, VCFHeaderLineType.Float, "RMS mapping quality of all reads in the site"));
        headerLines.add(new VCFInfoHeaderLine("MQ0", 1, VCFHeaderLineType.Integer, "Total number of mapping quality zero reads in the site"));
        headerLines.add(new VCFInfoHeaderLine("RD", 1, VCFHeaderLineType.Integer, "Depth of coverage for reference sample at site"));
        headerLines.add(new VCFFormatHeaderLine("AD", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Allelic depths for the ref and alt alleles in the order listed"));
        headerLines.add(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Read Depth (only filtered reads used for calling)"));
        headerLines.add(new VCFFormatHeaderLine("GQ", 1, VCFHeaderLineType.Float, "Genotype Quality"));
        headerLines.add(new VCFFormatHeaderLine("AL", 3, VCFHeaderLineType.Integer, "Allele count likelihood and the 5% confidence interval"));
        headerLines.add(new VCFInfoHeaderLine("Dels", 1, VCFHeaderLineType.Float, "Fraction of Reads Containing Spanning Deletions"));        headerLines.add(new VCFHeaderLine("PoolCallerWalker", "todo -- add parameter list here"));
        headerLines.add(new VCFFilterHeaderLine(Filters.LOW_QUAL.toString(), "Low quality"));
        headerLines.add(new VCFFilterHeaderLine(Filters.LOW_POWER.toString(), "Low confidence"));
        headerLines.add(new VCFFilterHeaderLine(Filters.LOW_REFERENCE_SAMPLE_DEPTH.toString(), "Not enough reference sample depth"));
        headerLines.add(new VCFFilterHeaderLine(Filters.NO_BASES_IN_REFERENCE_SAMPLE.toString(), "Reference sample is not covered at all"));
        headerLines.add(new VCFFilterHeaderLine(Filters.FILTERED_REFERENCE_SAMPLE_CALL.toString(), "Reference sample vcf was filtered"));
        headerLines.add(new VCFFilterHeaderLine(Filters.MAX_DELETION_FRACTION_EXCEEDED.toString(), "Site has more deletion fraction than threshold"));

        samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        samples.remove(referenceSampleName);
        // Set the number of samples in the pools ( - reference sample)
        nSamples = samples.size();

        if (DEBUG_IGNORE_LANES) {
            samples.clear();
            samples.add("Pool1");
            nSamples = 1;
        }
        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if ( !BaseUtils.isRegularBase(ref.getBase()) )
            return 0;

        // Get reference base from VCF or Reference
        VariantContext referenceSampleContext = tracker.getFirstValue(referenceSampleRod, context.getLocation());

        boolean filteredRefSampleCall = false;
        if (referenceSampleContext != null && referenceSampleContext.filtersWereApplied() && referenceSampleContext.isFiltered())
            filteredRefSampleCall = true;

        Collection<Allele> trueReferenceAlleles = getTrueAlleles(referenceSampleContext, ref);

        // If there is no true reference base in this locus, skip it.

        if (!context.hasBasePileup())
            return 0;

        VariantContext vcInput = UnifiedGenotyperEngine.getVCFromAllelesRod(tracker, ref, context.getLocation(), false, logger, UAC.alleles);

        if (UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES) {
            // see if there is an input call in GGA mode
            // exclude filtered of multiallelic sites
            if (vcInput == null || vcInput.isFiltered() || !vcInput.isBiallelic())
                return 0;

        }
        Site site = new Site(tracker, ref, context, referenceSampleName, trueReferenceAlleles, minQualityScore, maxQualityScore, phredScaledPrior, maxAlleleCount,
                UAC , minPower, minReferenceDepth, filteredRefSampleCall, DEBUG_IGNORE_LANES, logger, vcInput, samples);

        if (site.needToEmitCall()) {
            vcfWriter.add(site.getCallFromSite());
            return 1;
        }
        return 0;
     }

    public Long reduceInit() {
        return 0l;
    }

    public Long reduce(Integer value, Long sum) {
        return value + sum;
    }

    public Long treeReduce(Long lhs, Long rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone( Long c ) {

    }
}

