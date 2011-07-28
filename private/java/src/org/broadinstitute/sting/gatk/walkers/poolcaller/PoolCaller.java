package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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
 *      -T PoolCaller
 *      -R reference.fasta \
 *	    -I mySequences.bam \
 *	    -B:reference,VCF /humgen/gsa-hpprojects/dev/carneiro/repval/data/reference_sample.vcf \
 *	    -refsample NA12878
 *  </pre>
 *
 * @author Mauricio Carneiro
 * @since 5/4/11
 */

public class PoolCaller extends LocusWalker<Integer, Long> implements TreeReducible<Long> {

    @Output(doc="File to which variants should be written", required=true)
    protected VCFWriter vcfWriter = null;

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

    @Argument(shortName = "min_call_conf", fullName = "min_confidence_threshold_for_calling", doc="The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be called.", required = false)
    double minCallQual = 30.0;

    @Argument(shortName = "min_call_power", fullName = "min_power_threshold_for_calling", doc="The minimum confidence in the error model to make a call. Number should be between 0 (no power requirement) and 1 (maximum power required).", required = false)
    double minPower = 0.95;

    @Argument(shortName="ef", fullName="exclude_filtered_reference_sites", doc="Don't include in the analysis sites where the reference sample VCF is filtered. Default: false.", required=false)
    boolean EXCLUDE_FILTERED_REFERENCE_SITES = false;

    @Hidden
    @Argument(shortName = "dl", doc="DEBUG ARGUMENT -- treats all reads as coming from the same lane", required=false)
    boolean DEBUG_IGNORE_LANES = false;


    int nSamples;
    int maxAlleleCount;
    boolean USE_TRUTH_ROD;

    final String REFERENCE_ROD_NAME = "reference";
    final String TRUTH_ROD_NAME = "truth";


    /**
     * Returns the true bases for the reference sample in this locus. Homozygous loci will return one base
     * but heterozygous will return two bases (hence why it returns a collection).
     *
     * @param referenceSampleContext the variant context from the reference sample ROD track
     * @param ref the reference sequence context
     * @return the true bases for the reference sample.
     */
    private Collection<Byte> getTrueBases(VariantContext referenceSampleContext, ReferenceContext ref) {

        ArrayList<Byte> trueReferenceBase = new ArrayList<Byte>();

        // Site is not a variant, take from the reference
        if (referenceSampleContext == null) {
            trueReferenceBase.add(ref.getBase());
        }

        else if (referenceSampleContext.isIndel()) {
            return null; // TODO: add special treatment for extended events. For Now just skip these altogether.
        }

        // Site has a VCF entry -- is variant
        else {
            // Site is filtered, don't mess with it if option is set
            if (referenceSampleContext.isFiltered() && EXCLUDE_FILTERED_REFERENCE_SITES) {
                return null;
            }

            Genotype referenceGenotype = referenceSampleContext.getGenotype(referenceSampleName);
            List<Allele> referenceAlleles = referenceGenotype.getAlleles();
            for (Allele allele : referenceAlleles) {
                byte [] bases = allele.getBases();
                for (byte b : bases) {
                    trueReferenceBase.add(b);
                }
            }
        }
        return trueReferenceBase;
    }


    public void initialize() {

        // Set the number of samples in the pools ( - reference sample)
        nSamples = getToolkit().getSAMFileSamples().size() - 1;

        // Set the max allele count (defines the size of the error model array)
        maxAlleleCount = 2*nSamplesPerPool;

        // Look for the reference ROD and the optional truth ROD. If truth is provided, set the truth "test" mode ON.
        List<ReferenceOrderedDataSource> rods = getToolkit().getRodDataSources();
        if (rods.size() < 1) {
            throw new IllegalArgumentException("You must provide a reference ROD.");
        }
        boolean foundReferenceROD = false;
        boolean foundTruthROD = false;
        for (ReferenceOrderedDataSource rod : rods) {
            if (rod.getName().equals(REFERENCE_ROD_NAME)) {
                foundReferenceROD = true;
            }
            if (rod.getName().equals(TRUTH_ROD_NAME)) {
                foundTruthROD = true;
            }
        }
        if (!foundReferenceROD) {
            throw new IllegalArgumentException("You haven't provided a reference ROD. Note that the reference ROD must be labeled " + REFERENCE_ROD_NAME + ".");
        }
        USE_TRUTH_ROD = foundTruthROD;

        // Initialize the VCF
        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.add(new VCFInfoHeaderLine("AC", 1, VCFHeaderLineType.Integer, "Allele count in the site, number of alternate alleles across all pools"));
        headerLines.add(new VCFInfoHeaderLine("AF", 1, VCFHeaderLineType.Float, "Allele frequency in the site. Proportion of the alternate alleles across all pools"));
        headerLines.add(new VCFInfoHeaderLine("AN", 1, VCFHeaderLineType.Integer, "Total number of alleles in the site. Total number of chromosomes represented on this site across all pools"));
        headerLines.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Total depth in the site. Sum of the depth of all pools"));
        headerLines.add(new VCFInfoHeaderLine("MQ", 1, VCFHeaderLineType.Float, "RMS mapping quality of all reads in the site"));
        headerLines.add(new VCFInfoHeaderLine("MQ0", 1, VCFHeaderLineType.Integer, "Total number of mapping quality zero reads in the site"));
        headerLines.add(new VCFFormatHeaderLine("AD", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Allelic depths for the ref and alt alleles in the order listed"));
        headerLines.add(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Read Depth (only filtered reads used for calling)"));
        headerLines.add(new VCFFormatHeaderLine("GQ", 1, VCFHeaderLineType.Float, "Genotype Quality"));
        headerLines.add(new VCFFormatHeaderLine("AL", 3, VCFHeaderLineType.Integer, "Allele count likelihood and the 5% confidence interval"));
        headerLines.add(new VCFFilterHeaderLine("lowQual", "Low quality"));
        headerLines.add(new VCFFilterHeaderLine("lowConf", "Low confidence"));
        headerLines.add(new VCFHeaderLine("PoolCaller", "todo -- add parameter list here"));
        vcfWriter.writeHeader(new VCFHeader(headerLines, SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader())));
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        // Get reference base from VCF or Reference
        VariantContext referenceSampleContext = tracker.getVariantContext(ref, REFERENCE_ROD_NAME, context.getLocation());
        VariantContext truthContext = tracker.getVariantContext(ref, TRUTH_ROD_NAME, context.getLocation());
        Collection<Byte> trueReferenceBases = getTrueBases(referenceSampleContext, ref);

        // If there is no true reference base in this locus, skip it.
        if (trueReferenceBases == null)
            return 0;

        // Creates a site Object for this location
        SiteParameters siteParameters = new SiteParameters(context.getBasePileup(),
                                                referenceSampleName,
                                                trueReferenceBases,
                                                ref.getBase(),
                                                minQualityScore,
                                                maxQualityScore,
                                                phredScaledPrior,
                                                maxAlleleCount,
                                                minCallQual,
                                                minPower);
        Site site = DEBUG_IGNORE_LANES ? Site.debugSite(siteParameters) : new Site(siteParameters);

        VariantContext call = new VariantContext("PoolCaller",
                                                  ref.getLocus().getContig(),
                                                  ref.getLocus().getStart(),
                                                  ref.getLocus().getStop(),
                                                  site.getAlleles(),
                                                  site.getGenotypes(),
                                                  site.getNegLog10PError(),
                                                  site.getFilters(),
                                                  site.getAttributes());


        vcfWriter.add(call, ref.getBase());
        return 1;
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

