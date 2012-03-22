package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 3/13/12
 * Time: 2:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class PoolCallerUnifiedArgumentCollection extends UnifiedArgumentCollection {

    /**
     * Controls the model used to calculate the probability that a site is variant plus the various sample genotypes in the data at a given locus.
     */
 //   @Argument(fullName = "genotype_likelihoods_model", shortName = "glm", doc = "Genotype likelihoods calculation model to employ -- SNP is the default option, while INDEL is also available for calling indels and BOTH is available for calling both together", required = false)
    //@Override
    //public String GLmodel = PoolGenotypeLikelihoodsCalculationModel.Model.POOLSNP.name();

    /**
     * Controls the model used to calculate the probability that a site is variant plus the various sample genotypes in the data at a given locus.
     */
//    @Argument(fullName = "p_nonref_model", shortName = "pnrm", doc = "Non-reference probability calculation model to employ -- EXACT is the default option, while GRID_SEARCH is also available.", required = false)
    //public AlleleFrequencyCalculationModel.Model AFmodel = AlleleFrequencyCalculationModel.Model.EXACT;

    @Argument(fullName = "allReadsSP", shortName = "dl", doc = "expt", required = false)
    public boolean TREAT_ALL_READS_AS_SINGLE_POOL = false;

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

    public PoolCallerUnifiedArgumentCollection clone() {
        PoolCallerUnifiedArgumentCollection uac = (PoolCallerUnifiedArgumentCollection)super.clone();

        uac.GLmodel = GLmodel;
        uac.TREAT_ALL_READS_AS_SINGLE_POOL = TREAT_ALL_READS_AS_SINGLE_POOL;
        uac.referenceSampleRod = referenceSampleRod;
        uac.referenceSampleName = referenceSampleName;
        uac.nSamplesPerPool = nSamplesPerPool;
        uac.maxQualityScore = minQualityScore;
        uac.maxQualityScore = maxQualityScore;
        uac.phredScaledPrior = phredScaledPrior;
        uac.minPower = minPower;
        uac.minReferenceDepth = minReferenceDepth;
        uac.EXCLUDE_FILTERED_REFERENCE_SITES = EXCLUDE_FILTERED_REFERENCE_SITES;
        return uac;

    }
}
