package org.broadinstitute.sting.gatk.walkers.replication_validation;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.*;

/**
 * Implementation of the replication and validation framework with reference based error model
 * for pooled sequencing.
 *
 * The input should be a BAM file with pooled sequencing data where each pool is represented by
 * samples with the same barcode.
 *
 * A reference sample name must be provided and it must be barcoded uniquely.
 */
public class ReplicationValidationWalker extends LocusWalker<Integer, Long> implements TreeReducible<Long> {

    @Argument(shortName="refsample", fullName="reference_sample_name", doc="Reference sample name.", required=true)
    String referenceSampleName;

    @Argument(shortName="nchr", fullName="number_of_chromosomes", doc="Number of chromosomes per sample (in case you're not dealing with diploids). Default: 2.", required=false)
    int ploidy = 2;

    @Argument(shortName="maxac", fullName="max_allele_count", doc="Max number of alleles expected in a site. Smaller numbers process faster. Default: 2 * number of samples. ", required=false)
    int overrideMaxAlleleCount = -1;

    @Argument(shortName="minqs", fullName="min_quality_score", doc="Min quality score to consider. Smaller numbers process faster. Default: Q1.", required=false)
    byte minQualityScore= 1;

    @Argument(shortName="maxqs", fullName="max_quality_score", doc="Max quality score to consider. Smaller numbers process faster. Default: Q40.", required=false)
    byte maxQualityScore= 40;

    @Argument(shortName="prior", fullName="site_quality_prior", doc="Phred-Scaled prior quality of the site. Default: Q20.", required=false)
    byte phredScaledPrior = 20;

    @Argument(shortName="ef", fullName="exclude_filtered_reference_sites", doc="Don't include in the analysis sites where the reference sample VCF is filtered. Default: false.", required=false)
    boolean EXCLUDE_FILTERED_REFERENCE_SITES = false;

    @Hidden
    @Argument(shortName = "dl", doc="DEBUG ARGUMENT -- treats all reads as coming from the same lane", required=false)
    boolean DEBUG_IGNORE_LANES = false;

    @Output(doc="Write output to this file instead of STDOUT")
    PrintStream out;

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

        // If we ignore the lanes, then the Reference Sample will be included in the analysis
        if (DEBUG_IGNORE_LANES)
            nSamples++;

        // Set the max allele count (defines the size of the error model array)
        maxAlleleCount = (overrideMaxAlleleCount > 0) ? overrideMaxAlleleCount : nSamples*ploidy;

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
        Site site = new Site(context.getBasePileup(), referenceSampleName, trueReferenceBases, ref.getBase(), minQualityScore, maxQualityScore, phredScaledPrior, maxAlleleCount);

        // todo -- get site's alleleCountModel and write out VCF

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
        out.println(c);
    }
}

