package org.broadinstitute.sting.gatk.walkers.replication_validation;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.MathUtils;
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

    @Argument(shortName="maxqs", fullName="max_quality_score", doc="Max quality score to consider. Smaller numbers process faster. Default: Q40.", required=false)
    int maxQualityScore= 40;

    @Argument(shortName="prior", fullName="site_quality_prior", doc="Phred-Scaled prior quality of the site. Default: Q20.", required=false)
    byte defaultPrior= 20;

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
    double THETA = 0.001; // Human heterozygozity rate

    final String REFERENCE_ROD_NAME = "reference";
    final String TRUTH_ROD_NAME = "truth";


    /**
     * GATK Engine creates readgroups of the form XXX.Y.Z
     * XXX.Y is the unique lane identifier.
     *     Z is the id of the sample to make the read group id unique
     * This function returns the list of lane identifiers.
     *
     * @param readGroups readGroups A collection of read group strings (obtained from the alignment context pileup)
     * @return a collection of lane ids.
     */
    private Set<String> getLaneIDs(Collection<String> readGroups) {
        HashSet<String> result = new HashSet<String>();
        for (String rgid : readGroups) {
            result.add(getLaneID(rgid));
        }
        return result;
    }

    /**
     * GATK Engine creates readgroups of the form XXX.Y.Z
     * XXX.Y is the unique lane identifier.
     *     Z is the id of the sample to make the read group id unique
     *
     * @param readGroupID the read group id string
     * @return just the lane id (the XXX.Y string)
     */
    private String getLaneID (String readGroupID) {
        String [] parsedID = readGroupID.split("\\.");
        return parsedID[0] + "." + parsedID[1];
    }

    /**
     * Calculates he probability of the data (reference sample reads) given the phred scaled site quality score.
     * @param referenceSamplePileup reference sample pileup
     * @param refBases base from the reference sequence for this site
     * @param phredScaledPrior phred scaled expected site quality (prior)
     * @return an array of log10 probabilities of site qualities ranging from Q1-Q40.
     */
    private double[] buildErrorModel (ReadBackedPileup referenceSamplePileup, Collection<Byte> refBases, byte phredScaledPrior) {
        double [] model = new double[maxQualityScore+1];
        byte [] data = referenceSamplePileup.getBases();

        int coverage = data.length;
        int mismatches = getNumberOfMismatches(data, refBases);
        int matches = coverage - mismatches;

        for (byte q=0; q<=maxQualityScore; q++) {
            double probMismatch = MathUtils.phredScaleToProbability(q);
            model[q] = MathUtils.phredScaleToLog10Probability(phredScaledPrior) +
                       MathUtils.log10BinomialCoefficient(coverage, mismatches) +
                       mismatches * Math.log10(probMismatch) +
                       matches * Math.log10(1-probMismatch);
        }
        return model;
    }

    /**
     * Returns the number of mismatching bases in a pileup
     * @param data the bases of a pileup
     * @param refBases the true bases to compare to
     * @return number of bases in data that are different from all refBases
     */
    private int getNumberOfMismatches (byte[] data, Collection<Byte> refBases) {
        int mismatches = 0;
        for (byte seqBase : data) {
            if (!refBases.contains(seqBase))
                mismatches++;
        }
        return mismatches;
    }

    /**
     * Returns the number of mismatching bases in a pileup
     * @param data the bases of a pileup
     * @param refBase the true base to comare to
     * @return number of bases in data that are different from refBase
     */
    private int getNumberOfMismatches (byte[] data, Byte refBase) {
        ArrayList<Byte> refBases = new ArrayList<Byte>(1);
        refBases.add(refBase);
        return getNumberOfMismatches(data, refBases);
    }


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

    /**
     * The prior probability of an allele count being observed based solely on the human heterozygozity rate
     * and the number of samples
     *
     * @param ac given allele count
     * @return the prior probability of ac
     */
    private double log10PriorAC(int ac) {
        // prior probability for a given allele count is THETA/AC.
        if (ac > 0)
            return Math.log10(THETA/ac);

        // if allele count is 0, the prior is one minus the sum of all other possibilities.
        double result = 0.0;
        for (int i=1; i<=maxAlleleCount; i++)
            result += THETA/i;
        return Math.log10(1-result);
    }

    /**
     * Calculates the pool's probability for all possible allele counts. Calculation is based on the error model
     * generated by the reference sample on the same lane. The probability is given by :
     *
     * Pr(ac=j | pool, errorModel) = sum_over_all_Qs ( Pr(ac=j) * Pr(errorModel_q) * [ (n-j/2n) * (1-e) + (j/n)*e]^m * [(n-j/n)*e + (j/n) * (1-e)]^(1-m)
     *
     * where:
     *  n = number of chromosomes
     *  e = probability of an error at a given Q level (e.g. Q30 = 0.001, Q20 = 0.01, ...)
     *  m = number of mismatches
     *
     * @param pool
     * @param errorModel
     * @param refBase
     * @return
     */
    private double[] getPoolACProbabilityDistribution (ReadBackedPileup pool, double[] errorModel, byte refBase) {
        double [] result = new double[maxAlleleCount+1];
        int mismatches = getNumberOfMismatches(pool.getBases(), refBase);
        int matches = pool.size() - mismatches;

        for (int ac=0; ac<=maxAlleleCount; ac++) {
            result[ac] = log10pOfAC(ac, errorModel, matches, mismatches);
        }

        return result;
    }

    @Requires({"ac>=0", "errorModel != null", "matches >= 0", "mismatches >= 0" })
    @Ensures({"result <= 0", "! Double.isInfinite(result)", "! Double.isNaN(result)"})
    private double log10pOfAC(int ac, double[] errorModel, int matches, int mismatches) {
        double log10PAC = log10PriorAC(ac);

        // for each quality probability in the model, calculate the probability of the allele count = ac
        // we skip Q0 because it's meaningless.
        double [] acc = new double[maxQualityScore]; // we're skipping Q0 so we don't need maxQualityScore + 1 here.
        for (int i = 0; i < maxQualityScore; i++) {
            final int qual = i + 1;
            acc[i] = log10pOfACGivenQual(ac, qual, errorModel, matches, mismatches);
        }
        return log10PAC + MathUtils.log10sumLog10(acc);
    }

    @Requires({"ac>=0", "qual>0", "errorModel != null && qual < errorModel.length", "matches >= 0", "mismatches >= 0" })
    @Ensures({"result <= 0", "! Double.isInfinite(result)", "! Double.isNaN(result)"})
    private double log10pOfACGivenQual(int ac, int qual, double[] errorModel, int matches, int mismatches) {
        double p = (double) ac / maxAlleleCount;
        double q = 1 - p;
        double e = MathUtils.phredScaleToProbability((byte)qual);
        double x = Math.log10(q * (1-e) + p * e);
        double y = Math.log10(q * e + p * (1-e));
        return errorModel[qual] + matches * x + mismatches * y;
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

        ReadBackedPileup contextPileup = context.getBasePileup();
        Set<String> lanesInLocus = getLaneIDs(contextPileup.getReadGroups());
        for (String laneID : lanesInLocus) {

            // make a pileup for this lane
            ReadBackedPileup lanePileup;
            if (DEBUG_IGNORE_LANES) lanePileup = contextPileup;
            else lanePileup = contextPileup.getPileupForLane(laneID);

            Collection<String> samplesInLane = lanePileup.getSampleNames();

            // we can only analyze loci that have reads for the reference sample
            if (samplesInLane.contains(referenceSampleName)) {

                // build reference sample pileup
                ReadBackedPileup referenceSamplePileup = lanePileup.getPileupForSampleName(referenceSampleName);

                // Build error model
                double [] errorModel = buildErrorModel(referenceSamplePileup, trueReferenceBases, defaultPrior);

                // iterate over all samples (pools) in this lane except the reference
                samplesInLane.remove(referenceSampleName);
                for (String pool : samplesInLane) {

                    // get the pileup for the pool
                    ReadBackedPileup poolPileup;
                    if (DEBUG_IGNORE_LANES) poolPileup = lanePileup;
                    else poolPileup = lanePileup.getPileupForSampleName(pool);

                    double [] AC = getPoolACProbabilityDistribution(poolPileup, errorModel, ref.getBase());

                    // Debug AC distribution
                    System.out.println("\n\n" + "[" + ref.getLocus() + "] " + laneID + " - " + pool +
                            "\nNumber of Samples: " + nSamples +
                            "\nRefSample Size: " + referenceSamplePileup.getBases().length +
                            "\nRefSample MMs: " + getNumberOfMismatches(referenceSamplePileup.getBases(), trueReferenceBases) +
                            "\nRefQ: " + MathUtils.maxElementIndex(errorModel) +
                            "\nPool Size: " + poolPileup.size() +
                            "\nPool MMs: " + getNumberOfMismatches(poolPileup.getBases(), ref.getBase()) +
                            "\nPool AF: " + (double) getNumberOfMismatches(poolPileup.getBases(), ref.getBase())/poolPileup.size() +
                            "\nTrutn AN: " + truthContext.getAttribute("AN") +
                            "\nPool AC / Truth AC: " + MathUtils.maxElementIndex(AC) + " / " + truthContext.getAttribute("AC") + " / " + truthContext.isFiltered());

                    System.out.println("\nError Model: ");
                    for (double v : errorModel)
                        System.out.print(v + ", ");
                    System.out.println("\n");

                    System.out.println("AC Distribution: ");
                    for (double v : AC)
                        System.out.print(v + ", ");
                    System.out.println();

                    if (DEBUG_IGNORE_LANES) break;
                }
            }

            if (DEBUG_IGNORE_LANES) break;

            // todo: merge pools
            // todo: decide whether or not it's a variant
        }
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

