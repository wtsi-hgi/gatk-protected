package org.broadinstitute.sting.gatk.walkers.poolcaller;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.picard.util.MathUtil;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 7/21/11
 * Time: 2:21 PM
 *
 * This is a site based implementation of an Error Model. The error model is a probability
 * distribution for the site given the phred scaled quality.
 */
public class ErrorModel extends ProbabilityModel {
    private byte maxQualityScore;
    private byte minQualityScore;
    private byte phredScaledPrior;
    private double log10minPower;
    private int refDepth;
    private boolean hasData = false;
    /**
     * Calculates the probability of the data (reference sample reads) given the phred scaled site quality score.
     */
    public ErrorModel (byte minQualityScore, byte maxQualityScore, byte phredScaledPrior, ReferenceSample referenceSample, double minPower) {
        this.maxQualityScore = maxQualityScore;
        this.minQualityScore = minQualityScore;
        this.phredScaledPrior = phredScaledPrior;
        log10minPower = Math.log10(minPower);


        model = new double[maxQualityScore-minQualityScore+1];

        if (referenceSample.getPileup() == null /*|| referenceSample.getPileup().isEmpty()*/ ) {
            double p = MathUtils.phredScaleToLog10Probability((byte)(maxQualityScore-minQualityScore));
            for (byte q=minQualityScore; q<=maxQualityScore; q++) {
                int i = q - minQualityScore; // fill the array from 0 to (maxQualityscore - minQualityScore)
                // maximum uncertainty if there's no ref data at site
                model[i] = p;
            }
            this.refDepth = 0;
        }
        else {
            byte [] data = referenceSample.getPileup().getBases();
            hasData = true;
            int coverage = data.length;
            int matches = 0;
            for (Allele allele : referenceSample.getTrueAlleles()) {
                byte base = allele.getBases()[0];
                // todo: tmp, need to generalize for indels
                matches += MathUtils.countOccurrences(base, data);
            }
            int mismatches = coverage - matches;

            for (byte q=minQualityScore; q<=maxQualityScore; q++) {
                int i = q - minQualityScore; // fill the array from 0 to (maxQualityscore - minQualityScore)
                //model[i] = log10ProbabilitySiteGivenQual(q, coverage, matches, mismatches);
                model[i] = log10PoissonProbabilitySiteGivenQual(q,coverage, matches, mismatches);
            }
            this.refDepth = coverage;
        }
    }

    @Requires({
            "q >= minQualityScore",
            "q <= maxQualityScore",
            "coverage >= 0",
            "matches >= 0",
            "matches <= coverage",
            "mismatches >= 0",
            "mismatches <= coverage"
    })
    @Ensures({"result <= 0", "! Double.isInfinite(result)", "! Double.isNaN(result)"})
//todo -- memoize this function
    // returns log10(p^mismatches* (1-p)^matches * choose(coverage,mismatches) * prior)

    private double log10ProbabilitySiteGivenQual(byte q, int coverage, int matches, int mismatches) {
        double probMismatch = MathUtils.phredScaleToProbability(q);
        return MathUtils.phredScaleToLog10Probability(phredScaledPrior) +
                MathUtils.log10BinomialCoefficient(coverage, mismatches) +
                mismatches * Math.log10(probMismatch) +
                matches * Math.log10(1-probMismatch);
    }

    private double log10PoissonProbabilitySiteGivenQual(byte q, int coverage, int matches, int mismatches) {
        // same as   log10ProbabilitySiteGivenQual but with Poisson approximation to avoid numerical underflows
        double lambda = MathUtils.phredScaleToProbability(q) * (double )coverage;
        // log(e^-lambda*lambda^k/k!) = -lambda + k*log(lambda) - logfactorial(k)
        return Math.log10(lambda)*mismatches - lambda- MathUtils.log10Factorial(mismatches);
    }

    @Requires({"qual-minQualityScore <= maxQualityScore"})
    public double getErrorGivenQual (int qual) {
        int index = qual - minQualityScore;
        return model[index];
    }

    public byte getMaxQualityScore() {
        return maxQualityScore;
    }

    public byte getMinQualityScore() {
        return minQualityScore;
    }

    public int getReferenceDepth() {
        return refDepth;
    }
    public boolean hasData() {
        return hasData;
    }

    public double[] getErrorModelVector() {
        return model;
    }
@Requires({"maxAlleleCount >= 0"})
//todo -- memoize this function
    public boolean hasPowerForMaxAC (int maxAlleleCount) {
        int siteQ = (int) Math.ceil(MathUtils.probabilityToPhredScale((double) 1/maxAlleleCount));
        double log10CumSum = getCumulativeSum(siteQ);
        return log10CumSum < log10minPower;
    }
}
