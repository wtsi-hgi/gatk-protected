package org.broadinstitute.sting.gatk.walkers.poolcaller;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.MathUtils;

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

    /**
     * Calculates the probability of the data (reference sample reads) given the phred scaled site quality score.
     *
     * @param p a support object with all the necessary parameters for the error model. See class for a more detailed description
     */
    public ErrorModel (ErrorModelParameters p) {
        maxQualityScore = p.maxQualityScore;
        minQualityScore = p.minQualityScore;
        phredScaledPrior = p.phredScaledPrior;
        log10minPower = Math.log10(p.minPower);


        model = new double[maxQualityScore-minQualityScore+1];

        byte [] data = p.referenceSample.getPileup().getBases();
        int coverage = data.length;
        int matches = 0;
        for (byte base : p.referenceSample.getTrueBases()) {
            matches += MathUtils.countOccurrences(base, data);
        }
        int mismatches = coverage - matches;

        for (byte q=minQualityScore; q<=maxQualityScore; q++) {
            int i = q - minQualityScore; // fill the array from 0 to (maxQualityscore - minQualityScore)
            model[i] = log10ProbabilitySiteGivenQual(q, coverage, matches, mismatches);
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
    private double log10ProbabilitySiteGivenQual(byte q, int coverage, int matches, int mismatches) {
        double probMismatch = MathUtils.phredScaleToProbability(q);
        return MathUtils.phredScaleToLog10Probability(phredScaledPrior) +
                MathUtils.log10BinomialCoefficient(coverage, mismatches) +
                mismatches * Math.log10(probMismatch) +
                matches * Math.log10(1-probMismatch);
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

@Requires({"maxAlleleCount >= 0"})
//todo -- memoize this function
    public boolean hasPowerForMaxAC (int maxAlleleCount) {
        int siteQ = (int) Math.ceil(MathUtils.probabilityToPhredScale((double) 1/maxAlleleCount));
        double log10CumSum = getCumulativeSum(siteQ);
        return log10CumSum < log10minPower;
    }
}
