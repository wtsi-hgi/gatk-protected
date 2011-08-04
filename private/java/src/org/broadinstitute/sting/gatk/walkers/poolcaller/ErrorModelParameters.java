package org.broadinstitute.sting.gatk.walkers.poolcaller;

/**
 * A support class to facilitate future addition/removal of parameters to the Error Model class
 *
 * @author Mauricio Carneiro
 * @since 7/27/11
 */
public class ErrorModelParameters {
    public byte maxQualityScore;
    public byte minQualityScore;
    public byte phredScaledPrior;
    public ReferenceSample referenceSample;
    public double minPower;

    public ErrorModelParameters(byte minQualityScore, byte maxQualityScore, byte phredScaledPrior, ReferenceSample referenceSample, double minPower) {
        this.maxQualityScore = maxQualityScore;
        this.minQualityScore = minQualityScore;
        this.phredScaledPrior = phredScaledPrior;
        this.referenceSample = referenceSample;
        this.minPower = minPower;
    }
}
