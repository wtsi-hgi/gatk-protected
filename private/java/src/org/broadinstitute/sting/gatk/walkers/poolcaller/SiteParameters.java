package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.Collection;

/**
 * A support class to facilitate future addition/removal of parameters to the Site class
 *
 * @author Mauricio Carneiro
 * @since 7/27/11
 */
public class SiteParameters {
    public ReadBackedPileup sitePileup;
    public String referenceSampleName;
    public Collection<Byte> trueReferenceBases;
    public byte referenceSequenceBase;
    public byte minQualityScore;
    public byte maxQualityScore;
    public byte phredScaledPrior;
    public int maxAlelleCount;
    public double minCallQual;
    public double minPower;

    public SiteParameters(ReadBackedPileup sitePileup, String referenceSampleName, Collection<Byte> trueReferenceBases, byte referenceSequenceBase, byte minQualityScore, byte maxQualityScore, byte phredScaledPrior, int maxAlelleCount, double minCallQual, double minPower) {
        this.sitePileup = sitePileup;
        this.referenceSampleName = referenceSampleName;
        this.trueReferenceBases = trueReferenceBases;
        this.referenceSequenceBase = referenceSequenceBase;
        this.minQualityScore = minQualityScore;
        this.maxQualityScore = maxQualityScore;
        this.phredScaledPrior = phredScaledPrior;
        this.maxAlelleCount = maxAlelleCount;
        this.minCallQual = minCallQual;
        this.minPower = minPower;
    }
}
