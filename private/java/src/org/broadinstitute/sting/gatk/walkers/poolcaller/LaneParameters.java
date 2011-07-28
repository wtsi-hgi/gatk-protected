package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.Collection;

/**
 * A support class to facilitate future addition/removal of parameters to the Lane class
 *
 * @author Mauricio Carneiro
 * @since 7/27/11
 */
public class LaneParameters {
    public String name;
    public ReadBackedPileup lanePileup;
    public String referenceSampleName;
    public Collection<Byte> trueReferenceBases;
    public byte referenceSequenceBase;
    public byte minQualityScore;
    public byte maxQualityScore;
    public byte phredScaledPrior;
    public int maxAlleleCount;
    public double minCallQual;
    public double minPower;

    public LaneParameters(String name, ReadBackedPileup lanePileup, String referenceSampleName, Collection<Byte> trueReferenceBases, byte referenceSequenceBase, byte minQualityScore, byte maxQualityScore, byte phredScaledPrior, int maxAlleleCount, double minCallQual, double minPower) {
        this.name = name;
        this.lanePileup = lanePileup;
        this.referenceSampleName = referenceSampleName;
        this.trueReferenceBases = trueReferenceBases;
        this.referenceSequenceBase = referenceSequenceBase;
        this.minQualityScore = minQualityScore;
        this.maxQualityScore = maxQualityScore;
        this.phredScaledPrior = phredScaledPrior;
        this.maxAlleleCount = maxAlleleCount;
        this.minCallQual = minCallQual;
        this.minPower = minPower;
    }
}
