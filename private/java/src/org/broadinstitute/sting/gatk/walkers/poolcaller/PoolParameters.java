package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * A support class to facilitate future addition/removal of parameters to the Pool class
 *
 * @author Mauricio Carneiro
 * @since 7/27/11
 */
public class PoolParameters {
    public String name;
    public ReadBackedPileup pileup;
    public ErrorModel errorModel;
    public byte referenceSequenceBase;
    public int maxAlleleCount;
    public double minCallQual;

    public PoolParameters(String name, ReadBackedPileup pileup, ErrorModel errorModel, byte referenceSequenceBase, int maxAlleleCount, double minCallQual) {
        this.name = name;
        this.pileup = pileup;
        this.errorModel = errorModel;
        this.referenceSequenceBase = referenceSequenceBase;
        this.maxAlleleCount = maxAlleleCount;
        this.minCallQual = minCallQual;
    }
}
