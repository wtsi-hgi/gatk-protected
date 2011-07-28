package org.broadinstitute.sting.gatk.walkers.poolcaller;

/**
 * A support class to facilitate future addition/removal of parameters to the Allele Count Model class
 *
 * @author Mauricio Carneiro
 * @since 7/27/11
 */
public class AlleleCountModelParameters {
    public int maxAlleleCount;
    public ErrorModel errorModel;
    public int matches;
    public int mismatches;
    public double minCallQual;

    public AlleleCountModelParameters(int maxAlleleCount, ErrorModel errorModel, int matches, int mismatches, double minCallQual) {
        this.maxAlleleCount = maxAlleleCount;
        this.errorModel = errorModel;
        this.matches = matches;
        this.mismatches = mismatches;
        this.minCallQual = minCallQual;
    }
}
