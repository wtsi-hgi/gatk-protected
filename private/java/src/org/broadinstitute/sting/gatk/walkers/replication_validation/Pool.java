package org.broadinstitute.sting.gatk.walkers.replication_validation;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 7/21/11
 * Time: 1:29 PM
 *
 * This is a site based implementation of a pool.
 * A pool object will contain all the information pertinent to a pool in one given site
 */

public class Pool {
    private String name;
    private ReadBackedPileup pileup;
    private int maxAlleleCount;
    private AlleleCountModel ACModel;
    private int matches;
    private int mismatches;


    @Deprecated
    public Pool(String name, ReadBackedPileup pileup, ErrorModel errorModel, byte referenceSequenceBase, int maxAlleleCount) {
        this(new PoolParameters(name, pileup, errorModel, referenceSequenceBase, maxAlleleCount));
    }
    public Pool(PoolParameters p) {
        name = p.name;
        pileup = p.pileup;
        maxAlleleCount = p.maxAlleleCount;

        byte [] data = pileup.getBases();
        int coverage = data.length;
        matches = MathUtils.countOccurrences(p.referenceSequenceBase, data);
        mismatches = coverage - matches;

        ACModel = new AlleleCountModel(new AlleleCountModelParameters(maxAlleleCount, p.errorModel, matches, mismatches));
    }

    /**
     * @return the name of the pool (sample name of the pool in the bam file)
     */
    public String getName() {
        return name;
    }

    /**
     * @param name The name of the pool (sample name of the pool in the bam file)
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return The pileup for the pool in its location
     */
    public ReadBackedPileup getPileup() {
        return pileup;
    }

    /**
     * @param pileup The pileup for the pool in its location
     */
    public void setPileup(ReadBackedPileup pileup) {
        this.pileup = pileup;
    }

    /**
     * @return the maximum allele count of this pool
     */
    public int getMaxAlleleCount() {
        return maxAlleleCount;
    }

    /**
     * @param maxAlleleCount the maximum allele count of this pool
     */
    public void setMaxAlleleCount(int maxAlleleCount) {
        this.maxAlleleCount = maxAlleleCount;
    }

    /**
     * @return the number of bases in the pool pileup
     */
    public int size() {
        return pileup.size();
    }

    public AlleleCountModel getACModel() {
        return ACModel;
    }

}
