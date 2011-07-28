package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 7/21/11
 * Time: 3:41 PM
 *
 * A reference sample is a sample which true callset is known (from previous high coverage analysis)
 * and that is sequenced together with other pools to serve as a reference for estimating the site
 * specific error model.
 *
 * This is a site based implementation of the reference sample concept.
 */

public class ReferenceSample {
    private String name;
    private ReadBackedPileup pileup;
    private Collection<Byte> trueBases;

    /**
     * Creates a reference sample object
     */
    public ReferenceSample(ReferenceSampleParameters p) {
        name = p.name;
        pileup = p.pileup;
        trueBases = p.trueBases;
    }

    public String getName() {
        return name;
    }

    public ReadBackedPileup getPileup() {
        return pileup;
    }

    public Collection<Byte> getTrueBases() {
        return trueBases;
    }
}
