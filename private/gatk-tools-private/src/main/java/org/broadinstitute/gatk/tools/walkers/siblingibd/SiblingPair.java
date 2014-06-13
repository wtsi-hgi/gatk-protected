package org.broadinstitute.gatk.tools.walkers.siblingibd;


import org.broadinstitute.gatk.engine.samples.Sample;

/**
* Created by cwhelan on 5/27/14.
 *
 * Organizes a pair of siblings for IBD analysis.
*/
class SiblingPair {
    final Sample sib1;
    final Sample sib2;

    public SiblingPair(final Sample sib1, final Sample sib2) {
        this.sib1 = sib1;
        this.sib2 = sib2;
    }

    public String getName() {
        return sib1.getID() + "-" + sib2.getID();
    }

    @Override
    public String toString() {
        return getName();
    }
}
