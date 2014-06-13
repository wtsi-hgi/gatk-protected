package org.broadinstitute.gatk.tools.walkers.siblingibd;

import htsjdk.variant.variantcontext.VariantContext;

/**
* Created by cwhelan on 5/27/14.
*
 * The IBD state at locus for a particular sibling pair
*/
public class VariantIBDState {
    final VariantContext vc;
    final SiblingPair siblingPair;
    final IBDState ibdState;
    final double[] statePosteriors;

    public VariantIBDState(final VariantContext vc, final SiblingPair siblingPair, final IBDState ibdState,  final double[] statePosteriors) {
        this.statePosteriors = statePosteriors;
        this.ibdState = ibdState;
        this.siblingPair = siblingPair;
        this.vc = vc;
    }

    @Override
    public String toString() {
        return "VariantIBDState: " + vc.getChr() + "\t" + vc.getStart() + "\t" + siblingPair + "\t" + ibdState;
    }
}
