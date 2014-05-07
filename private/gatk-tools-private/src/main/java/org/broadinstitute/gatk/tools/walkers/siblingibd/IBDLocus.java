package org.broadinstitute.gatk.tools.walkers.siblingibd;

import htsjdk.variant.variantcontext.VariantContext;

/**
* Created by cwhelan on 5/20/14.
 *
 * A simple data structure to hold the the predicted IBD state at a given locus.
*/
public class IBDLocus {
    VariantContext vc;
    IBDState state;
    double[] statePosteriors;

    public IBDLocus(final VariantContext vc, final IBDState state, final double[] statePosteriors) {
        this.vc = vc;
        this.state = state;
        this.statePosteriors = statePosteriors;
    }
}
