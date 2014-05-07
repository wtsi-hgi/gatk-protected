package org.broadinstitute.gatk.tools.walkers.siblingibd;


import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;

/**
 * Created by cwhelan on 5/27/14.
 *
 * An interface to support varying models for predicting IBD state across the genome.
 */
public interface IBDStateModel {

    /**
     * Add a genotyped site to the model.
     *
     * @param vc
     * @param siblingPair
     * @param sib1Gt
     * @param sib2Gt
     * @return an empty list if we can't make any predictions right now; otherwise a list of sites at which we are ready to make IBD predictions
     */
    public List<VariantIBDState> addSite(final VariantContext vc, final SiblingPair siblingPair, final Genotype sib1Gt, final Genotype sib2Gt);

    /**
     * Call after all sites have been processed
     *
     * @return a list of any additional sites for which we can make IBD predictions
     */
    public List<VariantIBDState> finalizeModel();
}
