package org.broadinstitute.gatk.tools.walkers.siblingibd;


import htsjdk.variant.variantcontext.Genotype;
import org.broadinstitute.gatk.utils.exceptions.GATKException;

/**
 * Created by cwhelan on 5/8/14.
 *
 * Given the genotypes of the two parents and two siblings at a given locus,
 * this class represents the type of information about IBD state that can be
 * inferred. Based on the genotypes, each site can imply IBD0, IBD1, or IBD2,
 * or one of two of the IBD states, or leave the possibility of any of the three
 * IBD states.
 */
public enum IBDObservation {
    ZERO,
    ONE,
    TWO,
    ZERO_OR_ONE,
    ZERO_OR_TWO,
    ONE_OR_TWO,
    ZERO_OR_ONE_OR_TWO;

    private static boolean bothHom(final Genotype sib1Gt, final Genotype sib2Gt) {
        return sib1Gt.isHom() && sib2Gt.isHom();
    }

    private static boolean bothHet(final Genotype sib1Gt, final Genotype sib2Gt) {
        return sib1Gt.isHet() && sib2Gt.isHet();
    }

    /**
     * Given the four genotypes at a locus, return the type of IBD state observation
     *
     * @param p1 parental genotype 1
     * @param p2 parental genotype 2
     * @param s1 sibling genotype 1
     * @param s2 sibling genotype 2
     * @return the type of IBD observation indicated by the genotype configuration of the quartet
     */
    public static IBDObservation getIBDStateObservation(final Genotype p1, final Genotype p2, final Genotype s1, final Genotype s2) {
        if (bothHet(p1, p2)) {
            if (bothHom(s1, s2)) {
                if (s1.sameGenotype(s2)) {
                    return IBDObservation.TWO;
                } else {
                    return IBDObservation.ZERO;
                }
            } else if (bothHet(s1,s2)) {
                return IBDObservation.ZERO_OR_TWO;
            } else {
                // one sib het, the other hom
                return IBDObservation.ONE;
            }
        } else if (p1.isHet() && p2.isHom() || (p1.isHom() && p2.isHet())) {
            if (bothHet(s1, s2)) {
                return IBDObservation.ONE_OR_TWO;
            } else if (bothHom(s1,s2)) {
                return IBDObservation.ONE_OR_TWO;
            } else {
                // one sib het, the other hom
                return IBDObservation.ZERO_OR_ONE;
            }
        } else {
            if (! bothHom(p1,p2)) {
                throw new GATKException("Unexpected parental genotypes " + p1 + " and " + p2 + "; expect both to be het, both to be hom, or one het and one hom");
            }
            assert(bothHom(p1,p2));
            return IBDObservation.ZERO_OR_ONE_OR_TWO;
        }

    }

}
