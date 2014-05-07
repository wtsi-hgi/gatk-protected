package org.broadinstitute.gatk.tools.walkers.siblingibd;

/**
* Created by cwhelan on 5/1/14.
 *
 * Class to represent the state of unambiguous minor allele sharing at a given locus between two siblings.
*/
public enum SharedMinorAlleleClass {
    HOMVAR_HOMVAR,
    HETVAR_HETVAR,
    HOMVAR_HETVAR,
    HETVAR_HOMREF,
    HOMREF_HOMREF,
    HOMVAR_HOMREF
}
