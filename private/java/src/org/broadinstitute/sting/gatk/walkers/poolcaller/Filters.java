package org.broadinstitute.sting.gatk.walkers.poolcaller;

/**
 * Filters used in the Pool Caller framework
 *
 * @author Mauricio Carneiro
 * @since 7/28/11
 */
public enum Filters {
    LOW_QUAL,
    LOW_POWER,
    LOW_REFERENCE_SAMPLE_DEPTH,
    NO_BASES_IN_REFERENCE_SAMPLE,
    FILTERED_REFERENCE_SAMPLE_CALL
}
