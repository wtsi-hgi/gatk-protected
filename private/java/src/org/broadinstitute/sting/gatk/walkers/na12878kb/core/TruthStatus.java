package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

/**
 * Super-high level variant call status: true positive, false positive, or unknown
 *
 * User: depristo
 * Date: 11/4/12
 * Time: 6:47 AM
 */
public enum TruthStatus {
    /**
     * The site is a false positive
     */
    FALSE_POSITIVE("False positive"),

    /**
     * The site is a true positive
     */
    TRUE_POSITIVE("True positive"),

    /**
     * The truth status of the site is unknown
     */
    UNKNOWN("Unknown"),

    /**
     * We suspect the site is a false positive, but we aren't certain
     */
    SUSPECT("Suspect site"),

    /**
     * Information about the truth of the site is inconsistent
     */
    DISCORDANT("Discordant");

    final String statusText;

    private TruthStatus(String statusText) {
        this.statusText = statusText;
    }

    public boolean isFalsePositive() { return this == FALSE_POSITIVE; }
    public boolean isTruePositive() { return this == TRUE_POSITIVE; }
    public boolean isUnknown() { return this == UNKNOWN; }
    public boolean isSuspect() { return this == SUSPECT; }
    public boolean isDiscordance() { return this == DISCORDANT; }

    // ---------------------------------------------------------------------------
    //
    // The following code implements the makeConsensus function by a pairwise lookup table
    // that maps two TruthStatus types to its consensus result
    //
    // ---------------------------------------------------------------------------

    final static TruthStatus CONSENSUS_MATRIX[][] = new TruthStatus[TruthStatus.values().length][TruthStatus.values().length];
    private static void addToConsensusMatrix(final TruthStatus t1, final TruthStatus t2, final TruthStatus result) {
        CONSENSUS_MATRIX[t1.ordinal()][t2.ordinal()] = result;
        CONSENSUS_MATRIX[t2.ordinal()][t1.ordinal()] = result;
    }

    static {
        // by default everything is discordant
        for ( final TruthStatus t1 : TruthStatus.values() ) {
            for ( final TruthStatus t2 : TruthStatus.values() ) {
                addToConsensusMatrix(t1, t2, DISCORDANT);
            }
        }

        for ( final TruthStatus t : TruthStatus.values() ) {
            // combining two of the same type results in the same type
            addToConsensusMatrix(t, t, t);

            // combining t with unknown results in t
            addToConsensusMatrix(t, UNKNOWN, t);

            // combining t with SUSPECT results in SUSPECT
            addToConsensusMatrix(t, SUSPECT, SUSPECT);
        }

        // note this must occur after SUSPECT since SUSPECT + FALSE_POSITIVE = FALSE_POSITIVE
        for ( final TruthStatus t : TruthStatus.values() ) {
            // combining t with FALSE_POSITIVE results in FALSE_POSITIVE
            addToConsensusMatrix(t, FALSE_POSITIVE, FALSE_POSITIVE);
        }
    }

    /**
     * What is the consensus TruthStatus type results from combining this TruthStatus and otherStatus?
     *
     * For example, two TP sites are TP, one TP and one UNKNOWN is TP, TP + FP = DISCORDANT, etc
     *
     * @param otherStatus the TruthStatus of the other call
     * @return a consensus TruthStatus resulting from this and otherStatus
     */
    final TruthStatus makeConsensus(final TruthStatus otherStatus) {
        return CONSENSUS_MATRIX[this.ordinal()][otherStatus.ordinal()];
    }
}
