package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

/**
 * Is the site polymorphic in NA12878, monomorphic, or simply unknown?
 *
 * User: depristo
 * Date: 11/9/12
 * Time: 12:38 PM
 */
public enum PolymorphicStatus {
    POLYMORPHIC,
    MONOMORPHIC,
    DISCORDANT,
    UNKNOWN;

    final static PolymorphicStatus[][] consensusTransitions =  new PolymorphicStatus[PolymorphicStatus.values().length][PolymorphicStatus.values().length];
    private static void addTransition(final PolymorphicStatus s1, final PolymorphicStatus s2, final PolymorphicStatus result) {
        consensusTransitions[s1.ordinal()][s2.ordinal()] = result;
        consensusTransitions[s2.ordinal()][s1.ordinal()] = result;
    }

    static {
        // by default everything is discordant
        for ( final PolymorphicStatus s1 : PolymorphicStatus.values() )
            for ( final PolymorphicStatus s2 : PolymorphicStatus.values() )
                addTransition(s1, s2, DISCORDANT);

        for ( final PolymorphicStatus s : PolymorphicStatus.values() ) {
            // matching => same type
            addTransition(s, s, s);

            // unknown + X => X
            addTransition(s, UNKNOWN, s);
        }
    }

    /**
     * Given two PolymorphicStatus statuses (this and otherStatus) return a PolymorphicStatus reflecting their consensus opinion
     *
     * The consensus basically requires agreement between the statuses.  But if one is unknown the
     * status of the non-unknown is returned
     *
     * @param otherStatus
     * @return
     */
    public PolymorphicStatus makeConsensus(final PolymorphicStatus otherStatus) {
        return consensusTransitions[this.ordinal()][otherStatus.ordinal()];
    }

    public boolean isPolymorphic() { return this == POLYMORPHIC; }
    public boolean isMonomorphic() { return this == MONOMORPHIC; }
    public boolean isDiscordant() { return this == DISCORDANT; }
    public boolean isUnknown() { return this == UNKNOWN; }
}