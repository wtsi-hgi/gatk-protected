package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import com.google.java.contract.Ensures;
import org.broadinstitute.sting.gatk.report.GATKReport;

import java.util.Date;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.Map;

/**
 * Summarizes information about consensus callsets
 *
 * User: depristo
 * Date: 11/29/12
 * Time: 9:04 AM
 */
public class ConsensusSummarizer {
    private final static String CONSENSUS = "Consensus";

    /**
     * Collects information about a specific callset
     */
    final static class CallSetSummary {
        final String callSetName;

        final EnumMap<TruthStatus, Integer> countsByTruth = new EnumMap<TruthStatus, Integer>(TruthStatus.class);
        final EnumMap<PolymorphicStatus, Integer> countsByPoly = new EnumMap<PolymorphicStatus, Integer>(PolymorphicStatus.class);

        int nTPPoly = 0;
        int nNonSingletonTPPoly = 0;
        int nSingletonTPPoly = 0;
        int nSingletons = 0;
        int nSNPs = 0;
        int nIndels = 0;

        CallSetSummary(String callSetName) {
            if ( callSetName == null ) throw new IllegalArgumentException("callSetName cannot be null");

            this.callSetName = callSetName;
            for ( final PolymorphicStatus pStatus : PolymorphicStatus.values() )
                countsByPoly.put(pStatus, 0);
            for ( final TruthStatus tStatus : TruthStatus.values() )
                countsByTruth.put(tStatus, 0);
        }

        /**
         * Update the statistics about this callset from a consensus MongoVariantContext that includes it
         * @param mvc a non-null consensus call that includes callSetName in its supporting callsets
         */
        final void update(final MongoVariantContext mvc) {
            if ( mvc == null ) throw new IllegalArgumentException("mvc cannot be null");
            if ( ! mvc.getSupportingCallSets().contains(callSetName) && ! callSetName.equals(CONSENSUS) )
                throw new IllegalArgumentException("Trying to include MVC from the wrong callset " + mvc + " into summary for " + callSetName);

            countsByTruth.put(mvc.getType(), countsByTruth.get(mvc.getType()) + 1);
            countsByPoly.put(mvc.getPolymorphicStatus(), countsByPoly.get(mvc.getPolymorphicStatus()) + 1);


            if ( mvc.isPolymorphic() && mvc.getType() == TruthStatus.TRUE_POSITIVE ) {
                nTPPoly++;
                if ( mvc.isSingleCallset() ) {
                    nSingletonTPPoly++;
                } else
                    nNonSingletonTPPoly++;
            }

            if ( mvc.isSingleCallset() )
                nSingletons++;

            if ( mvc.getVariantContext().isSNP() )
                nSNPs++;
            else
                nIndels++;
        }

        @Ensures("result != null")
        public String getCallSetName() {
            return callSetName;
        }

        /**
         * Add information about counts and frequency for status to a GATK report
         * @param report the GATK report to add a row to
         * @param status the polymorphic status to report on
         */
        public void addRow(final GATKReport report, final PolymorphicStatus status) {
            report.addRow(callSetName, "polymorphic." + status + ".count", getN(status));
            report.addRow(callSetName, "polymorphic." + status + ".percent", getPercent(status));
        }

        /**
         * Add information about counts and frequency for status to a GATK report
         * @param report the GATK report to add a row to
         * @param status the truth status to report on
         */
        public void addRow(final GATKReport report, final TruthStatus status) {
            report.addRow(callSetName, "polymorphic." + status + ".count", getN(status));
            report.addRow(callSetName, "polymorphic." + status + ".percent", getPercent(status));
        }

        @Ensures("result >= 0")
        public int getN(final PolymorphicStatus status) { return countsByPoly.get(status); }

        @Ensures("result >= 0")
        public int getN(final TruthStatus status) { return countsByTruth.get(status); }

        @Ensures("result != null")
        public String getPercent(final PolymorphicStatus status) { return percent(countsByPoly.get(status)); }

        @Ensures("result != null")
        public String getPercent(final TruthStatus status) { return percent(countsByTruth.get(status)); }

        /**
         * What percent of all sites (nTotal) are included in this TPPoly sites?
         *
         * @param nTotalTPPoly the total number of TPPoly sites
         * @return the percent (as a formatted string) of this TP/poly sites of all polymorphic sites in the consensus
         */
        @Ensures("result != null")
        public String getPercentOfAllSites(final int nTotalTPPoly) {
            return percent(nTPPoly, nTotalTPPoly);
        }

        /**
         * How many total sites are present in this callset?
         * @return
         */
        @Ensures("result >= 0")
        public int getNSites() {
            int c = 0;
            for ( final PolymorphicStatus status : PolymorphicStatus.values() )
                c += countsByPoly.get(status);
            return c;
        }

        @Ensures("result != null")
        public String percent(final int n) {
            return percent(n, getNSites());
        }

        @Ensures("result != null")
        public String percent(final int n, final int d) {
            return d == 0 ? "NA" : String.format("%.2f", n / (d * 0.01));
        }

        @Ensures("result >= 0")
        public int getnSingletonTPPoly() { return nSingletonTPPoly; }
        @Ensures("result >= 0")
        public int getnSingletons() { return nSingletons; }

        @Ensures("result >= 0")
        public int getnSNPs() { return nSNPs; }
        @Ensures("result >= 0")
        public int getnIndels() { return nIndels; }

        @Ensures("result >= 0")
        public int getnTPPoly() {
            return nTPPoly;
        }

        @Ensures("result >= 0")
        public int getnNonSingletonTPPoly() {
            return nNonSingletonTPPoly;
        }
    }

    /** Total number of sites seen */
    int nSites = 0;

    /** total TP / poly sites are in the consensus? */
    int totalPolyTPSites = 0;

    /** The call set summaries for every callset observed in the consensus */
    final Map<String, CallSetSummary> callSetSummaries = new HashMap<String, CallSetSummary>();

    /**
     * Get a CallSetSummary object for name, creating one if it hasn't been see yet
     * @param name the call set's name
     * @return a CallSetSummary to be used with name
     */
    @Ensures("result != null")
    private CallSetSummary getSummary(final String name) {
        if ( name == null ) throw new IllegalArgumentException("Name cannot be null");

        if ( ! callSetSummaries.containsKey(name) )
            callSetSummaries.put(name, new CallSetSummary(name));
        return callSetSummaries.get(name);
    }

    /**
     * Update the summary statistics based on the consensus call consensus
     * @param consensus a consensus call as a MongoVariantContext
     */
    public void add(final MongoVariantContext consensus) {
        if ( consensus == null ) throw new IllegalArgumentException("Consensus cannot be null");

        nSites++;
        if ( consensus.getType() == TruthStatus.TRUE_POSITIVE && consensus.isPolymorphic() )
            totalPolyTPSites++;
        for ( final String callset : consensus.getSupportingCallSets() ) {
            final CallSetSummary summary = getSummary(callset);
            summary.update(consensus);
        }

        getSummary(CONSENSUS).update(consensus);
    }

    /**
     * Generate a summary GATKReport using the current data in this ConsensusSummarizer
     * @param detailed should the report include lots of details or just high-level information?
     * @return a GATKReport containing the information from this ConsensusSummarizer
     */
    @Ensures("result != null")
    public GATKReport summaryGATKReport(final boolean detailed) {
        final GATKReport report = GATKReport.newSimpleReportWithDescription("ConsensusSummary",
                "Date " + new Date(), "CallSet", "Variable", "Value");

        for ( final CallSetSummary summary : callSetSummaries.values() ) {
            final String name = summary.getCallSetName();
            report.addRow(name, "total.sites", summary.getNSites());
            report.addRow(name, "percent.of.all.tp.poly.sites", summary.getPercentOfAllSites(totalPolyTPSites));
            report.addRow(name, "n.SNPs", summary.getnSNPs());
            report.addRow(name, "n.Indels", summary.getnIndels());
            report.addRow(name, "poly.nonsingleton.count", summary.getnNonSingletonTPPoly());
            report.addRow(name, "poly.nonsingleton.percent", summary.percent(summary.getnNonSingletonTPPoly()));
            report.addRow(name, "poly.singletons.count", summary.getnSingletonTPPoly());
            report.addRow(name, "poly.singletons.percent", summary.percent(summary.getnSingletonTPPoly()));
            report.addRow(name, "singletons.count", summary.getnSingletons());
            report.addRow(name, "singletons.percent", summary.percent(summary.getnSingletons()));

            if ( detailed ) {
                for ( final PolymorphicStatus pStatus : PolymorphicStatus.values() )
                    summary.addRow(report, pStatus);

                for ( final TruthStatus tStatus : TruthStatus.values() )
                    summary.addRow(report, tStatus);
            }
        }

        return report;
    }

    /**
     * Get the total number of TP / polymorphic consensus sites
     * @return
     */
    @Ensures("result >= 0")
    public int getTotalPolyTPSites() {
        return totalPolyTPSites;
    }

    /**
     * Get the total number of consensus sites
     * @return
     */
    @Ensures("result >= 0")
    public int getnSites() {
        return nSites;
    }
}
