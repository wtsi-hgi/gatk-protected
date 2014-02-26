/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

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
        int nComplex = 0;

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

            if ( mvc.isComplexEvent() )
                nComplex++;
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
        public int getnComplex() { return nComplex; }

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
            report.addRow(name, "n.Complex", summary.getnComplex());
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
