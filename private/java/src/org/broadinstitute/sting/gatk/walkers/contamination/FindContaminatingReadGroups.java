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

package org.broadinstitute.sting.gatk.walkers.contamination;

import cern.jet.stat.Probability;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.NamedTable;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * FindContaminatingReadGroupsWalker lists read groups in a single-sample BAM file that appear
 * to be contaminants (meaning a read group that's not actually associated with the sample) by searching
 * for evidence of systematic underperformance at likely homozygous-variant sites.
 *
 * @author Kiran Garimella
 */
public class FindContaminatingReadGroups extends LocusWalker<Integer, Integer> {
    @Output
    private PrintStream out;

    @Argument(fullName="balance", shortName="bal", doc="The expected alternate allele balance for homozygous-variant sites", required=false)
    private Double BALANCE = 0.95;

    @Argument(fullName="limit", shortName="lim", doc="The pValue limit for which a read group will be deemed to be a contaminant", required=false)
    private Double LIMIT = 1e-9;

    @Argument(fullName="scaleForSample", shortName="scale", doc="the scale by which the pvalue limit should reduce for testing samples directly. "+
              "E.g. if a sample has three 1e-3 read groups, pvalue is 1e-9 -- significant; so the scale should reduce by some multiplicative factor"+
              "For each read group associated with the sample. Defaults to 1e-4 [1e-9 for 1 RG, 1e-13 for 2 RG, 1e-17 for 3, etc]", required=false)
    private Double SCALE = 1e-4;

    private UnifiedGenotyperEngine ug;
    private NamedTable altTable;
    private final double EPSILON = 1e-20;

    public void initialize() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.STANDARD_CONFIDENCE_FOR_CALLING = uac.STANDARD_CONFIDENCE_FOR_EMITTING = 50.0;
        ug = new UnifiedGenotyperEngine(getToolkit(), uac);

        altTable = new NamedTable();
    }

    /**
     * Identify likely homozygous-variant sites that are called as
     * heterozygous, so that we can isolate our inspection to these sites.
     *
     * @param tracker  the meta-data tracker
     * @param ref      information regarding the reference
     * @param context  information regarding the reads
     * @return true if this site is a suspicious het, false if otherwise
     */
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        int altCount = 0;
        int totalCount = 0;

        ReadBackedPileup pileup = context.getBasePileup();
        int refIndex = BaseUtils.simpleBaseToBaseIndex(ref.getBase());

        for (byte base : pileup.getBases() ) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex((char) base);

            if (baseIndex != refIndex) {
                altCount++;
            }
            totalCount++;
        }

        double altBalance = ((double) altCount)/((double) totalCount);

        if (altBalance > 0.70) {
            VariantCallContext ugResult = ug.calculateLikelihoodsAndGenotypes(tracker, ref, context).get(0);

            if (ugResult != null && ugResult.getNSamples() > 0) {
                return ugResult.getGenotype(0).isHet();
            }
        }

        return false;
    }

    /**
     * For each read group represented in the pileup, determine the fraction of bases supporting the alternate allele
     *
     * @param tracker  the meta-data tracker
     * @param ref      information regarding the reference
     * @param context  information regarding the reads
     * @return 1
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        NamedTable alleleCounts = new NamedTable();

        int refIndex = BaseUtils.simpleBaseToBaseIndex(ref.getBase());
        String colName = String.format("%s.%d", context.getContig(), context.getPosition());

        for (int i = 0; i < context.size(); i++) {
            GATKSAMRecord read = context.getReads().get(i);
            int offset = context.getOffsets().get(i);

            SAMReadGroupRecord rg = read.getReadGroup();
            int alleleIndex = BaseUtils.simpleBaseToBaseIndex((char) read.getReadBases()[offset]);

            alleleCounts.increment(rg.getReadGroupId(), (alleleIndex == refIndex) ? "ref" : "alt");
        }

        for (String rg : alleleCounts.getRowNames()) {
            double altCount = alleleCounts.get(rg, "alt");
            double refCount = alleleCounts.get(rg, "ref");

            altTable.set(rg, colName, altCount / (altCount + refCount));
        }

        return 1;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return null;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    /**
     * Perform the t-test and list the read groups that are significant underperformers.
     *
     * @param result  the number of suspicious sites we're inspecting (this argument is ignored)
     */
    public void onTraversalDone(Integer result) {
        //out.println("readgroup\tpvalue\tstatus\tbalances");
        out.printf("%-10s\t%-13s\t%-10s\t%-10s%n", "readgroup", "pvalue", "status", "balances");

        HashMap<String,Double> pvalByReadGroup = new HashMap<String,Double>();
        for (String rg : altTable.getRowNames()) {
            String balances = "";

            // Compute mean
            double sum = 0.0, total = 0.0;

            for (String locus : altTable.getColumnNames()) {
                double value = altTable.get(rg, locus);

                sum += value;
                total += 1.0;

                balances += String.format("%2.2f,", value);
            }

            double mean = sum/total;

            // Compute stdev
            double squareSumOfMeanDifferences = 0.0;

            for (String locus : altTable.getColumnNames()) {
                double value = altTable.get(rg, locus);

                squareSumOfMeanDifferences += Math.pow(value - mean, 2.0);
            }

            double stdev = Math.sqrt(squareSumOfMeanDifferences/total);

            // Compute standard error of the mean (SEM)
            double sem = stdev/Math.sqrt(total);

            // Compute test statistic t
            double t = (mean - BALANCE) / sem;

            // Degrees of freedom
            double dof = total - 1.0;

            // Compute pValue
            double pValue = Probability.studentT(dof, t);
            pValue = pValue < EPSILON ? EPSILON : pValue;
            pvalByReadGroup.put(rg,pValue);

            //out.printf("%s\t%e\t%s\t[%s]\n", rg, pValue, (pValue < LIMIT ? "aberrant" : "nominal"), balances);
            out.printf("%-10s\t%-13s\t%-10s\t[%-10s]\n",
                       rg,
                       String.format("%e", pValue),
                       (pValue < LIMIT ? "aberrant" : "nominal"),
                       balances);

            logger.debug(rg);
        }

        out.printf("%n%n%s%n","SECTION ON BADLY CONTAMINATED SAMPLES");
        out.printf("%s\t%s\t%s\t%s%n","sample","p-value","status","info");

        HashMap<String,List<String>> samplesToReadGroups = new HashMap<String,List<String>>();
        for ( SAMReadGroupRecord rec : getToolkit().getSAMFileHeader().getReadGroups() ) {
            if ( samplesToReadGroups.containsKey(rec.getSample()) ) {
                samplesToReadGroups.get(rec.getSample()).add(rec.getReadGroupId());
            } else {
                ArrayList<String> newList = new ArrayList<String>();
                newList.add(rec.getReadGroupId());
                samplesToReadGroups.put(rec.getSample(),newList);
            }
        }

        for ( String sample : samplesToReadGroups.keySet() ) {
            double p_value = 1;
            double limit = LIMIT;
            boolean containsAberrantReads = false;
            for ( String rg : samplesToReadGroups.get(sample) ) {
                double rg_pval = ( pvalByReadGroup.get(rg) == null ? 1 : pvalByReadGroup.get(rg) );
                p_value = p_value*rg_pval;
                containsAberrantReads = containsAberrantReads || rg_pval < LIMIT;
                limit = limit*SCALE;
                logger.debug(rg);
            }

            out.printf("%s\t%-13s\t%s\t%s%n", sample, String.format("%e",p_value), ( p_value < limit ? "aberrant" : "nominal"), ( containsAberrantReads ? "contains_aberrant_RG" : "no_aberrant_RG"));
        }
    }
}
