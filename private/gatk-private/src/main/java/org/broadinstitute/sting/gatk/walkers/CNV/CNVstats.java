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

package org.broadinstitute.sting.gatk.walkers.CNV;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.Formatter;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;


/**
 * Walks along all variant ROD loci, and tabulates the statistics of the CNVs detected.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})
@By(DataSource.REFERENCE_ORDERED_DATA)

public class CNVstats extends RodWalker<CNVstatistics, CNVstatistics> {

    @Output(doc = "File to which copy number counts should be written", required = true)
    protected PrintStream out;

    @Argument(fullName = "alleleCountsCopyNumberFreqs", shortName = "AC_CNF", doc = "File to which discovered allele copy and copy number frequencies should be written", required = false)
    private PrintStream alleleCountsCopyNumberFreqs = null;

    @Argument(fullName = "minFracPassGt", shortName = "minFracPassGt", doc = "Minimum fraction of callable genotypes required to report any genotypes at all", required = false)
    private double minFracPassGt = 0.0;

    /**
     * All CNV variants found in these VCF files will be analyzed
     */
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    public static final String CNV_TAG = "<CNV>";
    public static final String CN_FIELD = "CN";

    public static final String SVLEN_FIELD = "SVLEN";
    public static final String AC_FIELD = "AC";

    public static final int DIPLOID = 2;

    public void initialize() {
    }

    public CNVstatistics reduceInit() {
        return new CNVstatistics();
    }

    /**
     * For each site, calculate the CNV stats.
     *
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return dummy Integer
     */
    public CNVstatistics map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        logger.debug("REF:" + ref.getLocus());
        CNVstatistics stats = new CNVstatistics();

        for (VariantContext vc : tracker.getValues(variantCollection.variants, context.getLocation())) {
            if (vc.isSymbolic() && vc.isBiallelic()) {
                Allele altAll = vc.getAlternateAllele(0);
                if (altAll.isSymbolic() && altAll.getDisplayString().equals(CNV_TAG)) {
                    logger.debug("Found CNV at locus...");
                    stats.cnvDeclaredLoci++;

                    CopyNumberCounts cnc = new CopyNumberCounts();

                    boolean hasDiploidGt = false;
                    boolean hasNonDiploidGt = false;
                    for (final Genotype gt : vc.getGenotypes()) {
                        int copyNum = gt.getAttributeAsInt(CN_FIELD, -1);
                        if (copyNum != -1 && ! gt.isFiltered()) {
                            cnc.incrementCopyNumber(copyNum);

                            if (copyNum == DIPLOID)
                                hasDiploidGt = true;
                            else
                                hasNonDiploidGt = true;
                        }
                    }

                    double calledFreq = ((double) cnc.calledCount()) / vc.getNSamples();
                    if (calledFreq < minFracPassGt) { // reset data as if it did not appear
                        cnc.resetCounts();
                    }
                    else {
                        if (hasDiploidGt && hasNonDiploidGt) {
                            stats.diploidAndNonDiploidLoci++;
                        }
                        else {
                            if (hasDiploidGt)
                                stats.diploidOnlyLoci++;
                            if (hasNonDiploidGt)
                                stats.nonDiploidOnlyLoci++;
                        }
                    }

                    int cnvEnd = vc.getEnd();
                    int cnvLength = vc.getAttributeAsInt(SVLEN_FIELD, -1);
                    if (cnvLength != -1)
                        cnvEnd = vc.getStart() + cnvLength - 1;
                    GenomeLoc vcLoc = getToolkit().getGenomeLocParser().createGenomeLoc(vc.getChr(), vc.getStart(), cnvEnd, true);
                    out.print(vcLoc);

                    for (Map.Entry<Integer, Integer> copyNumEntry : cnc.entrySet()) {
                        out.print("\t" + copyNumEntry.getKey() + ":" + copyNumEntry.getValue());
                    }
                    out.println();

                    if (alleleCountsCopyNumberFreqs != null) {
                        int ac = vc.getAttributeAsInt(AC_FIELD, -1);
                        CopyNumberCounts.DeletionDuplicationFreqs freqs = cnc.deletionDuplicationFreqs();
                        double cnvCount = freqs.deletionFreq + freqs.duplicationFreq;

                        alleleCountsCopyNumberFreqs.println(vcLoc + "\t" + ac + "\t" + freqs.deletionFreq + "\t" + freqs.duplicationFreq + "\t" + cnvCount);
                    }
                }
            }
        }

        return stats;
    }

    public CNVstatistics reduce(CNVstatistics result, CNVstatistics total) {
        if (result == null)
            return total;

        return total.addIn(result);
    }

    /**
     * @param result statistics of CNV sites
     */
    public void onTraversalDone(CNVstatistics result) {
        System.out.println();
        System.out.println("--------------------------------------");
        System.out.println("CNV summary:");
        System.out.println("--------------------------------------");

        System.out.println("cnvDeclaredLoci: " + result.cnvDeclaredLoci);

        System.out.println();
        System.out.println("noGenotypesLoci: " + result.noGenotypesLoci());

        System.out.println();
        System.out.println("nonDiploidOnlyLoci: " + result.nonDiploidOnlyLoci);

        System.out.println();
        System.out.println("lociWithDiploid: " + result.lociWithDiploid());
        System.out.println("diploidOnlyLoci: " + result.diploidOnlyLoci);
        System.out.println("diploidAndNonDiploidLoci: " + result.diploidAndNonDiploidLoci);
        String onlyDiploidRateStr = percentageString(result.diploidOnlyLoci, result.lociWithDiploid());
        System.out.println("onlyDiploidRate = " + onlyDiploidRateStr + "%");

        System.out.println();
        int noDiploidGenotypes = result.noGenotypesLoci() + result.nonDiploidOnlyLoci;
        System.out.println("loci with no diploid genotypes: " + noDiploidGenotypes);
        String noDiploidGtRateStr = percentageString(noDiploidGenotypes, result.cnvDeclaredLoci);
        System.out.println("noDiploidGtRate = " + noDiploidGtRateStr + "%");
    }

    private static String percentageString(int numerator, int denominator) {
        int NUM_DECIMAL_PLACES = 2;

        return new Formatter().format("%." + NUM_DECIMAL_PLACES + "f", MathUtils.percentage(numerator, denominator)).toString();
    }
}


class CNVstatistics {
    protected int cnvDeclaredLoci, diploidOnlyLoci, nonDiploidOnlyLoci, diploidAndNonDiploidLoci;

    public CNVstatistics() {
        this.cnvDeclaredLoci = 0;
        this.diploidOnlyLoci = 0;
        this.nonDiploidOnlyLoci = 0;
        this.diploidAndNonDiploidLoci = 0;
    }

    public CNVstatistics addIn(CNVstatistics other) {
        this.cnvDeclaredLoci += other.cnvDeclaredLoci;
        this.diploidOnlyLoci += other.diploidOnlyLoci;
        this.nonDiploidOnlyLoci += other.nonDiploidOnlyLoci;
        this.diploidAndNonDiploidLoci += other.diploidAndNonDiploidLoci;

        return this;
    }

    public int noGenotypesLoci() {
        return cnvDeclaredLoci - (diploidOnlyLoci + nonDiploidOnlyLoci + diploidAndNonDiploidLoci);
    }

    public int lociWithDiploid() {
        return diploidOnlyLoci + diploidAndNonDiploidLoci;
    }
}

class CopyNumberCounts {
    private Map<Integer, Integer> copyNumToCountsMap;
    private int calledCount;

    public CopyNumberCounts() {
        this.copyNumToCountsMap = new TreeMap<Integer, Integer>();
        this.resetCounts();
    }

    public void incrementCopyNumber(int copyNum) {
        Integer count = copyNumToCountsMap.get(copyNum);
        if (count == null)
            count = 0;

        copyNumToCountsMap.put(copyNum, count + 1);
        calledCount++;
    }

    public Set<Map.Entry<Integer, Integer>> entrySet() {
        return copyNumToCountsMap.entrySet();
    }

    public int calledCount() {
        return calledCount;
    }

    public void resetCounts() {
        copyNumToCountsMap.clear();
        calledCount = 0;
    }

    class DeletionDuplicationFreqs {
        public double deletionFreq;
        public double duplicationFreq;

        public DeletionDuplicationFreqs() {
            this.deletionFreq = 0;
            this.duplicationFreq = 0;
        }
    }

    public DeletionDuplicationFreqs deletionDuplicationFreqs() {
        int total = 0;
        DeletionDuplicationFreqs freqs = new DeletionDuplicationFreqs();

        for (Map.Entry<Integer, Integer> copyNumEntry : this.entrySet()) {
            int copyNum = copyNumEntry.getKey();
            int count = copyNumEntry.getValue();

            if (copyNum < CNVstats.DIPLOID) {
                freqs.deletionFreq += count;
            }
            else if (copyNum > CNVstats.DIPLOID) {
                freqs.duplicationFreq += count;
            }

            total += count;
        }

        freqs.deletionFreq /= total;
        freqs.duplicationFreq /= total;

        return freqs;
    }
}