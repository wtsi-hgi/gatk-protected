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

package org.broadinstitute.sting.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.DuplicateWalker;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.duplicates.DuplicateComp;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

class MismatchCounter {
    long nObs = 0;
    long nMismatches = 0;

    public void inc(long incNObs, long incNMismatches) {
        nObs += incNObs;
        nMismatches += incNMismatches;
    }

    public void inc(boolean mismatchP) {
        inc(1, mismatchP ? 1 : 0);
    }


    public double mismatchRate() {
        return (double)nMismatches / nObs;
    }

    public byte empiricalQualScore() {
        return QualityUtils.probToQual(1 - mismatchRate(), 0);
    }

    public String headerString() {
        return "mismatchRate\tempiricalQ\tnObs\tnMismatches";
    }

    public String toString() {
        return String.format("%.10f\t%d\t%d\t%6d", mismatchRate(), empiricalQualScore(), nObs, nMismatches);
    }
}

class QualityTracker {
    final private int MAX_QUAL_SCORE = 100;
    MismatchCounter[][] mismatchesByQ = new MismatchCounter[MAX_QUAL_SCORE][MAX_QUAL_SCORE];

    public QualityTracker() {
        for ( int i = 0; i < MAX_QUAL_SCORE; i++ ) {
            for ( int j = 0; j < MAX_QUAL_SCORE; j++ ) {
                mismatchesByQ[i][j] = new MismatchCounter();
            }
        }
    }

    public void inc(int b1Qi, int b2Qi, boolean mismatchP, boolean orderDependent) {
        int b1Q = orderDependent ? b1Qi : Math.max(b1Qi, b2Qi);
        int b2Q = orderDependent ? b2Qi : Math.min(b1Qi, b2Qi);

        if ( b1Q > MAX_QUAL_SCORE ) throw new RuntimeException("Unexpectedly large base quality " + b1Q);
        if ( b2Q > MAX_QUAL_SCORE ) throw new RuntimeException("Unexpectedly large base quality " + b2Q);

        mismatchesByQ[b1Q][b2Q].inc(mismatchP);
    }

    public void inc(DuplicateComp dc, boolean orderDependent) {
        inc(dc.getQLarger(), dc.getQSmaller(), dc.isMismatchP(), orderDependent);
    }

    public int probMismatchQ1Q2(int q1, int q2) {
        double e1 = 1 - QualityUtils.qualToProb(q1);
        double e2 = 1 - QualityUtils.qualToProb(q2);
        double eMM = e1 * (1 - e2) + (1 - e1) * e2 - 1/3 * e1 * e2;
        return QualityUtils.probToQual(1 - eMM, 0.0);
    }
    
    public void printToStream(PrintStream out, boolean filterUnobserved) {
        out.printf("Q1\tQ2\tQmin\t%s%n", mismatchesByQ[0][0].headerString());
        for ( int i = 0; i < MAX_QUAL_SCORE; i++ ) {
            for ( int j = 0; j < MAX_QUAL_SCORE; j++ ) {
                MismatchCounter mc = mismatchesByQ[i][j];
                //System.out.printf("MC = %s%n", mc);
                if ( filterUnobserved && mc.nObs == 0 )
                    continue;
                out.printf("%d\t%d\t%d\t%s\t%n", i, j, probMismatchQ1Q2(i,j), mc.toString());
            }
        }
    }
}

public class DuplicateQualsWalker extends DuplicateWalker<List<DuplicateComp>, QualityTracker> {
    @Argument(fullName="filterUnobservedQuals", required=false, doc="Show only quality bins with at least one observation in the data")
    public boolean FILTER_UNOBSERVED_QUALS = false;

    @Argument(fullName="maxPairwiseCompsPerDupSet", required=false, doc="Maximumize number of pairwise comparisons to perform among duplicate read sets")
    public int MAX_PAIRSIZE_COMPS_PER_DUPLICATE_SET = 100;

    @Argument(fullName="combinedQuals", required=false, doc="Combine and assess pairwise base qualities")
    public boolean COMBINE_QUALS = false;

    @Argument(fullName="combineAllDups", required=false, doc="Combine and assess pairwise base qualities")
    public boolean COMBINE_ALL_DUPS = false;

    @Argument(fullName="orderDependent", required=false, doc="")
    public boolean orderDependent = false;

    @Argument(fullName="compareToUniqueReads", required=false, doc="If true, then we will compare only to unique (i.e., non-duplicated molecules) at the same duplicate site")
    public boolean compareToUniqueReads = false;

    @Argument(fullName="comparePairToSingleton", required=false, doc="If true, then we will compare a combined dup to a random other read in the duplicate set, not a combined pair itself")
    public boolean comparePairToSingleton = false;

    final boolean DEBUG = false;
    final private boolean ACTUALLY_DO_WORK = true;

    public void onTraversalDone(QualityTracker result) {
        result.printToStream(out, FILTER_UNOBSERVED_QUALS);
    }

    public QualityTracker reduceInit() {
        return new QualityTracker();
    }

    public QualityTracker reduce(List<DuplicateComp> dupComps, QualityTracker tracker) {
        for ( DuplicateComp dc : dupComps ) {
            tracker.inc(dc, orderDependent);
        }
        
        return tracker;
    }

    // Print out data for regression
    public List<DuplicateComp> map(GenomeLoc loc, AlignmentContext context, Set<List<SAMRecord>> readSets ) {
        //logger.info(String.format("%s has %d duplicates and %d non-duplicates", loc, duplicateReads.size(), uniqueReads.size()));
        List<DuplicateComp> pairwiseComps = new ArrayList<DuplicateComp>();

        // todo -- fixme -- the logic here is all wrong given new interface        
//        if ( ! ACTUALLY_DO_WORK )
//            return pairwiseComps;
//
//        if ( COMBINE_QUALS ) {
//            Pair<SAMRecord, SAMRecord> combinedReads = DupUtils.combinedReadPair( duplicateReads );
//            if ( combinedReads != null ) {
//                SAMRecord combined1 = combinedReads.first;
//                SAMRecord combined2 = combinedReads.second;
//
//                if ( comparePairToSingleton )
//                    pairwiseComps = addPairwiseMatches( pairwiseComps, combined1, duplicateReads.get(2), uniqueReads );
//                else
//                    pairwiseComps = addPairwiseMatches( pairwiseComps, combined1, combined2, uniqueReads );
//            }
//        } else {
//            int nComparisons = 0;
//            for ( SAMRecord read1 : duplicateReads ) {
//                for ( SAMRecord read2 : duplicateReads ) {
//                    if ( read1.hashCode() < read2.hashCode() && DupUtils.usableDuplicate(read1, read2) ) {
//                        // the hashcode insures we don't do A vs. B and B vs. A
//                        //System.out.printf("Comparing %s against %s%n", read1, read2);
//                        nComparisons++;
//                        pairwiseComps = addPairwiseMatches( pairwiseComps, read1, read2, uniqueReads );
//                        if ( nComparisons > MAX_PAIRSIZE_COMPS_PER_DUPLICATE_SET )
//                            break;
//                    }
//                }
//            }
//        }

        return pairwiseComps;
    }
    
    private List<DuplicateComp> addPairwiseMatches(List<DuplicateComp> comps,
                                                   SAMRecord read1, SAMRecord read2,
                                                   List<SAMRecord> uniqueReads ) {
        if ( compareToUniqueReads ) {
            // we want to compare to a read in the unique read set
            if ( uniqueReads.size() > 0 ) {                 // there's actually something to compare to
                SAMRecord uniqueRead = uniqueReads.get(0);  // might as well get the first one
                return pairwiseMatches(comps, read1, uniqueRead);
            } else {
                return comps;
            }
        } else {
            // default, just do read1 vs. read2
            return pairwiseMatches(comps, read1, read2);
        }
    }

    /**
     * Calculates the pairwise mismatches between reads read1 and read2 and adds the result to the comps list.
     * Doesn't contain any logic deciding what to compare, just does read1 and read2
     * 
     * @param comps
     * @param read1
     * @param read2
     * @return
     */
    private List<DuplicateComp> pairwiseMatches(List<DuplicateComp> comps, SAMRecord read1, SAMRecord read2 ) {
        byte[] read1Bases = read1.getReadBases();
        byte[] read1Quals = read1.getBaseQualities();
        byte[] read2Bases = read2.getReadBases();
        byte[] read2Quals = read2.getBaseQualities();

        for ( int i = 0; i < read1Bases.length; i++) {
            byte qual1 = read1Quals[i];
            byte qual2 = read2Quals[i];
            boolean mismatchP = ! BaseUtils.basesAreEqual(read1Bases[i], read2Bases[i]);
            DuplicateComp dc = new DuplicateComp(qual1, qual2, mismatchP);
            comps.add(dc);
        }

        return comps;
    }
}