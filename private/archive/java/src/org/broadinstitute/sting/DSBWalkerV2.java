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

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.CircularArray;
import org.broadinstitute.sting.utils.PrimitivePair;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Dec 12, 2009
 * Time: 2:25:44 PM
 * To change this template use File | Settings | File Templates.
 */
@Requires({DataSource.READS, DataSource.REFERENCE})

public class DSBWalkerV2 extends LocusWalker<Integer,Integer> {
//    @Argument(fullName="coverage",shortName="C",doc="Regions with coverage above specified threshold will be reported",required=true)
//    int COV_CUTOFF = 0;
//    @Argument(fullName="minLength",shortName="ml",doc="Only regions longer than the specified value will be reported",required=false)
//    int MINLENGTH_CUTOFF = 0;
    @Argument(fullName="windowSize",shortName="W",doc="Size of the sliding window",required=true)
    int WINDOW_SIZE = 100;
    @Argument(fullName="enrichmentCutoff",shortName="E",doc="Report windows with enrichment (signal/control) above this cutoff",required=true)
    double ENRICHMENT_CUTOFF = 5.0;
    @Argument(fullName="minSignal",shortName="ms",doc="Do not report windows with signal lower than this value "+
                "(this cutoff is secondary to enrichmentCutoff and guards against windows where control signal is 0 or too low,"+
                "so that control*enrichmentCutoff is too low to be convincing)",required=true)
    int MIN_SIGNAL = 10;

    private CircularArray<PrimitivePair.Int> signalWindow = null;
    private CircularArray<PrimitivePair.Int> controlWindow = null;
    private CircularArray<PrimitivePair.Int> signalStrandsWindow = null;
    private CircularArray<PrimitivePair.Int> controlStrandsWindow = null;

    private PrimitivePair.Long totalSignalCoverage = new PrimitivePair.Long();
    private PrimitivePair.Long totalControlCoverage = new PrimitivePair.Long();
    private PrimitivePair.Long totalSignalFwdStrands = new PrimitivePair.Long();
    private PrimitivePair.Long totalControlFwdStrands = new PrimitivePair.Long();

    private Set<String> signalReadGroups; // we are going to remember which read groups are stimulated tagged and which are unstimulated untagged in order to be able
    private Set<String> controlReadGroups ; // to properly assign the reads coming from a merged stream

    private long windowStart = -1;
    private long windowStop = -1;
    private int curContig = -1;
    private String curContigName = "";

    // the following variables are for buffering and merging windows :
    private long regionStart = -1;
    private long lastWindowStart = -1;
    private PrimitivePair.Int maxSignalReads = new PrimitivePair.Int();
    private PrimitivePair.Int minSignalReads = new PrimitivePair.Int();
    private PrimitivePair.Int maxControlReads = new PrimitivePair.Int();
    private PrimitivePair.Int minControlReads = new PrimitivePair.Int();
    private double minEnrichmentUnique;
    private double maxEnrichmentUnique;
    private double minEnrichmentNonUnique;
    private double maxEnrichmentNonUnique;
    private double minEnrichmentTotal;
    private double maxEnrichmentTotal;
    private double minUniqueSignalStrandBalance = 0.0;
    private double maxUniqueSignalStrandBalance = 0.0;
    private double minNonUniqueSignalStrandBalance = 0.0;
    private double maxNonUniqueSignalStrandBalance = 0.0;
    private double minUniqueControlStrandBalance = 0.0;
    private double maxUniqueControlStrandBalance = 0.0;
    private double minNonUniqueControlStrandBalance = 0.0;
    private double maxNonUniqueControlStrandBalance = 0.0;

    @Override
    public void initialize() {
        int nSams = getToolkit().getArguments().samFiles.size();

        if ( nSams != 2 ) {
            out.println("ERROR: two input bam files (signal and backround control) must be specified");
            System.exit(1);
        }
        List<Set<String>> readGroupSets = getToolkit().getMergedReadGroupsByReaders();
        signalReadGroups = readGroupSets.get(0);
//        System.out.println(signalReadGroups.size()+" read groups in signal");
        controlReadGroups = readGroupSets.get(1);
//        System.out.println(controlReadGroups.size()+" read groups in control");
        signalWindow = new CircularArray<PrimitivePair.Int>(WINDOW_SIZE);
        controlWindow = new CircularArray<PrimitivePair.Int>(WINDOW_SIZE);
        signalStrandsWindow = new CircularArray<PrimitivePair.Int>(WINDOW_SIZE);
        controlStrandsWindow = new CircularArray<PrimitivePair.Int>(WINDOW_SIZE);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

         ReadBackedPileup pileup = context.getPileup();
         List<SAMRecord> reads = pileup.getReads();

         // compute coverages at the current site:
         PrimitivePair.Int signalCov = new PrimitivePair.Int();
         PrimitivePair.Int controlCov = new PrimitivePair.Int();
         PrimitivePair.Int signalFwdStrands = new PrimitivePair.Int();
         PrimitivePair.Int controlFwdStrands = new PrimitivePair.Int();

         for ( SAMRecord r : reads ) {
            if ( signalReadGroups.contains( r.getReadGroup().getReadGroupId() ) ) {
                if ( r.getMappingQuality() == 0 ) {
                    signalCov.second++;
                    if ( ! r.getReadNegativeStrandFlag() ) signalFwdStrands.second++;
                }
                else {
                    signalCov.first++;
                    if ( ! r.getReadNegativeStrandFlag() ) signalFwdStrands.first++;
                }
            } else {
                if ( controlReadGroups.contains( r.getReadGroup().getReadGroupId() ) ) {
                    if ( r.getMappingQuality() == 0 ) {
                        controlCov.second++;
                        if ( ! r.getReadNegativeStrandFlag() ) controlFwdStrands.second++;
                    }
                    else {
                        controlCov.first++;
                        if ( ! r.getReadNegativeStrandFlag() ) controlFwdStrands.first++;
                    }
                } else {
                    throw new StingException("Read "+r+" belongs to unknown read group ("+r.getReadGroup()+")");
                }
            }
         }

         GenomeLoc loc = context.getLocation();

     //    if ( curContig != 0 ) System.out.println(loc+" "+signalCov.first+" "+signalCov.second+" "+controlCov.first+" "+controlCov.second);

         if ( loc.getContigIndex() != curContig || loc.getStart() >= windowStop+WINDOW_SIZE ) {
             // we jumped to the next contig, or we are on the same contig but the current position is
             // more than WINDOW_SIZE away from the current window's end (i.e. there's nothing to shift)
             checkCurrentWindow(true);

             if ( loc.getContigIndex() != curContig ) {
                 System.out.println("on contig "+loc.getContig());
             }
             curContig = loc.getContigIndex();
             curContigName = loc.getContig();
//             prevPos = loc.getStart();
             windowStart = loc.getStart();
             windowStop = windowStart + WINDOW_SIZE - 1;
             signalWindow.clear();
             controlWindow.clear();
             totalSignalCoverage.assignFrom( signalCov );
             totalControlCoverage.assignFrom( controlCov );
             totalSignalFwdStrands.assignFrom( signalFwdStrands );
             totalControlFwdStrands.assignFrom( controlFwdStrands );
             signalWindow.set(0,signalCov);
             controlWindow.set(0,controlCov);
             signalStrandsWindow.set(0,signalFwdStrands);
             controlStrandsWindow.set(0,controlFwdStrands);
             return 1;
         }

         // offset of the current position w.r.t. the start of the window:
         int offset = (int)(loc.getStart() - windowStart); 

         if ( offset >= WINDOW_SIZE ) {
             // if we are here, the current position is outside of the current window, but not
             // far enough so that we'd need to reinitialize the window from scratch (that was already checked above).
             //  Now we need to shift.

             // We are receiving covered positions in order, so we are guaranteed that everything prior to
             // the current position was already counted; if some elements of the windows are still nulls, it means
             // there was no coverage there

             int shift = offset - WINDOW_SIZE + 1;

             // scroll the window(s) base by base until the current position is inside the window. At each step
             // we will check if the window meets the requirements and should be printed out.
             for ( int i = 0 ; i < shift ; i++ ) {

                 // we are going to shift; check if the window as it is now is worth printing
                 checkCurrentWindow(false);

                 // discard coverage from the first element of the window (this element is about to be shifted out of scope)
                 if ( signalWindow.get(0) != null ) totalSignalCoverage.subtract(signalWindow.get(0));
                 if ( signalStrandsWindow.get(0) != null ) totalSignalFwdStrands.subtract(signalStrandsWindow.get(0));

                 if ( controlWindow.get(0) != null ) totalControlCoverage.subtract(controlWindow.get(0));
                 if ( controlStrandsWindow.get(0) != null ) totalControlFwdStrands.subtract(controlStrandsWindow.get(0));

                 // advnace window coordinates on the ref
                 windowStart++;
                 windowStop++;

                 // shift the data in the window(s):
                 signalWindow.shiftData(1);
                 controlWindow.shiftData(1);
                 signalStrandsWindow.shiftData(1);
                 controlStrandsWindow.shiftData(1);

                 offset--; // this is the new offset w.r.t. to the shifted window
             }

         }

         // at this point, either the current position was inside the current window, or it was outside,
         // but the window was already shifted
         totalSignalCoverage.add(signalCov);
         totalControlCoverage.add(controlCov);
         totalSignalFwdStrands.add(signalFwdStrands);
         totalControlFwdStrands.add(controlFwdStrands);
         signalWindow.set(offset,signalCov);
         controlWindow.set(offset,controlCov);
         signalStrandsWindow.set(offset,signalFwdStrands);
         controlStrandsWindow.set(offset,controlFwdStrands);
         return 1;
    }


    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum+value;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void onTraversalDone(Integer result) {
        printRegion();
        super.onTraversalDone(result);
    }

    /** Checks if the currently held window satisfies the conditions set up for significance, and invokes buffered printout if so.
     * If the parameter is set to true, printout of previously held region is forced, and the buffer is reinitialized with
     * the new window if it passes the cutoffs, or left empty.
     *
     */
    private void checkCurrentWindow(boolean force) {
        if ( force ) printRegion();
        if ( signalWindow.get(0) == null && controlWindow.get(0) == null ) return; // do not emit windows that start from empty cell; we will get them later
        if ( totalControlCoverage.first * ENRICHMENT_CUTOFF / 36.0 < MIN_SIGNAL ) { // control coverage zero or too low
            if ( totalSignalCoverage.first /28.0 > MIN_SIGNAL ) emitWindow(false); // require at least MIN_SIGNAL coverage for signal
            return;
        }

        // if we have decent coverage in control, just check for required enrichment in the signal
        if ( ((double)totalSignalCoverage.first/28.0) / (totalControlCoverage.first/36.0) > ENRICHMENT_CUTOFF ) emitWindow(false);
    }

    /** This is actually a delayed print command: it buffers the successive windows set for printout, merges the windows that
     * are close enough and prints only when a train of close-by windows has ended and next window received is far enough
     */
    private void emitWindow(boolean force) {

        if ( regionStart == -1 ) {
            resetBuffer();
            return;
        }

        if ( force || windowStart > lastWindowStart + WINDOW_SIZE ) {
            // new window is far enough from the region we were buffering: emit old region

            printRegion();
            resetBuffer();
            return;
        }

        // current window is too close (overlapping) with a previous one: we need to merge

        lastWindowStart = windowStart;
        maxSignalReads.first = Math.max(maxSignalReads.first, (int)Math.round(totalSignalCoverage.first/28.0));
        maxSignalReads.second = Math.max(maxSignalReads.second,(int)Math.round(totalSignalCoverage.second/28.0));
        minSignalReads.first = Math.min(minSignalReads.first, (int)Math.round(totalSignalCoverage.first/28.0));
        minSignalReads.second = Math.min(minSignalReads.second,(int)Math.round(totalSignalCoverage.second/28.0));
        maxControlReads.first = Math.max(maxControlReads.first,(int)Math.round(totalControlCoverage.first/36.0));
        maxControlReads.second = Math.max(maxControlReads.second,(int)Math.round(totalControlCoverage.second/36.0));
        minControlReads.first = Math.min(minControlReads.first,(int)Math.round(totalControlCoverage.first/36.0));
        minControlReads.second = Math.min(minControlReads.second,(int)Math.round(totalControlCoverage.second/36.0));
        maxEnrichmentUnique = Math.max(maxEnrichmentUnique,((double)totalSignalCoverage.first/28.0)/(totalControlCoverage.first/36.0));
        minEnrichmentUnique = Math.min(minEnrichmentUnique, ((double)totalSignalCoverage.first/28.0)/(totalControlCoverage.first/36.0));
        maxEnrichmentNonUnique = Math.max(maxEnrichmentNonUnique,((double)totalSignalCoverage.second/28.0)/(totalControlCoverage.second/36.0));
        minEnrichmentNonUnique = Math.min( minEnrichmentNonUnique,  ((double)totalSignalCoverage.second/28.0)/(totalControlCoverage.second/36.0) );
        maxEnrichmentTotal = Math.max( maxEnrichmentTotal, ((double)(totalSignalCoverage.first+totalSignalCoverage.second)/28.0)/
                ((totalControlCoverage.first+ totalControlCoverage.second)/36.0) );
        minEnrichmentTotal = Math.min( minEnrichmentTotal, ((double)(totalSignalCoverage.first+totalSignalCoverage.second)/28.0)/
                ((totalControlCoverage.first+ totalControlCoverage.second)/36.0) );


        maxUniqueSignalStrandBalance = Math.max(maxUniqueSignalStrandBalance,((double)totalSignalFwdStrands.first)/totalSignalCoverage.first);
        minUniqueSignalStrandBalance = Math.min(minUniqueSignalStrandBalance,((double)totalSignalFwdStrands.first)/totalSignalCoverage.first);
        maxNonUniqueSignalStrandBalance = Math.max(maxNonUniqueSignalStrandBalance,((double)totalSignalFwdStrands.second)/totalSignalCoverage.second);
        minNonUniqueSignalStrandBalance = Math.min(minNonUniqueSignalStrandBalance,((double)totalSignalFwdStrands.second)/totalSignalCoverage.second);
        maxUniqueControlStrandBalance = Math.max(maxUniqueControlStrandBalance,((double)totalControlFwdStrands.first)/totalControlCoverage.first);
        minUniqueControlStrandBalance = Math.min(minUniqueControlStrandBalance,((double)totalControlFwdStrands.first)/totalControlCoverage.first);
        maxNonUniqueControlStrandBalance = Math.max(maxNonUniqueControlStrandBalance,((double)totalControlFwdStrands.second)/totalControlCoverage.second);
        minNonUniqueControlStrandBalance = Math.min(minNonUniqueControlStrandBalance,((double)totalControlFwdStrands.second)/totalControlCoverage.second);
        

    }

    private void resetBuffer() {
        regionStart = windowStart;
        lastWindowStart = windowStart;
        maxSignalReads.first = (int)Math.round(totalSignalCoverage.first/28.0);
        maxSignalReads.second = (int)Math.round(totalSignalCoverage.second/28.0);
        minSignalReads.assignFrom(maxSignalReads);
        maxControlReads.first = (int)Math.round(totalControlCoverage.first/36.0);
        maxControlReads.second = (int)Math.round(totalControlCoverage.second/36.0);
        minControlReads.assignFrom(maxControlReads);
        minEnrichmentUnique = maxEnrichmentUnique = ((double)totalSignalCoverage.first/28.0)/(totalControlCoverage.first/36.0);
        minEnrichmentNonUnique = maxEnrichmentNonUnique = ((double)totalSignalCoverage.second/28.0)/(totalControlCoverage.second/36.0);
        minEnrichmentTotal = maxEnrichmentTotal = ((double)(totalSignalCoverage.first+totalSignalCoverage.second)/28.0)/
                ((totalControlCoverage.first+ totalControlCoverage.second)/36.0);

        minUniqueSignalStrandBalance = maxUniqueSignalStrandBalance = ((double)totalSignalFwdStrands.first)/totalSignalCoverage.first;
        minNonUniqueSignalStrandBalance = maxNonUniqueSignalStrandBalance = ((double)totalSignalFwdStrands.second)/totalSignalCoverage.second;
        minUniqueControlStrandBalance = maxUniqueControlStrandBalance = ((double)totalControlFwdStrands.first)/totalControlCoverage.first;
        minNonUniqueControlStrandBalance = maxNonUniqueControlStrandBalance = ((double)totalControlFwdStrands.second)/totalControlCoverage.second;
    }

    private void printRegion() {
        if ( regionStart == -1 ) return;
        out.print(curContigName+":"+regionStart+"-"+windowStop+"\t"+(windowStop-regionStart+1) +"\t"+
                minSignalReads.first+"-"+maxSignalReads.first+"\t"+
                minSignalReads.second+"-"+maxSignalReads.second+"\t"+
                minControlReads.first+"-"+maxControlReads.first+"\t"+
                minControlReads.second+"-"+maxControlReads.second+"\t");
        out.printf("%.2f-%.2f\t",minEnrichmentUnique,maxEnrichmentUnique);
        out.printf("%.2f-%.2f\t",minEnrichmentNonUnique,maxEnrichmentNonUnique);
        out.printf("%.2f-%.2f\t",minEnrichmentTotal,maxEnrichmentTotal);
        out.printf("%.2f-%.2f\t",minUniqueSignalStrandBalance,maxUniqueSignalStrandBalance);
        out.printf("%.2f-%.2f\t",minNonUniqueSignalStrandBalance,maxNonUniqueSignalStrandBalance);
        out.printf("%.2f-%.2f\t",minUniqueControlStrandBalance,maxUniqueControlStrandBalance);
        out.printf("%.2f-%.2f",minNonUniqueControlStrandBalance,maxNonUniqueControlStrandBalance);

        if ( minUniqueSignalStrandBalance > 0.75 || minUniqueSignalStrandBalance < 0.25 ) out.print("\tS_U_STRAND_FILTER");
        out.println();

        regionStart = -1; // to indicate that there is nothing left to print, the buffer is empty
    }
}
