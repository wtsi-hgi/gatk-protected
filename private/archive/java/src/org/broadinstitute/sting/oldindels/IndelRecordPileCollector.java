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

package org.broadinstitute.sting.indels;

import java.util.*;


import org.broadinstitute.sting.indels.Indel.IndelType;
import org.broadinstitute.sting.utils.*;
import net.sf.samtools.*;

/** Ultimately, this class is a splitter for a stream of alignment records. It detects putative indels, or
 * trains of sufficiently close indels, and sends the alignments two-way: those that do not overlap with any
 * detected indels or trains of indels, and those that do. The latters are emitted in finished piles of
 * all alignments that overlap with genomic interval of interest. This collector should be bound to
 * and driven by an alignment traversal engine that sends alignment records one by one.
 *
 * NOTE 1: alignments must be sent to the collector strictly in the order
 * of non-decreasing reference start position.
 *
 * NOTE 2: a train of indels is defined as a sequence of (putative) indels such that each pair of adjacent indels
 * is overlapped by at least one alignment (that alignment does not have to have both indels in it, but only
 * to span over both positions). A "genomic region of interest" is defined as the smallest interval
 * containing all indels in the train, and <i>all</i> the alignments that overlap with that region will be collected
 * into one pile. For instance, if reads of different length are present in the dataset, it is possible that two
 * adjacent indels are overlapped by a single longer read (which stitches them into the train), but there are
 * shorter reads that fall completely into the in-between region (so that they technically do not overlap with any
 * indels in the train). According to the above definition, these shorter reads will be still emitted into the pile,
 * since they overlap with the "region of interest".
 *
 * NOTE 3: due to performance/memory issues, the collector may refuse to assemble a pile over pathologically long
 * train of indels. In this case, it will keep detecting the indel train in order to be able to understand what is
 * going on and to recover later, but the reads will be sent to the "non-overlapping" output channel.
 *
 * In order to abstract and decouple the operation of emitting records, the collector expects to be bound to an
 * implementation of RecordEmitter interface. It is the emitter's implementation that decides what to do with
 * alignments of the two types (not related to indels vs. piles of alignments overlapping with indels). While
 * this collector has some delay between receiving an alignment and being able to decide which way it should go,
 * no records are ever discarded.
 *
 * Implementation note:
 *
 * In order to achive its goal, the collector has a state ('wait' or 'active') and always
 * keeps a variable size "backlog" pile of alignments that were sent to it most recently. In 'wait' state collector
 * has not detected any putative indels just yet. The backlog pile contains only alignments of "undecided fate": those
 * that still might overlap with an indel should it be detected in the future. All alignments that end before the
 * current position on the genome have their fate determined (as not overlapping with any indels) and emitted.
 * When an indel is encountered, the collector flips into the 'active' state and from that moment on keeps all
 * the alignments in the pile and collects information on the indels (their positions on the reference and numbers
 * of observations).
 *
 * Since only alignments are sorted, but not indels (an indel in a later read may occur closer
 * to its start and thus before a previously seen indel), and also because it is relatively difficult (TO_DO?)to break a
 * pile in the middle immediately when it becomes clear that two adjacent indels might have been overlapped by a
 * single read, but no such read ever surfaced, the collector is conservative at this stage and keeps
 * accumulating the pile (and indel train) until it moves sufficiently far away from the last indel seen (full
 * maximum read length is recommended). Then it switches back into wait state and performs post-processing
 * of the indel train and the collected pile: only at this stage the preliminary pile is closely examined and if
 * there are pairs of adjacent indels not spanned by any read, the pile is broken into smaller piles
 * that conform to the contract outlined above. These piles are directed into the RecordEmitter, and the
 * reads that fall in between the piles, if any (i.e. those that do not overlap with final indel trains
 * determined at the post-processing stage) are relabeled as "not interesting" and redirected
 * to the appropriate output channel.
 *
 * @author asivache
 *
 */
public class IndelRecordPileCollector implements RecordReceiver {

    private final int WAIT_STATE = 0;
    private final int ACTIVE_STATE = 1;
    private int INDEL_COUNT_CUTOFF = 2;
    
    private boolean avoiding_region; // some regions are too funky (contain very long indel trains)-
    // we will determine their span and report them,
                                     // but we won't be counting any indels there or building piles

    private List<SAMRecord> mRecordPile; // here we keep the records before we decide how we want to emit them
    private TreeSet<CountedObject<Indel> > mAllIndels; ///< individual indels encountered, with observation counts
    private int mLastContig ;    ///< keeps the index of the contig last alignment was on
    private int mLastStartOnRef; ///< keeps the start position of the last alignment
    private int mState; ///< WAIT_STATE or ACTIVE_STATE
    private int mIndelSeparation; ///< Indels that are farther away from one another than this value
    ///< will be emitted separately; trains of indels with less then
    ///< mIndelSeparation bases between each adjacent pair will be emitted 
    ///< as one pile.
	
    // we will build histograms (distributions) of encountered indel lengths on the fly
    private List<Integer> mIndelLengthHistI;
    private List<Integer> mIndelLengthHistD;

    private RecordReceiver defaultReceiver; // we will send there records that do not overlap with regions of interest
    private RecordPileReceiver indelPileReceiver; // piles over indel regions will be sent there

    private boolean controlRun = false;

    public String memStatsString() {
		String s = "mRecordPile: ";
		return s+mRecordPile.size() + " mAllIndels: "+mAllIndels.size() + " mLastContig=" +mLastContig + " mLastStartOnref="+mLastStartOnRef;
                //+" Bndries="+mIndelRegionStart +":"+ mIndelRegionStop;
    }
	
    public void setIndelCountAcceptanceCutoff(int n) { INDEL_COUNT_CUTOFF = n; }
    
    public IndelRecordPileCollector(RecordReceiver rr, RecordPileReceiver rp) throws java.io.IOException {
        mRecordPile = new LinkedList<SAMRecord>();
        mAllIndels = new TreeSet<CountedObject<Indel> >(
              new CountedObjectComparatorAdapter<Indel>(new IntervalComparator()));
        mLastContig = -1;
        mLastStartOnRef = -1;
        mIndelSeparation = 51;
        mIndelLengthHistI = new ArrayList<Integer>();
        mIndelLengthHistD = new ArrayList<Integer>();
        for ( int i = 0 ; i < 5 ; i++ ) {
            mIndelLengthHistI.add(0);
            mIndelLengthHistD.add(0);
        }
        defaultReceiver = rr;
        indelPileReceiver = rp;
        setWaitState();
    }

    /** Fully reinitializes wait state: clears record pile and indel list, resets flags and states.
     * Does not emit records, just clears/resets the variables.
     */
    private void setWaitState() {
    	
        mRecordPile.clear();
        mAllIndels.clear();
//        mIndelRegionStart = 1000000000;
//        mIndelRegionStop = -1;
        avoiding_region = false;
        mState = WAIT_STATE; // got to do this if we were in avoid_region state
    }

    public void setControlRun(boolean c) { controlRun = c; }

    
    /** A utility method: emits into nonindelReceiver and purges from the currently held SAM record pile
     * all the consequtive records with alignment end positions less than or equal to the specified
     * position <code>pos</code>, until the first record is encountered that does not meet this condition. Note that
     * there might be more alignments that end at or before <code>pos</code> later on in the pile, but
     * they will <i> nit</i> be emitted/removed by this method.
     * @param pos all leading records with alignments ending before or at this position will be purged from the pile,
     * up to the first record that does not end at or before pos.
     */
    protected void purgeRecordsEndingAtOrBefore(final long pos) {
        Iterator<SAMRecord> i = mRecordPile.iterator();
        while ( i.hasNext() ) {
            SAMRecord r = i.next();
            if ( r.getAlignmentEnd() <= pos ) {
                defaultReceiver.receive(r);
                i.remove();
            } else break;
        }
    }
	
	/** A utility method: purges from the currently held SAM record pile all the records with alignment
	 * start positions greater than or equal to the specified position <code>pos</code>
	 * @param pos all records with alignments starting at or after this position will be purged from the pile
	 */
    protected void purgeRecordsStartingAtOrAfter(final int pos) {
        Iterator<SAMRecord> i = mRecordPile.iterator();
        while ( i.hasNext() ) {
            SAMRecord r = i.next();
            if ( r.getAlignmentStart() >= pos ) {
                defaultReceiver.receive(r);
                i.remove();
            } else break;
        }
    }
    
    /** This method MUST be called when no more reads are left in order to enforce the collector to emit the current pile of reads
     * it is still holding.
     */
    public void close() {
    	emit();
    }

    /** This is the main interface method of the collector: it receives alignments, inspects them, detects indels,
     *  updates and purges the read pile it keeps and emits alignments as needed.
     * Depending on the state, the following behaviors are possible
     *
     * <ul>
     *  <li> If the collector is in wait state (no indels seen recently): all
     *       alignments that end prior to the start of currently inspected alignment can not overlap
     *       with any future indels, including those that may be present in the current alignment; these records
     *       get purged from the pile and emitted immediately. Current alignment gets added to the pile.
     *       If current alignment has indels, collector switches into 'active' state.
     *  <li> in active state: if the current alignment starts sufficiently far away from the last indel seen,
     *       examine the currently held pile closely, split into a few separate piles/indel trains if needed, emit and
     *       completely purge the pile, add alignment to the pile, switch to wait state if alignment has no indels or
     *       stay in active state if it does. Otherwise (alignment too close to last indel),
     *       just add alignment to the pile, since it is yet impossible to tell whether new indels are coming soon and
     *       indel train will need to be extended; if alignment does have indels of its own, add them
     *       to the current indel train
     * </ul>
     *
     * This method checks that records arrive in reference-sorted order and throws RuntimeException if out-of-order
     * record arrives.
     *
     * @param r
     * @throws RuntimeException
     */
    public void receive(final SAMRecord r) throws RuntimeException {
		
        if ( r.getReadUnmappedFlag() ) {
        	defaultReceiver.receive(r); // do not throw reads away even if they are of no use for us, keep them in the output bam....
        	return; // read did not align, nothing to do
        }

        if ( controlRun ) {
            defaultReceiver.receive(r);
            return;
        }

        int currContig = r.getReferenceIndex();
        int currPos = r.getAlignmentStart();
		
        if ( currContig < mLastContig ) throw new RuntimeException("SAM file is not ordered by contigs");
        if ( currContig == mLastContig && currPos < mLastStartOnRef ) throw new RuntimeException("SAM file is not ordered by start positions");
	
        if ( currContig > mLastContig ) {
            // we jumped onto a new contig; emit everything we might have been building and purge the piles:
            emit();
        } else { // still on the same contig:

            switch (mState) {
                // everything ending up to currPos is guaranteed to have no overlaps with indels yet to come
            case WAIT_STATE: purgeRecordsEndingAtOrBefore(currPos); break;
                
                // next indel can start only after currPos (whether it is in the current read or in the
                // reads yet to come). If it is far enough from the last indel we have seen, we can emit
            case ACTIVE_STATE: if ( currPos - mAllIndels.last().getObject().getStop() > mIndelSeparation ) emit(); break;
            default: throw new RuntimeException("Unknown state");
            }
        }

        // does nothing if alignment has no indels, otherwise adds the indels to the list and (re)sets state to 'active'
        extractIndelsAndUpdateState(r.getCigar(),currPos);

        if ( mState == ACTIVE_STATE && ( ! avoiding_region ) && ( mAllIndels.size() > 20 || mRecordPile.size() > 1000 ) ) {
        	avoiding_region = true;
        }

        if ( ! avoiding_region ) mRecordPile.add(r); // add new record if this is not some crazy region
        else defaultReceiver.receive(r); // if we do not want to or can not deal with a region, pass reads through; 
                                                                   // the pile we have already collected before discovering it's a bad region will be sent through on the next call to emit() 

        mLastContig = currContig;
        mLastStartOnRef = currPos;
		
    }

    /** Emits all reads from the currently held pile, cleans the pile and fully reinitializes wait state
     * (clears indel list etc).
     *
     * If the current state is 'wait', simply sends all the records from the pile to nonindelReceiver before
     * the cleanup. If the state is 'active', then performs final inspection of the pile built over a train of indels,
     * splits the train (and the pile) into multiple trains/piles as needed (i.e. if there are pairs of adjacent
     * indels that are not overlapped by any read), and emits the final piles of records into indelReceiver.
      */
    private void emit() {

        if ( mState == WAIT_STATE || avoiding_region ) {
  //          System.out.println("Emitting uninteresting pile");
          if ( avoiding_region )  {
                long start = mAllIndels.first().getObject().getStart();
                long stop = mAllIndels.last().getObject().getStop();
                System.out.println("Genomic region "+mLastContig+":"+ start + "-"+ stop +
                                   " was ignored: ");
                System.out.println("   "+mAllIndels.size() +" unique indels with average distance of "+
                                   ((double)(stop - start))/((double)mAllIndels.size()-1) +
                                   " bases between indels");
                System.out.println("   "+mRecordPile.size() +" reads collected before aborting");
            }

            // no indels or avoiding indels in bad region: send all records to defaultReceiver and clear the pile
            for ( SAMRecord r : mRecordPile ) {
                defaultReceiver.receive(r);
            }
            setWaitState();
            return;
        }
        
        // last minute cleanup:
        // at this stage we have all the indels collected conservatively (in a sense
        // that they can be farther away than it is needed) - this means that there actually
        // can be more than one pile in what we have stored. Also, we can still have gapless reads
        // at the ends of the piles that do not really overlap with indel sites.
	
  //     System.out.println("Emitting pile with indels ("+mRecordPile.size()+" reads, "+mAllIndels.size()+" indels)");
        if ( mAllIndels.size() == 0 ) throw new RuntimeException("Attempt to emit pile with no indels");
        
        HistogramAsNeeded(mAllIndels);


        // indels are in a sorted map, and reads were added to the pile in the order they were received (also sorted).
        // we will traverse the two collections in parallel and detect exactly where we can break the indel train into
        // subtrains
        Iterator<CountedObject<Indel> > i_iter = mAllIndels.iterator();

        // will keep list of indels and list of records, respectively, in one final train
        List< CountedObject<Indel> > finalTrain = new ArrayList<CountedObject<Indel>>();
        List< SAMRecord > finalPile = new ArrayList<SAMRecord>();

        long curr_stop = -1; // the rightmost stop position among all the alignments seen so far

        CountedObject<Indel> indel = i_iter.next(); // we checked that list of indels contains at least one element!

        SAMRecord record ;

        while ( indel != null ) {

            // first, if we just started new indel train, then emit into defaultReceiver all alignments
            // that end prior to the first indel in the train:
            if ( finalTrain.size() == 0 ) purgeRecordsEndingAtOrBefore(indel.getObject().getStart() - 1);

            finalTrain.add(indel);

            Iterator<SAMRecord> r_iter = mRecordPile.iterator();

            if ( r_iter.hasNext() ) record = r_iter.next();
            else record = null;

            // record now contains first alignment that ends in or after the indel, or null if there are no more records

            // now collect all the alignments that overlap with the current indel (start before or inside) and
            // record the rightmost alignment stop position:
            while ( record != null && record.getAlignmentStart() <= indel.getObject().getStop() ) {
                finalPile.add(record);
                r_iter.remove(); // remove from the original pile the record we just moved to the current final pile
                curr_stop = Math.max(curr_stop, record.getAlignmentEnd());
                if ( r_iter.hasNext() ) record = r_iter.next();
                else record = null;
            }

            // record is now the first alignment that starts after the indel, or null if there are no more records

            // we are done with current indel, get next one if any:
            if ( i_iter.hasNext() ) {
                indel = i_iter.next();
            } else indel = null;
            if ( indel == null || curr_stop < indel.getObject().getStart() ) {
                // if there are no more indels or
                // all alignments that overlapped with the previous indel ended before the current indel started,
                // this means that the current train and pile of reads overlapping with it are fully built
                // and can be emitted

                if ( shouldAcceptForOutput(finalTrain ) ) {
                     System.out.print("SITE: " + mLastContig+":"+ finalTrain.get(0).getObject().getStart() + "-" +
                             finalTrain.get(finalTrain.size()-1).getObject().getStop() + " " +
                             finalTrain.size() + " indels; ");
                     System.out.print(finalPile.size() + " reads in the pile;")  ;
                     System.out.println(formatRange(finalTrain));
                     indelPileReceiver.receive(finalPile);
                } else {
                    for ( SAMRecord r : finalPile ) {
                        defaultReceiver.receive(r);
                    }
                }
                finalPile.clear();
                finalTrain.clear();
                curr_stop = -1;
            } // ELSE: otherwise we have reads that overlap with both previous and current indel, so we just continue
              // with building the indel train
        }

        // we may still have reads in the original pile that start after the last indel:
        for ( SAMRecord r : mRecordPile ) {
            defaultReceiver.receive(r);
        }

        setWaitState();
    }
		
	
    /** Looks for indels in the cigar and, if finds any, updates list of indels in the current train ans sets
     * the state to 'active'. If cigar contains no indels, this method does not do anything (it does <i>not</i>
     * set state back to 'wait' either!). If this method finds any indels in the cigar, it first tries to find them
     * in the list of previously seen indels. If the indel was already seen before, its counter is updated (indels
     * are stored in the list as counted objects), oherwise indel is added to the list with initial count of 1.
     *
     * @param c alignment cigar; if it contains no indels, nothing will be done
     * @param start position, at which the alignment represented by cigar <code>c</code> starts on the reference
     */
    private void extractIndelsAndUpdateState(final Cigar c, final int start) {
		//
		// firstpos,lastpos span of the indel will be interpreted as follows:
		// any alignment that ends strictly before firstpos or starts strictly after lastpos
		// on the *reference* (not inclusive!) does not overlap with an indel; in the case of 
		// insertion it will result in firstpos > lastpos!
		//         lastpos
		//         |   firstpos 
		//         |   |
		//         v   v
		// ---------III----- Ref  Insertion: bases I are not in the ref; any alignment that starts
		//                        after lastpos or ends before firstpos *on the reference*
		//                        is completely over the reference bases to the right or to
		//                        the left, respectively, of the insertion site
		//
		//      firstpos
		//      | lastpos
		//      | |
		//      v v
		//------------------ Ref   Deletion: any alignment that ends before firstpos or starts after lastpos
		// -----DDD--- alignment   on the reference does not overlap with the deletion 
		int runninglength = start; // position on the original reference; start = alignment start position

        if ( c.numCigarElements() == 1 ) return; // most of the reads have no indels, save a few cycles by returning early

        for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
			
            final CigarElement ce = c.getCigarElement(i);
            Indel indel = null;

            switch(ce.getOperator()) {
    		case I: indel = new Indel(runninglength, ce.getLength(), IndelType.I); break;
    		case D: indel = new Indel(runninglength, ce.getLength(), IndelType.D);
    				runninglength += ce.getLength();
    				break;
    		case M: runninglength += ce.getLength(); break; // advance along the gapless block in the alignment
    		default :
    			throw new IllegalArgumentException("Unexpected operator in cigar string");
            }

            if ( indel == null ) continue; // element was not an indel, go grab next element...

            mState = ACTIVE_STATE; // this is redundant and will be executed unnecessarily many times, but it's cheap...

            CountedObject<Indel> indelWithCount = new CountedObject<Indel>(indel);
            CountedObject<Indel> found = mAllIndels.floor(indelWithCount);
            
            if ( indelWithCount.equals( found ) ) found.increment(); // we did find our indel, advance the counter
            else mAllIndels.add(indelWithCount); // this is a new indel. Add it.
        } // end for loop over all alignment cigar elements
		
    } // end extractIndels() method



    /** Counts the size of the passed <indel> argument into the appropriate size histogram 
     * 
     * @param indel size of this indel will be counted in
     */
    private void addToSizeHistogram(Indel indel) {
        // count this indel's size into the appropriate bin of the appropriate histogram
        // (we count insertions and deletions separately), resizing the histogram array if needed:
        List<Integer> histogram;
        if ( indel.getType() == Indel.IndelType.D ) {
            histogram = mIndelLengthHistD; 
        } else if ( indel.getType() == Indel.IndelType.I ) {
            histogram = mIndelLengthHistI; 
        } else {
            throw new RuntimeException("Indel of unknown type");
        }
        if( indel.getIndelLength() > histogram.size() ) {
            for ( int j = histogram.size() ; j < indel.getIndelLength() ; j++ ) histogram.add(0);
            histogram.set((int)indel.getIndelLength()-1, 1); // we are seeing this length for the first time, so count == 1
        } else {
            int n = histogram.get((int)indel.getIndelLength()-1);
            histogram.set((int)indel.getIndelLength()-1, n+1);
        }		
    }
	
    /** Adds sizes of the indels from the list that pass some filters to the histograms
     * 
     * @param indels collection of indels with counts
     */
    private void HistogramAsNeeded(Collection<CountedObject<Indel>> indels) {
        for ( CountedObject<Indel> o : indels ) {
            if ( o.getCount() >= 2 ) addToSizeHistogram(o.getObject());
        }
    }
	
    /** Returns true if an attempt should be made to clean alignments around the
     * specified indel train; currently, indel run is acceptable
     * if it contains at least one indel onbserved more than once.
     * @param indels list of indels with counts to check for being acceptable
     * @return true if the indel run has to be printed
     */
    private boolean shouldAcceptForOutput(List<CountedObject<Indel>> indels) {
        for ( CountedObject<Indel> o :  indels ) {
            if ( o.getCount() >= INDEL_COUNT_CUTOFF ) return true;
        }
        return false;
    }

    private String formatRange(List<CountedObject<Indel>> indels) {
        StringBuffer b = new StringBuffer();
        StringBuffer all = new StringBuffer();
	
        long min = 1000000000;
        long max = 0;
        
        all.append("; passing indels:");
        for ( CountedObject<Indel> o :  indels ) {
            if ( o.getCount() < 2 ) continue;
            all.append(" ");
            all.append(o.getObject().getIndelLength());
            if ( o.getObject().getIndelLength() < min ) min = o.getObject().getIndelLength();
            if ( o.getObject().getIndelLength() > max ) max = o.getObject().getIndelLength();
        }
        if ( max == 0 ) return new String(); // no passinf indels, return empty string
        
        b.append(" passing min length: ");
        b.append(min);
        b.append("; passing max length: ");
        b.append(max);
        b.append(all);
        return b.toString();
    }

    public void printLengthHistograms() {
        if ( mIndelLengthHistD.size() < mIndelLengthHistI.size() ) {
            for ( int i = mIndelLengthHistD.size(); i < mIndelLengthHistI.size(); i++ ) mIndelLengthHistD.add(0);
        }
        if ( mIndelLengthHistI.size() < mIndelLengthHistD.size() ) {
            for ( int i = mIndelLengthHistI.size(); i < mIndelLengthHistD.size(); i++ ) mIndelLengthHistI.add(0);
        }
        System.out.println("length n_insertions n_deletions");
        for ( int i = 0 ; i < mIndelLengthHistD.size(); i++  ) {
            System.out.println((i+1)+" "+mIndelLengthHistI.get(i)+" "+mIndelLengthHistD.get(i));
        }
    }

    /** Returns true iff the SAM record (or, strictly speaking, its cigar) has at least one insertion or deletion
     *
     * @param r record to analyze
     * @return true if cigar contains at least one I or D element, false otherwise
     */
//    private boolean hasIndel(SAMRecord r) {
//        Cigar c = r.getCigar();
//        for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
//            CigarOperator co = c.getCigarElement(i).getOperator();
//            if ( co.equals(CigarOperator.I) || co.equals(CigarOperator.D) ) {
//                // we got an indel!
//                return true;
//            }
//        }
//        return false;
//    }
}
