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

import org.broadinstitute.sting.utils.Interval;

/**  This class represents an indel as an interval with respect to the <i>original</i> reference and, in addition,
 * stores the indel type ( (I)nsertion or (D)eletion ) and can return meaningful event size (see below).
 * Depending on the indel type, the positions on the reference are:
 *    <ul>
 *    <li> Deletion ( e.g. deletion of ACAC from the ref: ACGTT[ACAC]TTTAG to ACGTT[]TTTAG) - start and stop of
 *         the interval are first and last deleted bases on the original reference (those in square brackets in the
 *         first sequence in the example)
 *    <li> Insertion ( e.g. insertion of GTGT into the ref: ACGTT{}TTTAG to ACGTT{GTGT}TTTAG) - start is the first
 *         position on the original reference <i>after</i> the insertion site (after the '}'), and stop is the last
 *         position <i>before</i> the insertion site (prior to '{').
 *    </ul>
 *
 * Given these definitions, the length of the interval, as returned by getLength() has the meaning of the length of
 * the event (affected bases) on the original reference: number of deleted bases for deletion and zero for insertion.
 * The length of the indel itself is returned by getIndelLength(), which is equal to getLength() for deletions and to
 * the actual number of inserted bases for insertions (while length on the reference, as returned by getLength() is zero).
 *
 * The overlaps are also meaningful with the above definitions: if an alignment to (or, in general, an interval on)
 * the original reference ends prior to <code>start</code>, or starts after <code>stop</code>, it does not overlap
 * with the indel event (neither spans over deleted region or contains any of the inserted bases).
 *
 */
public class Indel implements Interval {

    public static enum IndelType { I, D };
	
	private long mStart;
    private long mLength;
    private IndelType mType;

    /** Creates nBases-long indel at specified start position; the object will be unusable
     * until indel type is set.
     * @param start start position on the reference
     * @param nBases number of inserted or deleted bases
     */
    //public Indel(long start, long nBases) {
   //     mType=null;
   //     mStart=start;
   //     mLength=nBases;
   // }

    /** Creates nBases-long indel of the specified type (insertion or deletion), at specified start position.
     * @param start start position on the reference
     * @param nBases number of inserted or deleted bases
     * @param type Indel type: I or D.
     */
    public Indel(long start, long nBases, IndelType type) {
        mType=type;
        mStart=start;
        mLength=nBases;
    }

	/** Start coordinate on the reference; for deletions it is the position of the first deleted base,
	 * for insertions it is the first base after the insertion.
	 * This is the "left boundary" of the event on the original reference: every alignment that ends
     * befor this position on the reference does not overlap with the indel.
	 * @return indel's left boundary
	 */
	public long getStart() { return mStart; }

    /** Sets start position of the interval.
     *
     * @param s start coordinate
     */
	public void setStart(long s) { mStart = s; }

	/** Indel's stop coordinate on the reference; for deletions it is the position of the last deleted base,
	 * for insertions it is the last base before the insertion site (which makes it equal to getStart() - 1).
	 * This is the "right boundary" of the event: every alignment that starts after
     * this position on the reference
	 * does not overlap with the indel.
	 * @return indel's right boundary
	 */
	public long getStop() {
        if ( mType == IndelType.I ) return mStart - 1;
		else return mStart + mLength - 1;
	}

    /** This method is not supported in IndelInterval and will throw an exception. Use setIndelLength() instead.
     *
     * @param s stop coordinate
     */
    public void setStop(long s) {
        throw new UnsupportedOperationException("Method setStop(long) is not supported in IndelInterval");
    }

    /** Returns type of this indel ( I or D).
     *
     * @return I or D enum element
     */
	public IndelType getType() { return mType; }

    /** Sets the number of bases in this indel (i.e. the actual number of inserted or
     * deleted bases). Stop position will be always correctly computed based on the indel length and indel type.
     * @param nBases length of the indel (<i>not</i> the length of the event on the original reference!)
     */
	public void setIndelLength(long nBases) { mLength = nBases; }

    /** Returns actual number of inserted or deleted bases in the indel.
     *
     * @return number of bases (<i>not</i> the event length on the original reference).
     * @see #getLength()
     */
	public long getIndelLength() { return mLength; }

    /**
     * Returns true if this interval overlaps with i as judjed by getStart() and getStop() positions of the
     * two interval objects.
     *
     * @param i Another interval
     * @return true iff intervals overlap
     */

    public boolean overlapsP(Interval i) {
        return ! disjointP(i);  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Returns true if this interval does not overlap with i as judjed by getStart() and getStop() positions of the
     * two interval objects.
     *
     * @param i Another interval
     * @return true iff intervals do not overlap
     */
    public boolean disjointP(Interval i) {
        return i.getStop() < this.getStart() || i.getStart() > this.getStop();
    }

    /** Returns length of the region affected by the indel on the original reference. Note that an insertion
	 *  has length of 0.
     *  @return length of the event on the original, unmodified reference
	 */
	public long getLength() {
		if ( mType == IndelType.I ) return 0; 
		return mLength;
	}

    @Override
    public boolean equals(Object o) {
        if ( ! ( o instanceof Indel ) ) return false;
        Indel i = (Indel)o;
        return this.mType == i.mType && this.mStart == i.mStart && this.mLength == i.mLength ;
    }

    @Override
    public int hashCode() {
        return (int)( mStart << 6 + mStart + mLength );
    }
}
