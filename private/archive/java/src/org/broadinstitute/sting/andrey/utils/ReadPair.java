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

package org.broadinstitute.sting.gatk.walkers.andrey.utils;

import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;

/**
 * Created by IntelliJ IDEA.
* User: asivache
* Date: Aug 6, 2010
* Time: 6:18:01 PM
* To change this template use File | Settings | File Templates.
*/
public class ReadPair {

    public enum PairType {
        UNKNOWN,
        BOTH_UNMAPPED,
        ONE_UNMAPPED,
        PROPER,
        LEFT,
        RIGHT,
        OUTER,
        INTER
    };


    private SAMRecord end1 = null;
    private SAMRecord end2 = null;
    private PairType pType = PairType.UNKNOWN;
    private int leftStart = -1;
    private int rightStart = -1;
    private SAMRecord leftRead = null;
    private SAMRecord rightRead = null;



    /** Creates an empty read pair object */
    public ReadPair() {}

    /** Creates a read pair objects initialized with the specified read */
    public ReadPair(SAMRecord read) {
        addRead(read);
    }

    /** Returns name of the paired read (it is assumed that both individual reads in the pair share same name).
     *
     * @return
     */
    public String getName() { return ( end1 != null ? end1.getReadName() : (end2 != null ? end2.getReadName() : null) ); }

    /** Returns true if both ends are recorded in this read pair object. Note that because SAM records carry 
     * mate information, a pair can be (partially) initialized from one end. This method verifies that this is not the case
     * and both records are actually present.
     * @return
     */
    public boolean hasBothEnds() { return end1 != null && end2 != null ; }

    /** Returns true if this pair object was initialized with at least one end. Since SAM records carry mate information,
     * it is sometimes sufficient to have only one read (fragment end) actually recorded in the pair object, at which
     * point some useful information can be retrieved for the pair already.
     * @return
     */
    public boolean hasAnyData() { return end1 != null || end2 != null ; }

    /** Returns true if both ends in the pair are mapped. The pair object must be at least partially initialized (i.e.
     * it has to hold a reference to at least one end of the pair), otherwise an exception will be thrown.
     * @return
     */
    public boolean bothEndsMapped() {
        if ( pType == PairType.UNKNOWN ) throw new StingException("ReadPair object was not initialized yet, method can not be applied");

        if ( pType == PairType.BOTH_UNMAPPED || pType == PairType.ONE_UNMAPPED ) return false;
        return true;
    }

    /** Returns true if both ends in the pair are mapped uniquely. This method requires both ends being already registered
     * in this pair object (i.e. hasBothEnds() is true), otherwise an exception will be thrown.
     * @return
     */
    public boolean bothEndsUniquelyMapped() {
        if ( ! hasBothEnds() ) throw new StingException("Can not determine if both ends are uniquely mapped until both ends are recorded");
        return bothEndsMapped() && end1.getMappingQuality() > 0 && end2.getMappingQuality() > 0;
    }

    /** Returns true if this pair is in proper orientation, i.e. ---> <--- on the same contig */
    public boolean isProper() { return pType == PairType.PROPER; }

    /* Returns true if this pair is in outer orientation, i.e. <--- ---> on the same chromosome */
    public boolean isOuter() { return pType == PairType.OUTER; }

    /** Returns left (coordinate-wise) read in the pair. Both ends need to be mapped, and they should map
     * onto the same contig, otherwise an exception will be thrown.
      * @return
     */
    public SAMRecord getLeftRead() {
        if ( ! bothEndsMapped() || pType == PairType.INTER )
            throw new StingException("Left read can be identified only when both reads are mapped onto the same contig, and the are not for "+getName());
        if ( leftRead == null )
            throw new StingException("Left read is not recorded. Maybe we have not seen it yet? Pair: "+getName());
        return leftRead;
    }

    /** Returns right (coordinate-wise) read in the pair. Both ends need to be mapped, and they should map
     * onto the same contig, otherwise an exception will be thrown.
      * @return
     */
    public SAMRecord getRightRead() {
        if ( ! bothEndsMapped() || pType == PairType.INTER )
            throw new StingException("Right read can be identified only when both reads are mapped onto the same contig, and the are not for "+getName());
        if ( rightRead == null )
            throw new StingException("Right read is not recorded. Maybe we have not seen it yet? Pair: "+getName());
        return rightRead;
    }

    public SAMRecord getEnd1() { return end1; }
    public SAMRecord getEnd2() { return end2; }

    public PairType getPairType() { return pType ; }

    public void addRead(SAMRecord r) {
        if ( ! r.getReadPairedFlag() ) throw new StingException("Read "+r.getReadName() +" is unpaired");
        if ( r.getFirstOfPairFlag() ) {
            if ( end1 != null ) throw new StingException("Read "+r.getReadName()+" is first of pair and the pair already has first read recorded");
            end1 = r;
            if ( end2 != null && ! end1.getReadName().equals(end2.getReadName()) )
                    throw new StingException("The pair already has read "+end2.getReadName() +"; the read being added does not match by name ("+r.getReadName()+")" );
        } else {
            if ( r.getSecondOfPairFlag() ) {
                if ( end2 != null ) throw new StingException("Read "+r.getReadName()+" is second of pair and the pair already has second read recorded");
                end2 = r;
                if ( end1 != null && ! end1.getReadName().equals(end2.getReadName()) )
                        throw new StingException("The pair already has read "+end1.getReadName() +"; the read being added does not match by name ("+r.getReadName()+")" );
            } else {
                throw new StingException("The read "+r.getReadName()+" is marked as paired, but the first/second of pair flag is not set");
            }
        }
        setPairInfo(r);
    }

    /** If pair type has not been set yet, then sets it to <code>t</code>. Otherwise (pair type already set),
     *  just checks if the pair type is <code>t</t>. If it is, the method returns quietly; if it is not (inconsistency detected),
     *  throws an exception. 
     *
     */
    private void setCheckPairType(PairType t) {
        if ( pType != PairType.UNKNOWN ) {
            if ( pType != t )
                throw new StingException("In pair "+getName()+" two ends provide conflicting alignment information");
        } else pType = t;
    }

    private void setCheckLeftStart(int pos) {
        if ( leftStart >= 0  ) {
            if ( leftStart != pos )
                throw new StingException("In pair "+getName()+" two ends provide conflicting alignment information");
        } else leftStart = pos;
    }

    private void setCheckRightStart(int pos) {
        if ( rightStart >= 0  ) {
            if ( rightStart != pos )
                throw new StingException("In pair "+getName()+" two ends provide conflicting alignment information");
        } else rightStart = pos;
    }

    private void setPairInfo(SAMRecord read) {

        setCheckPairType(getPairType(read));

        // there is nothing left to do unless both ends are mapped onto the same contig:
        if ( pType == PairType.INTER ) return;

        if ( pType == PairType.ONE_UNMAPPED ) {
            // set putative left or right read depending on the orientation of the only mapped mate
            if ( ! AlignmentUtils.isReadUnmapped(read ) ) {
                // we can set left/right read only if it is the current read that is mapped; if we have the
                // unmapped mate, skip and wait for the mapped read to come!
                if ( read.getReadNegativeStrandFlag() ) {
                    setCheckRightStart(read.getAlignmentStart());
                    if ( rightRead != null ) throw new StingException("Right read was already set for the pair");
                    rightRead = read;
                } else {
                    setCheckLeftStart(read.getAlignmentStart());
                    if ( leftRead != null ) throw new StingException("Left read was already set for the pair");
                    leftRead = read;                    
                }
            }
            return;
        }

        // we are here if both ends are mapped and they map onto the same contig
        if ( read.getAlignmentStart() < read.getMateAlignmentStart() ) { //left/right = read/mate

            setCheckLeftStart(read.getAlignmentStart());
            setCheckRightStart(read.getMateAlignmentStart());

            if ( leftRead != null ) throw new StingException("Left read was already set for the pair");
            leftRead = read;
        } else {
            // left/right = mate/read

            setCheckLeftStart(read.getMateAlignmentStart());
            setCheckRightStart(read.getAlignmentStart());

            if ( rightRead != null ) throw new StingException("Right read was already set for the pair");
            rightRead = read;
        }
    }

    /** Returns pair type that describes this read and its mate. The alignment information for both the read itself
     * and its mate is taken from the read's sam record passed as the argument, so the mate information is expected to be
     * correctly set!
     * @param read
     * @return
     */
    public static PairType getPairType(SAMRecord read) {

        if ( AlignmentUtils.isReadUnmapped(read) ) {
            if ( AlignmentUtils.isMateUnmapped(read) ) return PairType.BOTH_UNMAPPED;
            else return PairType.ONE_UNMAPPED;
        }

        return getWouldBePairType(read,read.getReferenceIndex(),read.getAlignmentStart(),read.getReadNegativeStrandFlag());
    }

    /** Returns pair type that would describe this read and its mate, if this read mapped onto refId:start in orientation
     * given by rc (forward is rc=false, reverse is rc=true). The read's alignment information (if any,
     * unmapped reads are allowed) present in the SAM record is completely ignored by this method,
     * only mate's information is used.
     * @param read
     * @param refId
     * @param start
     * @param rc
     * @return
     */
    public static PairType getWouldBePairType(SAMRecord read, int refId, int start, boolean rc) {


        if ( AlignmentUtils.isMateUnmapped(read) ) return PairType.ONE_UNMAPPED ;

        // both read and mate are mapped:

        if ( refId != read.getMateReferenceIndex() ) return PairType.INTER;

        // both read and its mate map onto the same chromosome

        if ( start < read.getMateAlignmentStart() ) { //left/right = read/mate

            if ( rc ) {
                if ( read.getMateNegativeStrandFlag() ) return PairType.LEFT;
                else return PairType.OUTER;
            } else {
                if ( read.getMateNegativeStrandFlag() ) return PairType.PROPER;
                else return PairType.RIGHT;
            }
        } else {
             // left/right = mate/read

             if ( rc ) {
                 if ( read.getMateNegativeStrandFlag() ) return PairType.LEFT;
                 else return PairType.PROPER;
             } else {
                 if ( read.getMateNegativeStrandFlag() ) return PairType.OUTER;
                 else return PairType.RIGHT;
             }
        }
    }

    public int getLeftStart() {
        if ( ! hasAnyData() ) throw new StingException("ReadPair object was not initialized yet, method can not be applied");
        return leftStart;
    }

    public int getRightStart() {
        if ( ! hasAnyData() ) throw new StingException("ReadPair object was not initialized yet, method can not be applied");
        return rightStart;
    }

    public int getFragmentSize() {
        if ( ! hasBothEnds() ) throw new StingException("Can not determine fragment size: pair object does not have both ends yet");
        if ( ! bothEndsMapped() ) throw new StingException("Can not determine fragment size: both ends must be mapped");
        if ( pType != PairType.PROPER ) throw new StingException("The pais is not in proper orientation, can not determine fragment size");

        return getFragmentSize(leftRead,rightRead);
    }

    /** Given a read (that must belong to this pair), returns the other end in the pair if it is already
     * recorded, or null otherwise.
     * @param read
     * @return
     */
    public SAMRecord getOtherEnd(SAMRecord read) {
        if ( read.getFirstOfPairFlag() ) return end2;
        else {
            if ( read.getSecondOfPairFlag() ) return end1;
        }
        return null;
    }

    public static int getFragmentSize(SAMRecord left, SAMRecord right) {

        if ( left == null || right == null ||
                AlignmentUtils.isReadUnmapped(left) || AlignmentUtils.isReadUnmapped(right) ) {
            throw new StingException("No read (null) or unmapped read provided: fragment size is not defined");
        }
        if ( !left.getReferenceIndex().equals(right.getReferenceIndex()) ) {
            throw new StingException("Left/right reads map onto different contigs: fragment size is not defined");
        }

        int fragment_length = left.getReadLength(); // fragment is at least as long as the left read, duh!
        int leftEnd = left.getAlignmentEnd();
        int rightStart = right.getAlignmentStart();

        if ( rightStart > leftEnd ) {
            // if reads are not overlapping, fragment length is lengths of both reads plus the distance (gap) between
            // the reads. Note that if the sequence between the reads happens to have insirtions or deletions,
            // our estimation of the actual distance between the reads (on the fragment) is incorrect, but we
            // can not do better given just those reads. This estimation is, in particular, incorrect
            // for left reads ending with 'I' and/or right reads starting with 'I'
            //
            //    left               right
            //  -------->...gap...<--------     fragment = left+gap+right

            return left.getReadLength() + right.getReadLength() + (rightStart - leftEnd-1);
        }

        // if we are here, the reads do overlap; fragment length is lengths of the two reads less the overlap.
        // in this case we can compute the actual overlap between the reads (on the fragment) taking into
        // account indels, if any
        //
        //      left    ****     right
        //     ------------>                   ****=overlap; fragment = left+right - overlap
        //              <--------------
        //
        // with deletion:
        //
        //      left    **   **     right
        //     -----------ddd->                   ****=overlap; fragment = left+right - overlap
        //              <-ddd-------------         note that overlap != leftEnd - rightStart+1
        //                                         instead, overlap = leftEnd-rightStart+1- length(D)
        // with insertion:
        //
        //     left      *******    right        ******* = overlap; fragment = left+right - overlap
        //    -------------iii->                   note that overlap != leftEnd - rightStart +1
        //               <-iii--------------       instead, overlap = leftEnd - rightStart +1 + length(I)
        //                                         (since 'i' bases are NOT on the ref)

        int posOnRef = rightStart;
//            int posOnRightRead = 0;

        int overlap = leftEnd - rightStart + 1 ;

        for(CigarElement ce : left.getCigar().getCigarElements() ) {
            switch(ce.getOperator()) {
                case S:
                case H:
//                      posOnRightRead+=ce.getLength();
                    break;
                case I:
                    overlap += ce.getLength();
                    break;
                case D:
                case N:
                    overlap -= ce.getLength();
                case M:
                case EQ:
                case X:
                    posOnRef += ce.getLength();
                    break;
                default:
            }
            if ( posOnRef > leftEnd ) break; // we need to examine only overlapping part of the reads
        }
        return left.getReadLength() + right.getReadLength() - overlap;
    }
}
