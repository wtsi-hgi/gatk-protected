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
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pair;


public class MultipleAlignment implements Iterable<Integer>  {
	private static final int IMPOSSIBLE = 1000000000;
	private Map<Integer,Integer> index;    // maps external id of the sequence onto its index in the pile
	private List<String> seqs;             // sequences, in order they were added
	private List<Integer> ext_ids;         // external ids of the sequences, in order they were added to the pile
	private List<Integer>  alignment_offsets; // offset of seqs[i] w/respect to seqs[0] (i.e. in order the seqs were added)
	private int best_mm; // mismatch count
	private int next_mm; // next-best mismatch count
	private ConsensusSequence consensus;

	public MultipleAlignment() {
		index = new HashMap<Integer,Integer>();
		seqs = new ArrayList<String>();
		alignment_offsets = new ArrayList<Integer>();
		ext_ids = new ArrayList<Integer>();
        consensus = new ConsensusSequence(); // we use reference position 0, e.g. we hook onto the first read in the pile
	}

	public void clear() {
		seqs.clear();
		index.clear();
		alignment_offsets.clear();
		ext_ids.clear();
	}

	/** Adds  single sequence with id set to i. Pile must be empty, or IllegalStateException will be thrown
	 * 
	 * @param seq sequence to add
	 * @param i id of the sequence (can be use later to query the pile)
     * @see #add(String,int,int)
	 */
	public void add( String seq, int i ) throws IllegalStateException {
		if ( size() != 0 ) throw new IllegalStateException("Single sequence can be added to an empty pile only");
        add(seq,i,0);
	}
	
    /** Adds  single sequence with id set to i and places it at the specified offset wrt the first sequence
     * in this pile (i.e. wrt reference position 0).
     *
     * @param seq sequence to add
     * @param i id of the sequence (can be use later to query the pile)
     * @see #add(String,int)
     */
    public void add( String seq, int i, int offset ) throws IllegalStateException {
        index.put(i,index.size());
        ext_ids.add(i);
        seqs.add(seq);
        alignment_offsets.add(offset);
        consensus.addSequence(seq,offset);
    }

	public void add( PairwiseAlignment a) {
		if ( a.id1() == -1 || a.id2() == -1 ) throw new IllegalArgumentException("Attempt to add pairwise alignemnt with sequence ids not properly set");
		add(a,a.id1(),a.id2());
	}
	
	/** Adds pair of aligned sequences to the pile, with the external ids of the first and second sequences being i and j,
	 * respectively. Pairwise alignment can be always added to an empty pile. If the pile is non-empty, exactly
     * one of the sequences held by the pair-wise alignment should be already in the pile; this sequence (and the
     * pairwise alignment itself) will be used to stitch the other sequence to the pile. If either both or
	 * none of the specified ids are already in the pile, an IllegalStateException will be thrown.
	 * @param a
	 * @param i
	 * @param j
	 */
	public void add( PairwiseAlignment a, int i, int j ) throws IllegalStateException {
		if ( seqs.size() == 0 ) {
            add(a.getSequence1(),i,0);
            add(a.getSequence2(),j,a.getBestOffset2wrt1());
			return;
		}
		
		Integer first = index.get(i);
		Integer second = index.get(j);
		
		if ( first != null && second != null ) {
			throw new IllegalStateException("Attempt to add pairwise alignment for two sequences that are already in the pile");
		}

		if ( first == null && second == null ) {
			throw new IllegalStateException("Attempt to add pairwise alignment for two sequences none of which is already in the pile");
		}
		
		if ( second == null ) add(a.getSequence2(),j, a.getBestOffset2wrt1() + alignment_offsets.get( first ) );
		else add(a.getSequence1(),i, -a.getBestOffset2wrt1() + alignment_offsets.get( second ) );
	}

    /** Adds another pile of aligned sequences to this pile, stitching them together using specified pairwise alignment
     * p of the sequences with external ids i and j. One of the indices i, j must be in this pile, and the other in
     * the pile being added, otherwise an IllegalArgumentException is thrown. Sequence id's i and j MUST be the ids
     * of the first and second sequences in the pairwise alignment, in that order. Specified ids override
     * ids, if any, set for the sequences in the pairwise alignment; it is not checked whether the specified and
     * stored ids match. The piles can not overlap.
     */
    public void add(MultipleAlignment a, PairwiseAlignment p, int i, int j) {
        int off2; // offset of the first sequence in pile 'a' wrt the first sequence in this pile
        if ( this.contains(i) ) {
            if ( ! a.contains(j)) throw new IllegalArgumentException("Sequence is not in the pile");
            off2 = getOffsetById(i)+p.getBestOffset2wrt1()-a.getOffsetById(j);
        } else {
            if ( this.contains(j)) {
                if ( ! a.contains(i)) throw new IllegalArgumentException("Sequence is not in the pile");
                off2 = getOffsetById(j)-p.getBestOffset2wrt1()-a.getOffsetById(i);
            } else throw new IllegalArgumentException("Sequence is not in the pile");
        }
        // stitch sequences from a into this pile:
        for ( Integer id : a ) {
            if ( this.contains(id) ) throw new IllegalArgumentException("Attempt to add a pile that shares sequences with the current one");
            add(a.getSequenceById(id),id,off2+a.getOffsetById(id));
        }
    }


    /** Adds another pile of aligned sequences (a) to this pile, stitching them together using specified
     * pairwise alignment p. Sequence ids must be set in the pairwise alignment, and one of those ids
     * must be in this pile, and the other in the pile 'a' being added, otherwise an IllegalArgumentException
     * is thrown. If pairwise alignment does not have sequence ids set, IllegalArgumentException is thrown.
     * The piles can not overlap.
     */
    public void add(MultipleAlignment a, PairwiseAlignment p) {
        if ( p.id1() == -1 || p.id2() == -1 ) throw new IllegalArgumentException("Attempt to add MSA based on pairwise alignemnt with sequence ids not properly set");
        add(a,p,p.id1(),p.id2());
    }

	/** Returns sequence associated with the specified external id, or null if sequence with this external id is
     * not found in the pile
	 * 
	 * @param id query id
	 * @return sequence for specified id or null
	 */
	public String getSequenceById(int id) {
		if ( ! contains(id)) return null;
		return seqs.get(index.get(id));
	}
	
	/** Returns offset relative to the first sequence in the pile for sequence associated with the specified
     * external id. If sequence with specified id is not found in the pile, RuntimeException is thrown.
	 * 
	 * @param id query id
	 * @return offset for sequence with specified id
	 */
	public int getOffsetById(int id) {
		if ( ! contains(id) ) throw new RuntimeException("Specified id is not in the pile");
		return alignment_offsets.get(index.get(id));
	}

    /** Returns external id of the read the offsets of this multiple alignment are based upon (i.e. all the offsets
     * are specified wrt the base read).
     * @return
     */
    public int getBaseReadId() { return ext_ids.get(0); }

    /** Returns offset of the read specified by its external id wrt the start of the consensus sequence in this
     * multiple alignment (consenus sequence is a major vote union of all the reads in this alignment).
     * @param id
     * @return
     */
    public int getOffsetWrtConsensus(int id) {
        return getOffsetById (id)- consensus.getStartOffset();
    }

	/** Returns true if the alignment already contains sequence with the specified id.
	 * 
	 * @param id
	 * @return
	 */
	public boolean contains(int id) {
		return index.containsKey(id);
	}

    /** Returns number of mismatches between sequences i and j (external ids) in the currently held multiple alignment.
     * Will return 0 if sequences do not overlap. Will throw RuntimeException if any of the specified ids is not
     * found in the current pile. 
     * @param i id of the first sequence
     * @param j id of the second sequence
     * @return mismatch count
     *
     * */
	public int countMismatches(int i, int j) {
		return PairwiseAlignment.countMismatches(getSequenceById(i), getSequenceById(j), getOffsetById(j)-getOffsetById(i));
	}

	/** Returns the length of the overlapping region of the two sequences specified by their external ids i and j.
	 * 
	 * @return overlap size
	 */
	public int getOverlap(int i, int j) {
		if ( ! contains(i) || ! contains(j)  ) throw new RuntimeException("Sequence with specified id is not in MSA pile");
		int off = getOffsetById(j) - getOffsetById(i);
        int L;
		if ( off >= 0 ) L = Math.min(getSequenceById(i).length()-off, getSequenceById(j).length());
		else L = Math.min(getSequenceById(j).length()+off, getSequenceById(i).length());
		return ( L < 0 ? 0 : L );
	}
	
	/** Given the two sequence ids, one of which has to be already in the pile, returns the one that is not in the pile.
	 * 
	 * @param i sequence id
	 * @param j sequence id
	 * @return one of the input arguments that is not found in the pile
	 * @throws IllegalArgumentException when either both or none of the specified indices are in the pile
	 */
	public int selectExternal(int i, int j) {
		if ( contains(i) ) {
			if ( contains(j) ) throw new IllegalArgumentException("Can not select external when both indices are in the pile");
			return j;
		} else {
			if ( ! contains(j) ) throw new IllegalArgumentException("Attempt to select external when both indices are not in the pile");
			return i;
		}
	}

    /** Returns a string consisting of n spaces.
     *
     * @param n
     * @return
     */
    private String skipN(int n) {
        StringBuilder b=new StringBuilder();
        for ( int k = 0 ; k < n ; k++ ) b.append(' ');
        return b.toString();
    }

    /** Prints n spaces directly into the specified string builder.
     *
     * @param n
     * @param b
     */
    private void skipN(int n, StringBuilder b) {
        for ( int k = 0 ; k < n ; k++ ) b.append(' ');
    }

	/** Returns a (multiline) string that represents the alignment visually: the sequences are appropriately
	 *  shifted and ready for printout;  
	 */
	public String toString(boolean inorder, boolean dotprint) {

		StringBuilder b = new StringBuilder();
		java.util.Formatter frmt = new java.util.Formatter(b);
		
		if ( seqs.size() == 0 ) return b.toString();
		
        final int first_offset = -consensus.getStartOffset();

        final int msa_length = consensus.length();
        char[][] consensusString = new char[4][msa_length];

        for ( int i = 0 ; i < msa_length ; i++ ) {

            Pair<Character,Integer> base = consensus.baseWithCountAt(i-first_offset);
            consensusString[3][i] = base.first;
            int mm = consensus.coverageAt(i-first_offset) - base.second;
            if ( mm > 0 ) {
                consensusString[2][i] = '*';
                if ( mm > 9 ) consensusString[0][i] = Character.forDigit(mm/10,10);
                else consensusString[0][i] = ' ';
                consensusString[1][i] = Character.forDigit(mm%10,10);
            } else {
                consensusString[0][i] = consensusString[1][i] = consensusString[2][i] = ' ';
            }
        }

        b.append("    "); b.append(consensusString[0]); b.append('\n');
        b.append("    "); b.append(consensusString[1]); b.append('\n');
        b.append("    "); b.append(consensusString[2]); b.append('\n');
        b.append("    "); b.append(consensusString[3]); b.append('\n');

        Integer[] perm = null;
        if ( inorder ) perm = Utils.SortPermutation(alignment_offsets);
		
		for ( int i = 0 ; i < seqs.size() ; i++ ) {
            int index = (inorder ? perm[i] : i);
			frmt.format("%3d:", ext_ids.get(index));
            int pos = alignment_offsets.get(index)+ first_offset; // start position on the consensus sequence
			skipN(pos,b);
            String aSeq = seqs.get(index);
            if ( dotprint ) {
                for ( int j = 0 ; j < aSeq.length() ; j++, pos++ ) {
                    if ( Character.toUpperCase(aSeq.charAt(j)) ==
                             Character.toUpperCase(consensusString[3][pos]) ) b.append('.');
                    else b.append(aSeq.charAt(j));
                }
            } else b.append(aSeq);
			b.append('\n');
		}
//		b.append(best_mm+" mismatches, "+ next_mm + " next best, " + getOverlap() + " overlapping bases, distance=" + distance() + "\n");
		return b.toString();
	}

    public String getConsensus() {
        return consensus.getSequence();
    }

    public String toString() { return toString(true, false); }

	public int size() { return seqs.size(); }
	
	/** Returns an iterator over the id's of the sequences currently stored in the pile 
	 * 
	 * @return
	 */
	public Iterator<Integer> sequenceIdIterator() { return index.keySet().iterator(); }

    /** Returns an iterator over external seuqnce ids of the sequences stored in the pile, presenting them in
     * the order of ascending alignment offsets.
     * @return
     */
    public Iterator<Integer> sequenceIdByOffsetIterator() {
        final Integer[] perm = Utils.SortPermutation(alignment_offsets);
        return new Iterator<Integer>() {
            private int i = 0;
            public boolean hasNext() {
                return i < perm.length;
            }
            public Integer next() {
                return ext_ids.get(perm[i++]);
            }
            public void remove() {
                throw new UnsupportedOperationException("remove not supported");
            }
        }   ;

    }

	public Iterator<Integer> iterator() {
		return sequenceIdIterator();
	}


}
