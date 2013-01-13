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

import org.broadinstitute.sting.utils.Pair;

import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 21, 2009
 * Time: 4:25:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class ConsensusSequence {
    private List< int[] > coverage; // counts of observations of every nucleotide at each position
    private long referencePos;// (arbitrary) reference position; when adding sequences, their offsets should be wrt this position
    private int startOffset; // offset of the leftmost base of this consensus sequence wrt referencePos; negative=left, positive=right
    private final static int NUMBINS = 4;
    private final static char[] BASES = { 'A','C','G','T' };

    public ConsensusSequence(int refPos) {
        coverage = new ArrayList< int[] >();
        referencePos = refPos;
        startOffset = 0;
    }

    public ConsensusSequence() {
        this(0);
    }

    /** Adds sequence <code>se</code> to the consensus, in the sense that all bases of seq are counted
     * into observed coverage kept by the consensus. The position of the sequence is specified by the
     * <code>offset</code> with respect to the fixed reference position of the consensus (the latter does not
     * have to be consensus start), and if the sequence extends beyound the consensus on either end, the
     * consensus will be extended appropriately to accomodate the full sequence.
     * @param seq nucleotide sequence ('ACGT...')
     * @param offset position of the start of the sequence relative to the fixed reference position of the consensus
     */
    public void addSequence(String seq, int offset) {
        // if sequence starts before than the currently held consensus oes, extend consensus to the left
        if ( offset < startOffset ) {
            coverage.addAll(0,instantiateCoverageList(startOffset-offset));
            startOffset = offset;
        }
        // if the sequence ends beyound the currently held consensus, extend consensus to the right
        if ( offset + seq.length() > startOffset + coverage.size() ) {
            coverage.addAll( instantiateCoverageList(offset+seq.length() - startOffset - coverage.size()) );
        }

        // count bases from the sequence into the coverage
        int posOnConsensus = offset - startOffset;
        for ( int i = 0 ; i < seq.length() ; i++, posOnConsensus++ ) {
            char base = Character.toUpperCase(seq.charAt(i));
            if ( base == 'N') continue;
            coverage.get(posOnConsensus)[baseToInt(base)]++;
        }
    }

    /** Removes sequence <code>seq</code> from the consensus. More exactly, 1 will be subtracted from current
     * observation counts kept by the consensus for each observed base at every position of the sequence. The
     * position of the sequence is specified by the <code>offset</code> with respect to the reference position
     * of the consensus. NOTE: this method is unchecked and does not verify that the sequence being subtracted
     * was indeed previously added to the consensus and/or that the consenus does accomodate full length of
     * the sequence. If it is not the case, the results can be unpredictable or assert failure may occur.
     *
     * @param seq nucleotide sequence ('ACGT...')
     * @param offset position of the start of the sequence relative to the fixed reference position of the consensus
     */
    public void removeSequence(String seq, int offset) {
        assert offset >= startOffset :
                "Attempt to remove from consensus a sequence that starts prior to consenus start";
        assert (offset+seq.length() < startOffset + coverage.size()) :
                "Attempt to remove from consensus a sequence that extends beyond consensus end";
        // subtract sequence bases from the coverage
        int posOnConsensus = offset - startOffset;
        for ( int i = 0 ; i < seq.length() ; i++, posOnConsensus++ ) {
            char base = Character.toUpperCase(seq.charAt(i));
            if ( base == 'N') continue;
            coverage.get(posOnConsensus)[ baseToInt(base) ]--;
        }
    }

    /** Returns offset of the start of consensus sequence with respect to the reference position the
     * consensus is pinned to.
     * @return
     */
    public int getStartOffset() { return startOffset; }

    /** Returns the length (number of bases) of the consensus sequence.
     *
     * @return
     */
    public int length() { return coverage.size(); }

    /** Returns the "distance" (score measuring the agreement) from the currently held consensus sequence to
     * the specified sequence <code>seq</code> starting at position <code>offset</code> wrt consenus reference position.
     * @param seq
     * @param offset
     * @return
     */
    public double distance(String seq, int offset) {
        int posOnConsensus; //  index into the currently held consensus sequence
        int i ; // index into the passed sequence argument
        if ( offset < startOffset ) {
            posOnConsensus = 0;
            i = startOffset - offset;
        } else {
            i = 0 ;
            posOnConsensus = offset - startOffset;
        }
        // stop position on the passed sequence (can be less than sequence length if consensus stops prematurely)
        int stop = Math.min(offset+seq.length(), startOffset+coverage.size() ) - offset;

        for (  ; i < stop ; i++, posOnConsensus++ ) {
            int base = baseToInt(Character.toUpperCase(seq.charAt(posOnConsensus)));
            int [] cov = coverage.get(posOnConsensus);
            int totalcov = cov[0]+cov[1]+cov[2]+cov[3];

        }
        return 0.0;
    }

    /** Returns consensus base at the specified offset wrt the consesus sequence's reference position.
     * Specified offset must be within the span of currently held consensus sequence. Consensus base is the
     * one with the maximum count of observations. If two different nucleotides were observed exactly the
     * same number of times (and that number is greater than the number of observations for othe nucleotides),
     * the "lesser" one, (order being ACGT) will be returned. If coverage at specified position is zero, 'N' will
     * be returned.
     * @param offset
     * @return
     */
    public char baseAt(int offset) {
        assert offset >= startOffset && offset < startOffset + coverage.size() : "Offset out of bounds";
        int [] cov = coverage.get(offset-startOffset);
        int total_cov = cov[0] + cov[1] + cov[2] + cov[3];
        int bmax = 0;
        char base = 'N';
        for ( int z = 0; z < 4 ; z++ ) {
            if ( cov[z] > bmax ) {
                bmax = cov[z];
                base = BASES[z];
            }
        }
        return base;        
    }

    /** Returns consensus base at the specified offset together with its observation count.
     *
     * @param offset
     * @return
     * @see #baseAt(int)
     */
    public Pair<Character,Integer> baseWithCountAt(int offset) {
        assert offset >= startOffset && offset < startOffset + coverage.size() : "Offset out of bounds";
        int [] cov = coverage.get(offset-startOffset);
        int total_cov = cov[0] + cov[1] + cov[2] + cov[3];
        int bmax = 0;
        char base = 'N';
        for ( int z = 0; z < 4 ; z++ ) {
            if ( cov[z] > bmax ) {
                bmax = cov[z];
                base = BASES[z];
            }
        }
        return new Pair<Character,Integer>(base,bmax);
    }

    /** Returns total coverage (all observations regardless of what base what observed) at position
     * specified by offset with respect to the conensus' reference position. offset does not have to be within
     * the bounds of the currently kept consensus sequence, if it falls outside, a 0 will be silently returned.
     * @param offset
     * @return
     */
    public int coverageAt(int offset) {
        if ( offset < startOffset || offset >= startOffset + coverage.size() ) return 0;
        int [] cov = coverage.get(offset-startOffset);
        return cov[0]+cov[1]+cov[2]+cov[3];
    }

    /** Returns consesus sequence as a astring of bases (ACGTN); N will be returned for positions with zero
     * coverage.
     * @return
     */
    public String getSequence() {
        char [] b = new char[coverage.size()];
        for ( int i = 0 ; i < b.length ; i++ ) {
            b[i] = baseAt(i+startOffset);
        }
        return new String(b);
    }

    private List<int[]> instantiateCoverageList(int n) {
        List< int[] > subseq = new ArrayList<int[] >(n);
        for ( int i = 0 ; i < n ; i++ ) subseq.add(new int[NUMBINS]);
        return subseq;
    }

    private int baseToInt(char c) {
        int base;
        switch( Character.toUpperCase(c) ) {
            case 'A': base = 0; break;
            case 'C': base = 1; break;
            case 'G': base = 2; break;
            case 'T': base = 3; break;
            case 'N': base = -1; break;
            default : throw new IllegalArgumentException("Sequence can contain only ACGTN symbols");
        }
        return base;
    }
}
