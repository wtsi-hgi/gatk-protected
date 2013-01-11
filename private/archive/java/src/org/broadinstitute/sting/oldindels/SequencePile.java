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

import java.util.List;
import java.util.ArrayList;
import net.sf.samtools.*;

public class SequencePile {
	private List<MSAColumn> mSeqGrid;
	private StringBuilder mRefGrid;
    private StringBuilder headerGrid;
	private int mDepth;
	private List<Boolean> mSeqRC;

	
	public SequencePile(String ref) {
		mRefGrid = new StringBuilder( ref );
        headerGrid = new StringBuilder();
        for ( int i = 0; i < ref.length(); i++ ) headerGrid.append(' ');
		mSeqGrid = new ArrayList<MSAColumn>();
		for ( int i = 0 ; i < mRefGrid.length(); i++ ) {
			mSeqGrid.add(new MSAColumn());
		}
		mDepth = 0;
		mSeqRC = new ArrayList<Boolean>();
	}
	
	/** Adds to the pile nucleotide sequence <seq> that aligns at zero-based position <refpos>
		relative to the original reference stretch the pile is built upon; the detailed alignment
		of the sequence to that reference stretch is specified by the <cigar>. 
		
		@param seq nucleotide sequence
		@param isRC true indicates that RC of the <seq> is being aligned
		@param cigar specification of the alignment of the sequence <seq> to the reference
		@param refpos 0-based position of the alignment with respect to the original stretch of the reference
			that was passed to the pile's constructor. Either <pos> or <pos>+sequence_length can be outside of 
			the pile's boundaries, the SequencePile class will deal with such situations correctly. 
		*/
	public void addAlignedSequence(String seq, boolean isRC, Cigar cigar, int refpos) {

		String alignedSeq = seq ;
//		if ( isRC ) {
//			alignedSeq = ReverseComplement(seq);
//		} else alignedSeq = seq;
		mSeqRC.add(isRC);

        // will hold actual position on the grid; reference can have insertions on the grid,
        // so position on the grid where we should start placing the read is not refpos!
		int pos = 0;
		for ( int i = 0 ; i < refpos ; i++ ) { // i is the position on the original reference
			// if we got some insertions on the reference prior to refpos, we need to count them in:
			while( mRefGrid.charAt(pos) == '+' ) {
				mSeqGrid.get(pos).add(' '); // add additional spaces in the line that will hold sequence seq
				pos++;
			}
			mSeqGrid.get(pos).add(' '); // fill with ' ' to the left of the read
			pos++;
		}
		
		// we reached start position of the alignment on the reference grid

		int readpos = 0; // position on the read
		
		for ( int i = 0 ; i < cigar.numCigarElements() ; i++ ) {

			final CigarElement ce = cigar.getCigarElement(i);
			
			switch(ce.getOperator()) {
    		case I: // read has an insertion
    				for ( int j = 0 ; j < ce.getLength() ; j++ ) {
						if ( pos >= mRefGrid.length() ) break;
						if ( pos >= 0 ) { 
							if ( mRefGrid.charAt(pos) !='+' ) {  // there was no insertion here yet: add it now!
								mRefGrid.insert(pos, '+');
                                headerGrid.insert(pos,'+');
								MSAColumn c = new MSAColumn();
                                // reads up to the previous depth (prior to adding current read) did not
                                // have an insertion here, so we insert '*' into all of them:
								for ( int k = 0 ; k < mDepth ; k++ ) {
									if ( mSeqGrid.get(pos-1).charAt(k) == ' ') c.add(' ');
									else c.add('*');
								}
								mSeqGrid.add(pos, c); // finally, add the base from the current read
							}
							mSeqGrid.get(pos).add(alignedSeq.charAt(readpos));
						}
						readpos++;
						pos++;
    				}
    				break;
    		case D: // read has a deletion
    				for ( int j = 0 ; j < ce.getLength() ; j++ ) {
						while( pos < mRefGrid.length() && mRefGrid.charAt(pos) == '+' ) { // skip insertions on the ref
							mSeqGrid.get(pos).add('*');
							pos++;
						}    					
						if ( pos >= mRefGrid.length() ) break;
						mSeqGrid.get(pos).add('-'); // mark deletion
                        headerGrid.setCharAt(pos,'-');
						pos++;
    				}
    				break;
    		case M: 
    				for ( int j = 0 ; j < ce.getLength() ; j++ ) {
						// if ref has an insertion, but the read does not: skip the insertion and continue with "gapless" alignment
						while( pos < mRefGrid.length() && mRefGrid.charAt(pos) == '+' ) {
							mSeqGrid.get(pos).add('*');
							pos++;
						}
						if ( pos >= mRefGrid.length() ) break;
						mSeqGrid.get(pos).add(alignedSeq.charAt(readpos));
                        if ( Character.toUpperCase(alignedSeq.charAt(readpos)) !=
                                Character.toUpperCase(mRefGrid.charAt(pos))
                                && headerGrid.charAt(pos)== ' ') headerGrid.setCharAt(pos,'*');
						pos++;
						readpos++;
    				}
    				break; 
    		default : throw new IllegalArgumentException("Unknown cigar element");
			}
		}
		for ( int i = pos ; i < mRefGrid.length() ; i++ ) { // i is the position on the modified reference
			mSeqGrid.get(i).add(' '); // fill with ' ' to the left of the read
		}
		mDepth++;
	}
	
	public String format() {
		StringBuffer b = new StringBuffer();
		b.append("  ");
		b.append(mRefGrid);
		b.append('\n');
		
		try {
		for ( int i = 0 ; i < mDepth; i++ ) {
			if ( mSeqRC.get(i).booleanValue() ) b.append("<-");
			else b.append("->");
			for ( int j = 0 ; j < mRefGrid.length() ; j++) {
				b.append(mSeqGrid.get(j).charAt(i));
			}
			b.append('\n');
		}
		} catch (Exception e) {}
		return b.toString();
	}
	
	private String ReverseComplement(String s) {
		StringBuffer b = new StringBuffer();
		char [] data = s.toCharArray();
		for ( int i = data.length - 1 ; i >= 0 ; i-- ) b.append(BaseComplement(data[i]));
		return b.toString();
	}
	
	private char BaseComplement(char b) {
		switch ( b ) {
		case 'A' : return 'T';
		case 'C': return 'G'; 
		case 'G': return 'C';
		case 'T': return 'A';
		default: throw new IllegalArgumentException(b + " is not a DNA base");
		}
	}

    public void colorprint() { colorprint(false); }

    public void dotprint(boolean printId) {

        String skip = null;
        if ( printId ) skip = new String("     ");
        else skip = new String("  ");

        System.out.print(formatHeader(skip));
        System.out.print(skip);
        System.out.println(mRefGrid);

        try {
        for ( int i = 0 ; i < mDepth; i++ ) {
            if ( printId ) System.out.printf("%3d",i);
            if ( mSeqRC.get(i).booleanValue() ) System.out.print("<-");
            else System.out.print("->");
            for ( int j = 0 ; j < mRefGrid.length() ; j++) {
                char seqbase = mSeqGrid.get(j).charAt(i);
                char refbase = mRefGrid.charAt(j);
                if ( isBase(refbase) && isBase(seqbase) &&
                        Character.toUpperCase(refbase) ==
                        Character.toUpperCase(seqbase) ) {
                    if ( mSeqRC.get(i) ) System.out.print(',');
                    else System.out.print('.');
                }
                else System.out.print(seqbase);
            }
            System.out.print('\n');
        }
        } catch (Exception e) {}
    }


	public void colorprint(boolean printId) {

        String skip = null;
        if ( printId ) skip = new String("     ");
		else skip = new String("  ");

        System.out.print(formatHeader(skip));
        System.out.print(skip);
        System.out.println(mRefGrid);
 
		try {
		for ( int i = 0 ; i < mDepth; i++ ) {
            if ( printId ) System.out.printf("%3d",i);
			if ( mSeqRC.get(i).booleanValue() ) System.out.print("<-");
			else System.out.print("->");
			for ( int j = 0 ; j < mRefGrid.length() ; j++) {
				char seqbase = mSeqGrid.get(j).charAt(i);
				char refbase = mRefGrid.charAt(j);
				if ( isBase(refbase) && isBase(seqbase) &&
                        Character.toUpperCase(refbase) !=
                                Character.toUpperCase(seqbase) ) System.out.print("\033[31m"+seqbase+"\033[30m");
				else System.out.print(seqbase);
			}
			System.out.print('\n');
		}
		} catch (Exception e) {}
	}

    private String formatHeader(String leadString) {
        char [][] mm_strings = new char[2][mRefGrid.length()];
        for ( int i = 0 ; i < mRefGrid.length() ; i++ ) {
            int count = 0;
            char refC = mRefGrid.charAt(i);
            MSAColumn col = mSeqGrid.get(i);
            if ( refC == '+' ) {
                 // count number of observations for insertion
                for ( int j = 0 ; j < col.size() ; j++ ) {
                     if ( col.charAt(j) != '*' && col.charAt(j) != ' ') count++;
                }
            } else {
                if ( headerGrid.charAt(i) == '-' ) {
                    // count number of observations for deletion
                    for ( int j = 0 ; j < col.size() ; j++ ) {
                         if ( col.charAt(j) == '-' ) count++;
                    }
                } else {
                    if ( headerGrid.charAt(i) == '*') {
                        for ( int j = 0 ; j < col.size() ; j++ ) {
                             if ( col.charAt(j)!=' ' &&
                                     Character.toUpperCase(col.charAt(j)) !=
                                     Character.toUpperCase(refC) ) count++;
                        }
                    }
                }
            }
            if ( count > 9 ) mm_strings[0][i] = Character.forDigit(count/10,10);
            else mm_strings[0][i] = ' ';
            if ( count > 0 ) mm_strings[1][i] = Character.forDigit(count%10,10);
            else mm_strings[1][i] = ' ';
        }

        StringBuilder b = new StringBuilder();
        b.append(leadString);
        b.append(mm_strings[0]);
        b.append('\n');
        b.append(leadString);
        b.append(mm_strings[1]);
        b.append('\n');
        b.append(leadString);
        b.append(headerGrid);
        b.append('\n');
        return b.toString();
    }

	private boolean isBase(char b) {
        b = Character.toUpperCase(b);
		return ( b=='A' ||b == 'C' || b=='G' || b=='T' || b=='N');
	}
}
