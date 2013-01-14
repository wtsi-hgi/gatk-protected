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

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class PairwiseAlignment {
		    private static final int IMPOSSIBLE = 1000000000;
			private String s1;
			private String s2;
			private int i1; // (external) id of the first sequence
			private int i2; // (external) id of the second sequence
			private int alignment_offset; // offset of s2 w/respect to s1
			private int best_mm; // mismatch count
			private int next_mm; // next-best mismatch count
			
			/** Initializes the alignment with pair of sequences (that will be immediately aligned) and
			 * stores their specified external ids id1, id2.
			 * @param is1 first nucleotide sequence (pre-indexed)
			 * @param is2 second nucleotide sequence (pre-indexed)
			 * @param id1 external id of the first sequence
			 * @param id2 external id of the second sequence
			 */
			public PairwiseAlignment(IndexedSequence is1, IndexedSequence is2, int id1, int id2 ) {
				s1 = new String(is1.getSequence());
				s2 = new String(is2.getSequence());
				i1 = id1;
				i2 = id2;
				best_mm = IMPOSSIBLE;
				next_mm = IMPOSSIBLE; 
				align(is1,is2);
			}
			
			/** Initializes the alignment with pair of sequences (that will be immediately aligned) and
			 * sets their external ids to -1. Such un-annotated pairwise alignment can not be added to MultipleAlignment.
			 *
			 */
			public PairwiseAlignment(IndexedSequence is1, IndexedSequence is2) {
				this(is1,is2,-1,-1);
			}
			
			/**
			 * Returns offset of sequence 2 with respect to sequence 1 in the best alignment
			 * @return positive offset if s2 is shifted right (starts later) wrt s1, or negative offset
			 *                 if s2 is shifted left (starts earlier) wrt s1
			 */
			public int getBestOffset2wrt1() { return alignment_offset; }

            /** Returns offset of the sequence j wrt sequence i in the best pairwise alignment found.
             *
             * @param i extrenal id of a sequence, must be one of the sequences kept by this alignment
             * @param j extrenal id of a sequence, must be one of the sequences kept by this alignment
             * @return offset of 2nd arg (j) wrt to the first arg (i)
             */
            public int getBestOffset2wrt1(int i, int j ) {
               if ( i == i1 && j == i2 ) return alignment_offset;
               else if ( i == i2 && j == i1 ) return -alignment_offset;
               throw new RuntimeException("Specified sequence id not found in the alignment");
            }

			public String getSequence1() { return s1; }
			public String getSequence2() { return s2; }
            public String getSequenceById(int i) {
                if ( i == i1 ) return s1;
                else if ( i == i2 ) return s2;
                throw new RuntimeException("Specified sequence id not found in the alignment");
            }
			public int id1() { return i1;}
			public int id2() { return i2;}
			
			/** Returns mismatch count in the best alignment found.
			 * 
			 * @return count of mismatches or impossibly large number of no mismatches were found
			 */
			public int getBestMMCount() { return best_mm; }
			
			/** Returns the number of mismatches in the next-best alignment found
			 * 
			 * @return next-best count of mismatches or impossibly large number if at most one alignment
			 * was ever found (that one would make the best then)
			 */
			public int getNextBestMMCount() { return next_mm; }
			
			/** Returns the length of the overlapping region of sequences s1 and s2 in the best alignment found, or -1 if
			 *   sequences do not align.
			 * 
			 * @return overlap size; can not be smaller than the size of the kmer used in IndexedSequence arguments the
			 * alignment was built from
			 */
			public int getOverlap() {
				if ( ! alignmentExists() ) return -1;
				if ( alignment_offset >= 0 ) {
					return Math.min(s1.length()-alignment_offset, s2.length());
				} else {
					return Math.min(s2.length()+alignment_offset, s1.length());
				}
			}

            public static int getOverlap(String seq1, String seq2, int offset2wrt1) {
                int L ;
                if ( offset2wrt1 >= 0 ) {
                    L = Math.min(seq1.length()-offset2wrt1, seq2.length());
                } else {
                    L = Math.min(seq2.length()+offset2wrt1, seq1.length());
                }
                return ( L < 0 ? 0 : L );
            }
			
			/** Returns true if at least one alignment, no matter how bad, was found between the two sequences
			 * (i.e. the sequences have at least one kmer in common).
			 */
			public boolean alignmentExists() { return best_mm < IMPOSSIBLE; }
			
			public void align(IndexedSequence is1, IndexedSequence is2) {
				
				Set<Integer> offsets = new HashSet<Integer>() ; // possible offsets of s2 wrt s1 as suggested by matching kmers
				for ( Map.Entry<Short,List<Integer>> e : is1 ) { // for each kmer in s1
					List<Integer> kmer_offsets_2 = is2.getOffsets(e.getKey());
					if ( kmer_offsets_2 == null ) continue; // uh-oh, kmer is not found in the other sequence
					for ( Integer i1 : e.getValue() ) {
						for ( Integer i2 : kmer_offsets_2 ) {
							offsets.add(i1-i2); // offset of seq 2 wrt seq1 as suggested by the currently inspected  occurences of the same kmer e.getKey() in both sequences 
						}
					}
				}
				// we have now a collection of distinct s1-s2 offsets seeded by matching kmers.
				// lets extend these kmer matches and count mismatches:
				
				for ( Integer trial_offset : offsets ) {
					int mm_cnt = countMismatches(is1.getSequence(), is2.getSequence(), trial_offset,next_mm+1);
//                    if ( (i1==4||i1==8) && i2==18) {
//                        if ( i1== 18 ) System.out.print("to " + i2+" : ");
//                        else  System.out.print("to " + i1+" : ");
//                        System.out.println("offset="+trial_offset.toString()+
//                                "; mm=" + countMismatches(is1.getSequence(),is2.getSequence(),trial_offset)+
//                                "(mm_cnt="+mm_cnt+")"+
//                                "; dist="+distance(is1.getSequence(),is2.getSequence(),trial_offset)+
//                                "; overlap="+getOverlap(is1.getSequence(),is2.getSequence(),trial_offset));
//                    }
                    // save current offset if alignment at this offset has fewer mismatches tham everything we've
                    // seen so far, or if it has same number of mismatches but has larger overlap (i.e. distance
                    // between sequences is smaller)
					if ( mm_cnt < best_mm ||
                           ( ( mm_cnt == best_mm ) &&
                            getOverlap(is1.getSequence(),is2.getSequence(),alignment_offset) <
                                    0.8*getOverlap(is1.getSequence(),is2.getSequence(),trial_offset) ) ) {
//                        if ( (i1==4||i1==8) && i2==18) System.out.println("Saved offset "+trial_offset.toString());
						alignment_offset = trial_offset;
 						next_mm = best_mm;
						best_mm = mm_cnt;
					} else {
						if ( mm_cnt < next_mm ) next_mm = mm_cnt;
					}
				}
			}
			
			public static int countMismatches(String seq1, String seq2, int offset2wrt1) {
				int pos1 = ( offset2wrt1 >= 0 ? offset2wrt1 : 0 );
				int pos2 = ( offset2wrt1 >= 0 ? 0 : -offset2wrt1 );  
				int cnt = 0;
				while ( pos1 < seq1.length() && pos2 < seq2.length() ) {
					if ( Character.toUpperCase(seq1.charAt(pos1++)) ==
                            Character.toUpperCase(seq2.charAt(pos2++)) ) continue;
					cnt++; // found mismatch						
				}
				return cnt;
			}

			public static int countMismatches(String seq1, String seq2, int offset2wrt1, int maxerr) {
				int pos1 = ( offset2wrt1 >= 0 ? offset2wrt1 : 0 );
				int pos2 = ( offset2wrt1 >= 0 ? 0 : -offset2wrt1 );  
				int cnt = 0;
				while ( pos1 < seq1.length() && pos2 < seq2.length() && cnt < maxerr ) {
					if ( Character.toUpperCase(seq1.charAt(pos1++)) ==
                            Character.toUpperCase(seq2.charAt(pos2++)) ) continue;
					cnt++; // found mismatch						
				}
				return cnt;
			}
			
			/** Returns a (multiline) string that represents the alignment visually: the sequences are appropriately
			 *  shifted and ready for printout; the pairwise alignment is followed by a stats line 
			 */
			public String toString() {
				StringBuffer b = new StringBuffer();
				int skip1 = ( alignment_offset >= 0 ? 0 : -alignment_offset );
				int skip2 = ( alignment_offset >=0 ? alignment_offset : 0 );
				for ( int k = 0 ; k < skip1 ; k++ ) b.append(' ');
				b.append(s1);
				b.append('\n');
				for ( int k = 0 ; k < skip2 ; k++ ) b.append(' ');
				b.append(s2);
				b.append('\n');
				b.append(best_mm+" mismatches, "+ next_mm + " next best, " + getOverlap() + " overlapping bases, distance=" + distance() + "\n");
				return b.toString();
			}

			public double distance() {
				int L = getOverlap();
				if ( L <=0 ) return 1e100;
				double l = ( best_mm==0? 1.0 : (double)best_mm + Math.sqrt((double)best_mm) );
				return ( l / (double)L );	
			}

            public static double distance(String seq1, String seq2, int offset2wrt1) {
                int L = getOverlap(seq1,seq2,offset2wrt1);
                if ( L <= 0 ) return 1e100;
                int mm = countMismatches(seq1,seq2,offset2wrt1);
                double l = ( mm == 0 ? 1.0 : (double)mm + Math.sqrt((double)mm) );
                return ( l / (double) L );
            }

}

	
