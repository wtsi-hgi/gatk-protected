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

package org.broadinstitute.sting.gatk.walkers.andrey.tools;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.ArrayList;
import java.util.List;

/**
 * Iterates synchronously over two SAM files. At each iteration returs alignments with the same read name (in the order
 * the read names appear in the files). Alignment(s) from the first/second SAM file will be stored as the first/second
 * element of the pair, respectively. Multiple alignments (alternative placements) for a given read are allowed in both
 * input files. If only one of the files have alignment(s) for a given read name, the returned
 * pair will contain an empty list in the element corresponding to the other file. To enable this sort of traversal
 * synchronized by read names, the input SAM files must be sorted by read name. Constructor of this class verifies
 * that this is the case: SAM file  headers must report either 'queryname' sorting order, or no sorting order
 * (to allow the code to work with not-fully compliant 3rd party tools that do not set header flags properly; a warning
 * will be currently printed to stdout in this case); if sorting order in either of the files is set to "coordinate",
 * an exception will be thrown.
 */
    public class ParallelSAMIterator implements CloseableIterator< Pair< List<SAMRecord>, List<SAMRecord> > > {
        private SAMFileReader reader1;
        private SAMFileReader reader2;
        PushbackIterator<SAMRecord> i1;
        PushbackIterator<SAMRecord> i2;
        List<SAMRecord> alignments1;
        List<SAMRecord> alignments2;

        public ParallelSAMIterator(SAMFileReader r1, SAMFileReader r2) {
            reader1 = r1;
            reader2 = r2;
            checkSortOrder(r1,"End 1");
            checkSortOrder(r2, "End 2");
            i1 = new PushbackIterator(r1.iterator());
            i2 = new PushbackIterator(r2.iterator());
            alignments1 = nextGroup(i1); // pre-read next set of alignments
            alignments2 = nextGroup(i2);
        }


        /**
         * Returns <tt>true</tt> if the iteration has more elements. (In other
         * words, returns <tt>true</tt> if <tt>next</tt> would return an element
         * rather than throwing an exception.)
         *
         * @return <tt>true</tt> if the iterator has more elements.
         */
        public boolean hasNext() {
            return alignments1.size() > 0  || alignments2.size() > 0;
        }

        /**
         * Returns the next element in the iteration.
         *
         * @return the next element in the iteration.
         * @throws java.util.NoSuchElementException
         *          iteration has no more elements.
         */
        public Pair< List<SAMRecord>, List<SAMRecord> > next() {
            Pair< List<SAMRecord>, List<SAMRecord> > result;

            if ( alignments1.size() == 0 ) {
                // no more alignments left for end1
                result =  new Pair< List<SAMRecord>, List<SAMRecord> >(alignments1,alignments2);
                alignments2 = nextGroup(i2);
                return result;
            }
            if ( alignments2.size() == 0 ) {
                // no more alignments left for end2
                result =  new Pair< List<SAMRecord>, List<SAMRecord> >(alignments1,alignments2);
                alignments1 = nextGroup(i1);
                return result;
            }
            // next group of alignments is held for both ends. Check the read names:
            String end1Name = alignments1.get(0).getReadName();
            String end2Name = alignments2.get(0).getReadName();

            int cmp = end1Name.compareTo(end2Name);
            if ( cmp < 0 ) {
                // end1 goes before end2; return end1 with empty list for corresponding end2 and read next end1
                result = new Pair< List<SAMRecord>, List<SAMRecord> >(alignments1,new ArrayList<SAMRecord>());
                alignments1 = nextGroup(i1);
            } else {
                if ( cmp > 0 ) {
                    // end2 goes before end1; return end2 with empty list for corresponding end1 and read next end2
                    result = new Pair< List<SAMRecord>, List<SAMRecord> >(new ArrayList<SAMRecord>(),alignments2);
                    alignments2 = nextGroup(i2);
                } else {
                    // end 1 and end2 have the same read name => we got a mate pair:
                    result = new Pair< List<SAMRecord>, List<SAMRecord> >(alignments1, alignments2);
                    alignments1 = nextGroup(i1);
                    alignments2 = nextGroup(i2);
                }
            }
            return result;
        }

        /**
         * Removes from the underlying collection the last element returned by the
         * iterator (optional operation).  This method can be called only once per
         * call to <tt>next</tt>.  The behavior of an iterator is unspecified if
         * the underlying collection is modified while the iteration is in
         * progress in any way other than by calling this method.
         *
         * @throws UnsupportedOperationException if the <tt>remove</tt>
         *                                       operation is not supported by this Iterator.
         * @throws IllegalStateException         if the <tt>next</tt> method has not
         *                                       yet been called, or the <tt>remove</tt> method has already
         *                                       been called after the last call to the <tt>next</tt>
         *                                       method.
         */
        public void remove() {
            throw new UnsupportedOperationException("ParallelSAMIterator does not support remove() operation.");
        }

        public void close() {
            reader1.close();
            reader2.close();
        }

        /**
         * Read next alignment, and all immediately following ones that share same read name with the first;
         * return them all as a list.
         * @param i
         * @return
         */
        private List<SAMRecord> nextGroup(PushbackIterator<SAMRecord> i) {
            List<SAMRecord> result = new ArrayList<SAMRecord>();
            String readName ;

            if ( ! i.hasNext() ) return result; // nothing left
            SAMRecord r = i.next();
            readName = r.getReadName();
            result.add(r);

            while ( i.hasNext() ) {
                r = i.next();
                if ( ! r.getReadName().equals(readName) ) {
                    i.pushback(r);
                    break;
                }
                result.add(r);
            }
            return result;
        }

        /**
         * Utility method: checks that the sorting order in the input file is right
        *
        * @param reader sam file reader
        * @param fileName name of the file the reader is associated with. Used only to create more intelligible warning/exception messages,
        * you can actually pass any string here.
        */
        private void checkSortOrder(SAMFileReader reader, String fileName) {

            if  ( reader.getFileHeader() == null ) {
                System.out.println("WARNING: File "+fileName+" has no header. Assuming that file is sorted by read name.");
            }

            switch ( reader.getFileHeader().getSortOrder() ) {
            case coordinate:
               throw new RuntimeException("File "+fileName+" is sorted by coordinate. Sort it by read name first.");
            case unsorted:
               System.out.println("WARNING: file "+fileName+" has sorting order tag set to 'unsorted'. "+
                        "Assuming that it is sorted by read name.");
                break;
            case queryname: break; // good, that's what we need
                 default: throw new RuntimeException("File "+fileName + ": unknown sorting order ("+
                        reader.getFileHeader().getSortOrder()+")");
            }

        }

    }
