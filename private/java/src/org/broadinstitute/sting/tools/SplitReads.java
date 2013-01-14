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

package org.broadinstitute.sting.tools;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.samtools.*;

import java.io.File;
import java.util.*;

public class SplitReads extends CommandLineProgram {
    @Usage(programVersion="1.0") public String USAGE = "Splits reads: extracts sub-sequences of the specified length(s) from left "+
            "and/or right ends of all the reads into the specified output bam file(s). For the reads in the input that are mapped, "+
            "the subsequences in the output bam(s) will have appropriately adjusted alignment positions and chopped cigars.";
    @Option(shortName="I",
            doc="Input file (bam or sam) with read sequences to split.",
            optional=false)
    public File IN = null;
    @Option(shortName="E", doc="Read end to select, 1=left, 2=right; default: select both ends.",
    		optional=true) public List<Integer> READ_ENDS = new ArrayList<Integer>();
    @Option(shortName="N", doc="Number of bases to keep in the corresponding segment of the read. "+
    			"Synchronized with READ_ENDS argument; if single number is given, all selected segments (ends) will have specified length.",
    		optional=false) public List<Integer> LENGTH = new ArrayList<Integer>();
    @Option(shortName="S", doc="Read name for each segment (read end) will be set as original read name followed by the corresponding suffix." +
    			"Synchronized with READ_ENDS argument and must have the same number of entries if specified (note that default READ_ENDS is a list of (1,2). "+
    			"By default, suffixes are empty strings, i.e. all segments have the same name(s) as the original read." , optional=true) public List<String> SUFFIXES = new ArrayList<String>();
    @Option(shortName="O",optional=false, doc="Each read end will be sent into the corresponding file " +
    		"(synchronized with READ_ENDS). If only one file name is specified, all read segments will be printed into that file."
    		) public List<File> OUTPUT_BAMS = new ArrayList<File>();
    @Option(shortName="U", doc="Split and output only unmapped reads; mapped reads will be ignored.",
    		optional=true) public boolean UNMAPPED = false;


    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new SplitReads().instanceMain(argv));
    }

    protected int doWork() {

    	// if read ends are not specified explicitly on the cmd line, set default 1,2 (both ends)
    	if ( READ_ENDS.size() == 0 ) {
    		READ_ENDS.add(1);
    		READ_ENDS.add(2);
    	}

    	for ( Integer i : READ_ENDS) {
    		if ( ! i.equals(1) && ! i.equals(2)) throw new RuntimeException("Unknown value specified for READ_ENDS: "+i);
    	}
    	
    	// if suffixes are not specified, set them to "", ""
    	if ( SUFFIXES.size() == 0 ) {
    		for ( Integer i : READ_ENDS) {
    			SUFFIXES.add( "" );
    		}
    	} else {
    		// or make sure that the number of suffixes matches the number of ends
    		if ( SUFFIXES.size() != READ_ENDS.size() ) throw new RuntimeException("Number of suffixes specified must be equal to the number of read ends requested."+
    				"Passed: "+ READ_ENDS.size() +" READ_ENDS and " + SUFFIXES.size() + " SUFFIXES arguments.");
    	}
    	
    	if ( LENGTH.size() == 1 ) {
    		// if only one length is specified, apply it to all ends:
    		LENGTH = Collections.nCopies(READ_ENDS.size(), LENGTH.get(0));
    	}
    	
    	if ( LENGTH.size() != READ_ENDS.size() ) throw new RuntimeException("Number of lengths specified must be equal to the number of read ends requested."+
				"Passed: "+ READ_ENDS.size() +" READ_ENDS and " + LENGTH.size() + " LENGTH arguments.");
    	
    	if ( READ_ENDS.size() != OUTPUT_BAMS.size() && OUTPUT_BAMS.size() != 1 )
            throw new RuntimeException("Number of output files must be either one, or equal to the number of read ends requested."+
				"Passed: "+ READ_ENDS.size() +" READ_ENDS and " + OUTPUT_BAMS.size() + " OUTPUT_BAMS arguments.");

        SAMFileReader inReader = new SAMFileReader(IN);

        List<SAMFileWriter> outWriters = new ArrayList<SAMFileWriter>(OUTPUT_BAMS.size());
        for ( File outName : OUTPUT_BAMS ) {
            outWriters.add(new SAMFileWriterFactory().makeSAMOrBAMWriter(inReader.getFileHeader(), true, outName)) ;
        }


        for ( SAMRecord read : inReader ) {

            if ( UNMAPPED && ! read.getReadUnmappedFlag() ) continue;
            
            for ( int i = 0 ; i < READ_ENDS.size(); i++ ) {

                SAMRecord newRecord = null;
                try {
                    newRecord = (SAMRecord)read.clone();
                } catch (CloneNotSupportedException e) {
                    throw new RuntimeException("Clone not supported by SAMRecord implementation");
                }

                final int whichEnd = READ_ENDS.get(i);
                final int length = LENGTH.get(i);
                String name = read.getReadName();
                if ( length > read.getReadLength() ) throw new RuntimeException("Read "+name+" is shorter than the specified length ("+read.getReadLength()+"<"+length+")");
                int start = 0 , stop = 0; // [start, stop) : segment of the read to be selected; coordinates are wrt read sequence; half-open 0 based
                switch ( whichEnd ) {
                case 1: start = 0 ; stop = start + LENGTH.get(i); break;
                case 2: stop = read.getReadLength() ; start = stop - LENGTH.get(i); break;
                }

                newRecord.setReadBases(Arrays.copyOfRange(read.getReadBases(),start,stop));
                newRecord.setBaseQualities(Arrays.copyOfRange(read.getBaseQualities(), start, stop));
                newRecord.setReadName(name+ SUFFIXES.get(i));
                if ( read.getReadUnmappedFlag() ) {
                    //newRecord.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                } else {
                    newRecord.setAlignmentStart(read.getAlignmentStart()+start);
                    newRecord.setCigar( chopCigar(read.getCigar(), start, length ));
                }

                if ( outWriters.size() > 1 ) outWriters.get(i).addAlignment(newRecord);
                else outWriters.get(0).addAlignment(newRecord);
            }

        }

        inReader.close();
        for ( SAMFileWriter w : outWriters ) w.close();

        return 0;
    }
    

	/**
	 * Returns new cigar representing segment of the alignment that starts at position <code>start</code> (0-based)
	 * with respect to the start of the original cigar and covers <code>length</code> bases on the original read the
	 * <code>origCigar</code> corresponds to (i.e. I elements count, but D do not).
	 * @param origCigar
	 * @param start
	 * @param length
	 * @return
	 */
	private Cigar chopCigar( Cigar origCigar, int start, int length ) {

		int elementEnd = 0; // next base after the end of the current cigar element on the read
		
		Cigar newCigar = new Cigar();
		
		Iterator<CigarElement> elements = origCigar.getCigarElements().iterator();

        if ( ! elements.hasNext() ) System.out.println("CIGAR HAS NO ELEMENTS!");

		CigarElement ce = null;
		
		while ( elementEnd <= start ) { // if we did not reach the start of selected segment yet:
//            System.out.println("INIT: start="+start+"; length="+length+"; elementEnd="+elementEnd);
			ce = elements.next();
			switch ( ce.getOperator() ) {
			case N:  //
			case D : // read misses bases wrt the ref, nothing to count on the read
				break;
			case I:
			case M:
            case EQ:
            case X:
            case S:
			case H: // all these elements are real bases on the read. Skip them completely if 
				    // 'start' is past them, or crop if it is inside:
				elementEnd += ce.getLength(); // 1 base past end of the current element on the read

			}
		}
		// at this point we are guaranteed that ce is the element that contains 'start' position;
		// now we start adding cigar elements:

		// add manually first element, since we need only a part of it after 'start':
		newCigar.add( new CigarElement(Math.min(elementEnd-start, length), ce.getOperator()) );

		int selectionEnd = start + length;
//		    System.out.println(origCigar.toString()+": start="+start+"; length="+length+"; selectionEnd="+selectionEnd+"; elementEnd="+elementEnd);
		while ( elementEnd < selectionEnd ) {
			ce = elements.next();
			switch ( ce.getOperator() ) {
			case N:  //
			case D : // read misses bases wrt the ref, nothing to count on the read, but the element has to be added:
				newCigar.add( new CigarElement(ce.getLength(), ce.getOperator()) );				
				break;
			case I:
			case M:
            case EQ:
            case X:
            case S:
			case H: // all these elements are real bases on the read. Add them and count them 
				    // making sure that the last element gets cropped if needed:
				elementEnd += ce.getLength(); // 1 base past end of the current element on the read
				if ( elementEnd > selectionEnd ) { // this is the last element we have to consider and it needs to be cropped:
					newCigar.add( new CigarElement(ce.getLength() -  elementEnd + selectionEnd , ce.getOperator()) );									
				} else {
					newCigar.add( new CigarElement(ce.getLength(), ce.getOperator()) );									
				}
			}
			
		}
		return newCigar;
		
	}
	
}
