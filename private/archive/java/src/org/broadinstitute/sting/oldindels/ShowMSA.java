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

import java.io.File;

import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class ShowMSA extends CommandLineProgram {

    // Usage and parameters
    @Usage(programVersion="1.0") public String USAGE = "Prints MSA into stdout\n";
    @Option(shortName="I", doc="SAM or BAM file with alignment data") public File INPUT_FILE;
    @Option(shortName="L", doc="Contig:Start-Stop or Contig:poslocation of the window to draw") public String LOCATION;
    @Option(shortName="W", doc="Number of bases on each side of specified position if LOCATION is in Contig:pos format; ignored otherwise", optional=true) public Integer WINDOW;
    @Option(shortName="R", doc="Reference fastb file") public File REF_FILE;
    @Option(shortName="P", doc="If true, then any read (partially) overlapping with the specified region will be shown. "+
    		"Otherwise (default), only reads fully contained in the specified interval are shown", optional=true) public Boolean PARTIAL;
    @Option(doc="Error counting mode: MM - count mismatches only, ERR - count errors (arachne style), MG - count mismatches and gaps as one error each") public String ERR_MODE;
    @Option(doc="Maximum number of errors allowed (see ERR_MODE)") public Integer MAX_ERRS;
    @Option(shortName="F",doc="Format: PILE - show alignment, FASTA - print sequences in fasta",optional=true) public String OUT_FORMAT;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new ShowMSA().instanceMain(argv));
    }
    
	protected int doWork() {
		
		if ( ! ERR_MODE.equals("MM") && ! ERR_MODE.equals("MG") && ! ERR_MODE.equals("ERR") ) {
			System.out.println("Unknown value specified for ERR_MODE");
			return 1;
		}

		if ( PARTIAL == null ) PARTIAL = new Boolean(false);
		if ( OUT_FORMAT == null ) OUT_FORMAT=new String("PILE");
		
		if ( ! OUT_FORMAT.equals("PILE") && ! OUT_FORMAT.equals("FASTA")) {
			System.out.println("OUT_FORMAT can only have values PILE or FASTA");
			return 1;
		}
		
		if ( ! INPUT_FILE.exists() ) {
			System.out.println("Specified INPUT_FILE does not exist");
			return 1;
		}

		if ( ! REF_FILE.exists() ) {
			System.out.println("Specified REF_FILE does not exist");
			return 1;
		}

		if ( LOCATION.indexOf(':') == -1 ) {
			System.out.println("LOCATION should follow Contig:Start-Stop or Contig:Pos format");
			return 1;
		}
		String[] s1 = LOCATION.split(":");
		int contig;
		try {
			contig = Integer.valueOf(s1[0]);
		}	catch (NumberFormatException e) {
			System.out.println("LOCATION: contig must be specified as an integer");
			return 1;
		}
		
		if ( s1.length != 2 ) {
			System.out.println("LOCATION should follow Contig:Start-Stop or Contig:Pos format");
			return 1;
		}
		
		String s2[] = s1[1].split("-");
		if ( s2.length > 2 ) {
			System.out.println("LOCATION should follow Contig:Start-Stop or Contig:Pos format");
			return 1;
		}
		int left, right;
		if ( s2.length == 2 ) {
			try {
				left = Integer.valueOf(s2[0]);
				right = Integer.valueOf(s2[1]);
			}	catch (NumberFormatException e) {
				System.out.println("LOCATION: window boundaries should be specified as integers");
				return 1;
			}
		} else {
			int pos = 0;
			try {
				pos = Integer.valueOf(s2[0]);
			}	catch (NumberFormatException e) {
				System.out.println("LOCATION: position on the contig should be specified as an integer");
				return 1;
			}
			if (WINDOW == null ) {
				System.out.println("WINDOW must be specified when LOCATION specifies a single poisiton (Contig:Pos)");
				return 1;
			}
			left = pos - WINDOW.intValue();
			right = pos+WINDOW.intValue();
		}
		
		
		String ref_contig ;
		
		try {
		    ReferenceSequenceFileWalker mRefReader =
                    new ReferenceSequenceFileWalker(ReferenceSequenceFileFactory.getReferenceSequenceFile(REF_FILE));
			ref_contig = mRefReader.get(contig).toString(); // reload ref
		} catch (Exception e) {
			System.out.println("Failed to read reference sequence from " + REF_FILE);
			return 1;
		}

		SAMFileReader reader ;
		try {
			reader = new SAMFileReader(INPUT_FILE);
		} catch ( Exception e) {
			System.out.println(e.getMessage());
			return 1;
		}
		
		SequencePile msa=null;
		
		if ( OUT_FORMAT.equals("PILE")) {
			msa = new SequencePile(ref_contig.substring(left-1, right));
		} else {
			System.out.println(">reference "+contig+":"+left+"-"+right);
			System.out.println(ref_contig.substring(left-1, right));
		}
		
		for( SAMRecord r : reader ) {
			if ( r.getReadUnmappedFlag() ) continue;
			if ( r.getReferenceIndex() < contig ) continue;
			if ( r.getReferenceIndex() > contig ) break;
			if ( r.getAlignmentEnd() < left ) continue;
			if ( r.getAlignmentStart() >= right ) break;
			if ( ! PARTIAL && ( r.getAlignmentStart() < left || r.getAlignmentEnd() >= right ) ) continue;
			
        	int err = -1;
        	if ( ERR_MODE.equals("MM")) err = numMismatches(r);
        	else if ( ERR_MODE.equals("ERR")) err = numErrors(r);
        	else if ( ERR_MODE.equals("MG")) err = numMismatchesGaps(r);
        	if ( err > MAX_ERRS ) continue;

        	if ( OUT_FORMAT.equals("PILE") ) {
        		msa.addAlignedSequence(r.getReadString(), r.getReadNegativeStrandFlag(), r.getCigar(), r.getAlignmentStart() - left);
        	} else {
        		System.out.print(">read "+r.getReadName());
        		if ( r.getReadNegativeStrandFlag() ) System.out.println("(rc)");
        		else System.out.println("(fw)");
        		System.out.println(r.getReadString());
        	}
		}
		
		if ( OUT_FORMAT.equals("PILE") ) msa.colorprint();
////			System.out.println(msa.format());
		
		return 0;
	}

	/** This method is a HACK: it is designed to work around the current bug in NM tags created  at CRD 
	 * 
	 * @param r SAM record that must specify an alignment
	 * @return number of errors (number of mismatches plus total length of all insertions/deletions
	 * @throws RuntimeException if cigar contains any elements other than M,I,D
	 */
	private static int numErrors(SAMRecord r) throws RuntimeException {
		
		// NM currently stores the total number of mismatches in all blocks + 1
		int errs = numMismatches(r);
		
		// now we have to add the total length of all indels:
		Cigar c = r.getCigar();
		for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
			CigarElement ce = c.getCigarElement(i);
			switch( ce.getOperator()) {
			case M : break; // we already have correct number of mismatches
			case I : 
			case D :
					errs += ce.getLength();
					break;
			default: throw new RuntimeException("Unrecognized cigar element");
			}
		}
		return errs;
	}

	/** This method is a HACK: it is designed to work around the current bug in NM tags created  at CRD 
	 * 
	 * @param r SAM record that must specify an alignment
	 * @return number of errors (number of mismatches plus total number of all insertions/deletions (each insertion or
	 * deletion will be counted as a single error regardless of the length)
	 * @throws RuntimeException if cigar contains any elements other than M,I,D
	 */
	private static int numMismatchesGaps(SAMRecord r) throws RuntimeException {
		
		// NM currently stores the total number of mismatches in all blocks + 1
		int errs = numMismatches(r);
		
		// now we have to add the total length of all indels:
		Cigar c = r.getCigar();
		for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
			CigarElement ce = c.getCigarElement(i);
			switch( ce.getOperator()) {
			case M : break; // we already have correct number of mismatches
			case I : 
			case D :
					errs++;
					break;
			default: throw new RuntimeException("Unrecognized cigar element");
			}
		}
		return errs;
	}
	
	
	/** This method is a HACK: it is designed to work around the current bug in NM tags created  at CRD */
	private static int numMismatches(SAMRecord r) throws RuntimeException {
		
		// NM currently stores the total number of mismatches in all blocks + 1
		return ((Integer)r.getAttribute("NM")).intValue() - 1;
		
	}
	
}
