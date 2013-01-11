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
import java.util.Map;
import java.util.HashMap;


import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.picard.reference.ReferenceSequence;

import net.sf.samtools.*;
import org.broadinstitute.sting.utils.*;

public class IndelInspectorMain extends CommandLineProgram {

    // Usage and parameters
    @Usage(programVersion="1.0") public String USAGE = "Investigates indels called in the alignment data\n";
    @Option(shortName="I", doc="SAM or BAM file for calling",optional=true) public File INPUT_FILE;
    @Option(shortName="L",doc="Genomic interval to run on, as contig[:start[-stop]]; whole genome if not specified", optional=true) public String GENOME_LOCATION;
    @Option(shortName="V",doc="Verbosity level: SILENT, PILESUMMARY, ALIGNMENTS", optional=true) public String VERBOSITY_LEVEL;
    @Option(doc="Output file (sam or bam) for non-indel related reads and indel reads that were not improved (see OUTF)") public String OUT1; 
    @Option(doc="Output file (sam or bam) for improved (realigned) indel related reads") public String OUT2;
    @Option(doc="Output file (sam or bam) for indel related reads that fail to realign", optional = true ) public String OUTF;
    @Option(doc="[paranoid] If true, all reads that would be otherwise picked and processed by this tool will be saved, unmodified, into OUT1", optional=true) public Boolean CONTROL_RUN;
    @Option(doc="Error counting mode: MM - mismatches only (from sam tags), MC - mismatches only doing actual mismatch count on the fly (use this if tags are incorrectly set); ERR - errors (arachne style: mm+gap lengths), MG - count mismatches and gaps as one error each") public String ERR_MODE;
    @Option(doc="Maximum number of errors allowed (see ERR_MODE)") public Integer MAX_ERRS;
    @Option(shortName="R", doc="Reference fasta or fasta.gz file") public File REF_FILE;
    @Option(doc="Ignore reads that are longer than the specified cutoff (not a good way to do things but might be necessary because of performance issues)", optional=true) public Integer MAX_READ_LENGTH;
    @Option(doc="Realignment will be attempted around trains of indels with at least one indel observed COUNT_CUTOFF times or more",optional=true) public Integer COUNT_CUTOFF;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new IndelInspectorMain().instanceMain(argv));
    }
    
    protected int doWork() {

        int discarded_cigar_count = 0;
        int discarded_long_read_count = 0;
        int discarded_maxerr = 0;
        int reads_accepted = 0;
        int reads_with_indels_accepted = 0;

        ReferenceSequenceFileWalker reference = new ReferenceSequenceFileWalker(
                    REF_FILE
            );

        if ( reference.getSequenceDictionary() == null ) {
            System.out.println("No reference sequence dictionary found. Abort.");
        }

        GenomeLocParser.setupRefContigOrdering(reference.getSequenceDictionary());
        GenomeLoc location = null;
        if ( GENOME_LOCATION != null ) {
            location = GenomeLocParser.parseGenomeLoc(GENOME_LOCATION);
        }
        
        if ( COUNT_CUTOFF == null ) COUNT_CUTOFF = 2;
        
        if ( ! ERR_MODE.equals("MM") && ! ERR_MODE.equals("MG") && ! ERR_MODE.equals("ERR") && ! ERR_MODE.equals("MC") ) {
            System.out.println("Unknown value specified for ERR_MODE: "+ERR_MODE);
            return 1;
        }

        final SAMFileReader samReader = new SAMFileReader(getInputFile(INPUT_FILE,"/broad/1KG/"));
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        //        setContigOrdering(samReader);


        if ( MAX_READ_LENGTH == null ) MAX_READ_LENGTH = 1000000000;
         
        ReferenceSequence contig_seq = null;

        IndelRecordPileCollector col = null;
        PassThroughWriter ptWriter = new PassThroughWriter(OUT1,samReader.getFileHeader());
        PassThroughWriter ptFailedWriter = null;
        if ( OUTF != null ) ptFailedWriter = new PassThroughWriter(OUTF,samReader.getFileHeader());
        PileBuilder pileBuilder = null;
        if ( CONTROL_RUN == null ) CONTROL_RUN=false;
        if ( ! CONTROL_RUN ) pileBuilder = new PileBuilder(OUT2,samReader.getFileHeader(), ptFailedWriter == null? ptWriter : ptFailedWriter);

        try {
            if ( CONTROL_RUN ) col = new IndelRecordPileCollector(ptWriter, new DiscardingPileReceiver() );
            else col = new IndelRecordPileCollector(ptWriter, pileBuilder );
        } catch(Exception e) { System.err.println(e.getMessage()); }
        if ( col == null ) return 1; 

        col.setControlRun(CONTROL_RUN);
        col.setIndelCountAcceptanceCutoff(COUNT_CUTOFF);

        if ( ! CONTROL_RUN ) {
            if ( VERBOSITY_LEVEL == null ) VERBOSITY_LEVEL = new String("SILENT");
            if ( VERBOSITY_LEVEL.toUpperCase().equals("SILENT")) pileBuilder.setVerbosity(PileBuilder.SILENT);
            else if ( VERBOSITY_LEVEL.toUpperCase().equals("PILESUMMARY") ) pileBuilder.setVerbosity(PileBuilder.PILESUMMARY);
            else if ( VERBOSITY_LEVEL.toUpperCase().equals("ALIGNMENTS") ) pileBuilder.setVerbosity(PileBuilder.ALIGNMENTS);
            else {
                System.out.println("Unrecognized VERBOSITY_LEVEL setting.");
                return 1;
            }
        }

        String cur_contig = null;
       long t=0,tc=System.currentTimeMillis(); // time
        boolean done_printing = false;
        
        for ( SAMRecord r : samReader ) {

            if ( r.getReadUnmappedFlag() ) {  continue; }  
            if ( r.getReferenceName() != cur_contig) {
                cur_contig = r.getReferenceName();
                System.out.println("Contig "+cur_contig);
                // if contig is specified and we are past that contig, we are done:
                if ( location != null && GenomeLocParser.compareContigs(cur_contig, location.getContig()) == 1 ) break;
                if ( location == null || GenomeLocParser.compareContigs(cur_contig, location.getContig()) == 0 ) {
                	if ( location != null ) System.out.println("Time spent to scroll input bam file to the specified chromosome: "+ ((System.currentTimeMillis()-tc)/1000) + " seconds.");
                	tc = System.currentTimeMillis();
                    contig_seq = reference.get(r.getReferenceIndex());
                    t = System.currentTimeMillis();
                    String refstr = new String(contig_seq.getBases());
                    if (!CONTROL_RUN) pileBuilder.setReferenceSequence(refstr);
                    System.out.println("Contig "+cur_contig+" (index="+r.getReferenceIndex()+") loaded in "+ ((t-tc)/1000) +" seconds; length="+contig_seq.getBases().length+" tst="+contig_seq.toString());
                }
            }

            // if contig is specified and we did not reach it yet, skip the records until we reach that contig:
            if ( location != null && GenomeLocParser.compareContigs(cur_contig, location.getContig()) == -1 ) continue;

            if ( location != null && r.getAlignmentEnd() < location.getStart() ) continue;

            if ( location != null && ! done_printing ) {
            	System.out.println("Time spent to scroll input bam file to the specified location on the chromosome: " + ((System.currentTimeMillis()-t)/1000)+" seconds.");
            	done_printing = true;
            }
            // if stop position is specified and we are past that, stop reading:
            if ( location != null && r.getAlignmentStart() > location.getStop() ) break;

            //    if ( cur_contig.equals("chrM") || GenomeLoc.compareContigs(cur_contig,"chrY") > 0 ) continue; // skip chrM and unplaced contigs for now

            // we currently do not know how to deal with cigars containing elements other than M,I,D, so 
            // let's just skip the reads that contain those other elements (clipped reads?)
            Cigar c = r.getCigar();
            boolean cigar_acceptable = true;
            boolean has_indel = false;
            
            for ( int z = 0 ; z < c.numCigarElements() ; z++ ) {
                CigarElement ce = c.getCigarElement(z);
                switch ( ce.getOperator() ) {
                case M: break;
                case I:
                case D: has_indel = true; break;
                default: 
                	cigar_acceptable = false;
                }
            }
            if ( ! cigar_acceptable ) {
               	discarded_cigar_count++;
               	continue;
            }
            
            if ( r.getReadLength() > MAX_READ_LENGTH ) {
            		discarded_long_read_count++;
            		continue;
            }

            int err = -1;
/*
            System.out.println("MM:     "+numMismatches(r));
            System.out.println("direct: "+numMismatchesDirect(r,contig_seq));
            System.out.print("  ");
            for ( int i = r.getAlignmentStart() - 1 ; i < r.getAlignmentEnd() ; i++ ) System.out.print((char)contig_seq.getBases()[i]);
            System.out.println();
            System.out.println((r.getReadNegativeStrandFlag()?"<-":"->")+r.getReadString());
            System.out.println("cigar: "+r.getCigarString());
            System.out.println();
            if (counter++ == 20 ) break;
            continue;
*/

            if ( ERR_MODE.equals("MM")) err = numMismatches(r,contig_seq);
            else if ( ERR_MODE.equals("MC") ) err = AlignmentUtils.numMismatches(r,contig_seq);
            else if ( ERR_MODE.equals("ERR")) err = numErrors(r,contig_seq);
            else if ( ERR_MODE.equals("MG")) err = numMismatchesGaps(r,contig_seq);
            if ( err > MAX_ERRS.intValue() ) {
            	discarded_maxerr++;
            	continue;
            }
            
            reads_accepted++;
            if ( has_indel ) reads_with_indels_accepted++;
            //        	counter++;
            //        	if ( counter % 1000000 == 0 ) System.out.println(counter+" records; "+col.memStatsString());
            col.receive(r);

        }
        
        if ( ! CONTROL_RUN ) {
            pileBuilder.printStats();
            pileBuilder.close();
        }
        System.out.println("done.");
        System.out.println("Discarded reads with non-M,I,D cigar elements: "+ discarded_cigar_count);
        System.out.println("Discarded long reads (above "+MAX_READ_LENGTH+" bp): "+ discarded_long_read_count);
        System.out.println("Discarded reads with error counts above "+MAX_ERRS+ ": "+ discarded_maxerr);
        System.out.println("Reads passed to realigner: "+ reads_accepted+" total; "+reads_with_indels_accepted+" with indel(s)");
        System.out.println();
        col.printLengthHistograms();
        samReader.close();
        ptWriter.close();
        if ( ptFailedWriter != null ) ptFailedWriter.close();
        return 0;
    }
	
	/** This method is a HACK: it is designed to work around the current bug in NM tags created  at CRD 
	 * 
	 * @param r SAM record that must specify an alignment
	 * @return number of errors (number of mismatches plus total length of all insertions/deletions
	 * @throws RuntimeException
	 */
    private static int numErrors(SAMRecord r, ReferenceSequence refseq) throws RuntimeException {
		
		// NM currently stores the total number of mismatches in all blocks + 1
        int errs = numMismatches(r,refseq);
		
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
	 * @throws RuntimeException
	 */
    private static int numMismatchesGaps(SAMRecord r,ReferenceSequence refseq) throws RuntimeException {
		
		// NM currently stores the total number of mismatches in all blocks + 1
        int errs = numMismatches(r,refseq);
		
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
    private static int numMismatches(SAMRecord r, ReferenceSequence refseq) throws RuntimeException {
		
        // NM currently stores the total number of mismatches in all blocks + 1
        Integer i = (Integer)r.getAttribute("NM");
        if ( i == null ) return AlignmentUtils.numMismatches(r,refseq);
        return ((Integer)r.getAttribute("NM")).intValue() - 1;
		
    }

    /** Trivial utility method that goes some distance trying to ensure that the input file is there;
     * the only purpose is reducing clutter in main(). Receives a default
     * input file argument, does a few checks (e.g. that it is non-null and exists), if they fail tries
     * to fire up a file chooser dialog using start_folder as initial directory, etc.
     * @param default_arg some "default" input file; if it is non-null and exists, nothing else will be done,
     *        and the same default_arg objetc will be returned; otherwise the method will try to ask for a "better" input.
     * @param start_folder should file open dialog be fired up, it will initially display this directory.
     * @return File object that is not null and does exist (there is no check that it is a valid SAM/BAM file though).
     */
    private File getInputFile(File default_arg, String start_folder) {
        File f = default_arg;
        if ( f==null || ! f.exists() ) {
            JFileChooser fc = new JFileChooser(start_folder);
            FileNameExtensionFilter ff = new FileNameExtensionFilter("SAM and BAM files","sam","bam");
            fc.setFileFilter(ff);
            fc.setFileSelectionMode(JFileChooser.FILES_ONLY);

            int ret = fc.showOpenDialog(null);
            f = fc.getSelectedFile();
            if ( ret != JFileChooser.APPROVE_OPTION ) {
                System.out.println("No input file specified. Exiting...");
                System.exit(1);
            }
        }

        if ( f == null || ! f.exists() ) {
            System.out.println("SAM or BAM input file must be specified. Exiting...");
            System.exit(1);
        }

        return f;
    }

    /** Auxiliary method to remove some clutter from main(); gets called only once and tries to get
     * contig ordering from the header provided by opened SAM reader; if no header info is available
     * falls back to default ordering; whichever ordering is used, it is set for GenomeLoc class.
     * @param r sam reader to get header from
     */
    private void setContigOrdering(SAMFileReader r) {
        SAMFileHeader h = r.getFileHeader();
        if ( h == null ) {
            System.out.println("No header found in SAM file, falling back to default contig ordering");
            setDefaultContigOrdering();
            return;
        }
        GenomeLocParser.setupRefContigOrdering(h.getSequenceDictionary());
    }

    private void setDefaultContigOrdering() {
        Map<String,Integer> rco = new HashMap<String,Integer>();
        rco.put("chrM",0);
        for ( int i = 1 ; i <= 22 ; i++ ) rco.put(Integer.toString(i),i);//rco.put("chr"+i,i);
        rco.put("chrX",23);
        rco.put("chrY",24);
    }
}
