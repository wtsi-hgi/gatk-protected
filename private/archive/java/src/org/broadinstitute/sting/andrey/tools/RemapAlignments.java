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

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;

import java.io.File;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;


public class RemapAlignments extends CommandLineProgram {

    // Usage and parameters
    @Usage(programVersion="1.0") public String USAGE = "Remaps custom-reference (e.g. transcriptome) alignments onto the genomic reference\n";
	@Option(shortName="M", 
			doc="Map file: from the reference the reads were aligned to, to the master reference the alignments should be remapped to. "+
            "In other words, for each custom-reference contig C this map must provide a (possibly disjoint) list of intervals "+
            "on the target reference, onto which C maps base-by-base. ",
			optional=false)
	public File MAP_FILE = null;
	@Option(shortName="I", 
			doc="Input file (bam or sam) with alignments to be remapped",
			optional=false)
	public File IN = null;
	@Option(shortName="O", 
			doc="File to write remapped reads to.", 
			optional=false)
	public File OUT = null;
	@Option(shortName="R", 
			doc="Target reference to remap alignments onto.",
			optional=false)
	public File REFERENCE = null;
	@Option(
			doc="If a read has multiple alignments that are exactly the same after remapping, "+
			"then keep only one copy of such alignment in output file. Multiple alignments that are "+
			"not equivalent after remapping are not affected by this flag. "+
			"Multiple alignments for the same query must be grouped on adjacent lines of the input file to be detected "+
            "(i.e. input file must be sorted by read name), " +
			"otherwise REDUCE will have no effect.", 
			optional=true)
	public boolean REDUCE = false;


	private GenomicMap map = null;
	private String lastReadName = null;
	private int totalReads = 0;
	private int totalRecords = 0;
	private int badRecords = 0;
	private int totalUnmappedReads = 0;
	private int writtenRecords = 0;

	private Set<SAMRecord> remappedReads = null;
	private SAMFileWriter writer = null;
	private SAMFileReader reader = null;
	
	private static int [] g_log_n; // copied from bwa
	
	
    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new RemapAlignments().instanceMain(argv));
    }
    
    protected int doWork() {
    			
    	g_log_n = new int[256];
    	for (int i = 1; i < 256; ++i) g_log_n[i] = (int)(4.343 * Math.log(i) + 0.5);
    	
    	reader = new SAMFileReader(IN);
    	reader.setValidationStringency(ValidationStringency.SILENT);
		SAMFileHeader oldHeader = reader.getFileHeader();
		if ( oldHeader == null ) throw new RuntimeException("Failed to retrieve SAM file header from the input bam file");
		
		if ( REDUCE && oldHeader.getSortOrder() != SortOrder.queryname ) 
			System.out.println("WARNING: Input file is not sorted by query name, REDUCE may have no effect. Sort order: "
					+oldHeader.getSortOrder());
		
		remappedReads = new TreeSet<SAMRecord>(new AlignmentComparator());
		
		SAMFileHeader h = new SAMFileHeader();
		
		for ( Entry<String, String> attr : oldHeader.getAttributes() ) h.setAttribute(attr.getKey(), attr.getValue());
		h.setGroupOrder(oldHeader.getGroupOrder());
		h.setProgramRecords(oldHeader.getProgramRecords());
		h.setReadGroups(oldHeader.getReadGroups());
		
		if ( oldHeader.getSortOrder() == SortOrder.queryname ) {
			h.setSortOrder(SortOrder.queryname);
		} else {
			h.setSortOrder(SortOrder.unsorted);
		}
		
		ReferenceSequenceFileWalker reference = new ReferenceSequenceFileWalker(REFERENCE);

        if ( reference.getSequenceDictionary() == null ) {
        	System.out.println("No reference sequence dictionary found. Aborting.");
        	reader.close();
        	System.exit(1);
        }
		
		h.setSequenceDictionary(reference.getSequenceDictionary());
        GenomeLocParser genomeLocParser = new GenomeLocParser(reference.getSequenceDictionary());

		map = new GenomicMap(10000);
		map.read(genomeLocParser,MAP_FILE);
		System.out.println("Map loaded successfully: "+map.size()+" contigs");
				
		
		writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(h, true, OUT);
		
		for ( SAMRecord read : reader ) {
			
			
			if ( map.remapToMasterReference(read,h,true) == null ) {
				badRecords++;
				continue;
			}
			if ( AlignmentUtils.isReadUnmapped(read) ) totalUnmappedReads++;

            // destroy mate pair mapping information, if any (we will need to reconstitute pairs after remapping both ends):
            read.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            read.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
//				if ( read.getReadPairedFlag() ) System.out.println("PAIRED READ!!");

			totalRecords++;
			
			if ( totalRecords % 1000000 == 0 ) System.out.println(totalRecords + " valid records processed");
			

			if ( ! read.getReadName().equals(lastReadName) ) {
				totalReads++;
				lastReadName = read.getReadName();
			
						
				if ( REDUCE ) {
					
					updateCountsAndQuals(remappedReads);
					
					for ( SAMRecord r : remappedReads ) {
						writer.addAlignment(r); // emit non-redundant alignments for previous query
						writtenRecords++;
					}
					remappedReads.clear(); 
				}
			} 
			if ( REDUCE ) remappedReads.add(read); 
			else {
				writer.addAlignment(read);
				writtenRecords++;
			}
		}

		// write remaining bunch of reads:
		if ( REDUCE ) {
			updateCountsAndQuals(remappedReads);
			for ( SAMRecord r : remappedReads ) {
				writer.addAlignment(r); // emit non-redundant alignments for previous query
				writtenRecords++;
			}
		}
		
		System.out.println("Total valid records processed: "+totalRecords);
		System.out.println("Incorrect records (alignments across contig boundary) detected: "+badRecords + 
				" (discarded and excluded from any other stats)");
		System.out.println("Total reads processed: "+totalReads);
		System.out.println("Total mapped reads: "+(totalReads-totalUnmappedReads));
		System.out.println("Average hits per mapped read: "+((double)(totalRecords-totalUnmappedReads))/(totalReads-totalUnmappedReads));
		System.out.println("Records written: "+writtenRecords);
		System.out.println("Average hits per mapped read written (after reduction): "
				+((double)(writtenRecords-totalUnmappedReads))/(totalReads-totalUnmappedReads));
		reader.close();
		writer.close();
		return 0;
	}
	
    class AlignmentComparator implements Comparator<SAMRecord> {

    	public int compare(SAMRecord r1, SAMRecord r2) {
    		if ( r1.getReferenceIndex() < r2.getReferenceIndex() ) return -1; 
    		if ( r1.getReferenceIndex() > r2.getReferenceIndex() ) return  1;
    		if ( r1.getAlignmentStart() < r2.getAlignmentStart() ) return -1;
    		if ( r1.getAlignmentStart() > r2.getAlignmentStart() ) return 1;
    		return r1.getCigarString().compareTo(r2.getCigarString());
    	}
    	
    }

    private void updateCountsAndQuals(Set<SAMRecord> reads) {
    	if ( reads.size() == 1 ) {
    		SAMRecord r = reads.iterator().next();
 
        	// technically, if edit distance of the read is equal to max_diff used in alignments, 
        	// we should have set 25... 
                if ( AlignmentUtils.isReadUnmapped(r) ) {
                    r.setMappingQuality(0);
                } else {
                    r.setMappingQuality(37);
                    r.setAttribute("X0", Integer.valueOf(1));
                    r.setAttribute("X1", Integer.valueOf(0));
                }
    		r.setNotPrimaryAlignmentFlag(false);
    		
    	} else {
    		
    		// we have multiple alignments for the read
    		// need to figure out how many best vs inferior alignments are there:
    		int minNM = 1000000;
    		int cnt = 0; // count of best alignments
            Iterator<SAMRecord> it = reads.iterator();
            int n = reads.size(); // total number of (alternative) alignments for the given read.
            boolean canComputeMapQ = true;
    		while ( it.hasNext() ) {
                SAMRecord r = it.next();
                if ( AlignmentUtils.isReadUnmapped(r) && n > 1) {
                    // we do not want to keep unmapped records in the set unless it's the last and only record!
                    it.remove();
                    n--; // one less alignment left in the current group of alignments
                    continue;
                }
                if ( ! canComputeMapQ ) continue; // some reads were missing NM attribute, so do not bother - we can not compute MapQ
                Object attr = r.getAttribute("NM");
                if ( attr == null ) {
                    canComputeMapQ = false; // can not recompute qualities!
                    continue;
                } else {
    			    int nm;
                    if ( attr instanceof Short ) nm = ((Short)attr).intValue();
                    else if ( attr instanceof Integer ) nm = ((Integer)attr).intValue();
                    else throw new RuntimeException("NM attribute is neither Short nor Integer, don't know what to do.");
    			    if ( nm < minNM  ) {
    				    minNM = nm;
    				    cnt = 1;
    			    } else if ( nm == minNM ) cnt++;
                }
    		}

            if ( n == 1 ) {
                SAMRecord r = reads.iterator().next() ;
                if (AlignmentUtils.isReadUnmapped(r) ) {
                // special case: we are left with a single unmapped alignment
                    r.setAttribute("X0", new Integer(0));
                    r.setAttribute("X1", new Integer(0));
                    return;
                }
            }

    		// now reset counts of available alignments and mapping quals (if we can) in every alignment record:
    		for ( SAMRecord r : reads ) {
    			
    			int cnt2 = reads.size() - cnt; // count of inferior alignments
    			
    	   		r.setAttribute("X0", new Integer(cnt));   
        		r.setAttribute("X1", new Integer(cnt2));

                if ( ! canComputeMapQ ) continue; // not all reads had NM field, so we can not recompute MapQ

        		if ( cnt2 > 255 ) cnt2 = 255; // otherwise we will be out of bounds in g_log_n

                int nm_attr;
                Object attr =  r.getAttribute("NM");
                if ( attr instanceof Short ) nm_attr = ((Short)attr).intValue();
                else if ( attr instanceof Integer ) nm_attr = ((Integer)attr).intValue();
                else throw new RuntimeException("NM attribute is neither Short nor Integer, don't know what to do.");
    			if ( nm_attr == minNM ) { 
    				
    				// one of the best alignments:

    				r.setNotPrimaryAlignmentFlag(false);
    				if ( cnt == 1 ) {    					
    					// single best alignment; additional inferior alignments will only affect mapping qual
    					r.setMappingQuality( 23 < g_log_n[cnt2] ? 0 : 23 - g_log_n[cnt2] ); // this recipe for Q is copied from bwa
    				} else {
    					r.setMappingQuality(0); // multiple best alignments - mapping quality is 0
    				}
    			} else {
    				
    				// secondary alignment ( we know we hold a better one)
    				r.setNotPrimaryAlignmentFlag(true);
    				r.setMappingQuality(0); // ??? should we set 0 for secondary??
    			}
    		}
    	}
    	
    }
    
/*    
    private int bwa_approx_mapQ(SAMRecord r, int max_diff) {
    	int c1 = (Integer)r.getExtendedAttribute("X0");
    	int c2 = (Integer)r.getExtendedAttribute("X1");
    	int mm = (Integer)r.getExtendedAttribute("NM");
    	if ( c1 > 0 ) return 0;
    	if ( c1 == 0 ) return 23;
    	if ( mm == max_diff ) return 25;
    	return 0;
    }
*/
}


