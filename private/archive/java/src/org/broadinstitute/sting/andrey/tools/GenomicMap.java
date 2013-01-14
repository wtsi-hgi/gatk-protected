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

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.variant.utils.BaseUtils;

import java.io.*;
import java.util.*;

public class GenomicMap implements Iterable<Map.Entry<String, Collection<GenomeLoc> > >{
	
	private Map<String,Collection<GenomeLoc> > map;
	
	/** Creates new empty genomic map preallocated to handle <code>initialContig</code> contigs.
	 * 
	 * @param initialContigs
	 */
	public GenomicMap(int initialContigs) {
		map = new HashMap<String,Collection<GenomeLoc> >(initialContigs);
	}
	
	/** Creates new empty genomic map */
	public GenomicMap() {
		this(1000);
	}

	/** Adds custom contig to the map, as a collection of intervals on the master reference.
	 * 
	 * @param name name of the custom contig; can not be null
	 * @param c mapping of the custom contig sequence onto intervals on the master reference
	 */
	public void addCustomContig(String name, Collection<GenomeLoc> c) {
		if ( name == null ) throw new ReviewedStingException("Custom contig name can not be null");
		if ( map.containsKey(name)) throw new ReviewedStingException("Custom contig "+name+" already exists");
		map.put(name, c);
	}
	
	/** Returns mapping of the specified custom contig onto the custom reference. \
	 * If no such contig exists, returns null.
	 * 
	 * @param name
	 * @return
	 */
	public Collection<GenomeLoc> getContigMapping(String name) { return map.get(name); }
	
	/** Read genomic map from specified Arachne multimap file. Format: 
	 * contig_id start stop contig_id start stop ... # name ...
	 * where start, stop are 0 based, closed intervals  
	 * @param f
	 */
	public void readArachne(SAMSequenceDictionary sequenceDictionary,GenomeLocParser genomeLocParser,File f) {
		
		try {
			BufferedReader reader = new BufferedReader( new FileReader(f) );
			
			String line = null;
			while( ( line = reader.readLine() ) != null ) {
				String[] halves = line.split("#",2);
				if ( halves.length < 2 ) 
					throw new UserException.MalformedFile(f, "Line: "+line+"\nin map file "+f+"\n does not contain contig name");
				
				int p1 = 0;
				for ( ;  p1 < halves[1].length() && Character.isWhitespace(halves[1].charAt(p1) ); p1++ ); 
				// p1 is now index of first non-space
				int p2 = p1;
				for ( ; p2 < halves[1].length() && ! Character.isWhitespace(halves[1].charAt(p2) ); p2++ );
				// p2 is index of first whitespace after first word
				
				if ( p1 == p2 ) 
					throw new UserException.MalformedFile(f, "Line: "+line+"\n in map file "+f+"\nNo contig name found after '#'");
				
				String name = halves[1].substring(p1, p2);
								
				String[] coord_parts = halves[0].split("\\s");
				if ( coord_parts.length % 3 != 0 ) 
					throw new UserException.MalformedFile(f, "Line: "+line+"\n in map file "+f+"\nNumber of coordinate fields is not a multiple of 3");
				
				List<GenomeLoc> segments = new ArrayList<GenomeLoc>( coord_parts.length / 3 );
				
				for ( int i = 0 ; i < coord_parts.length ; i += 3 ) {
					// Arachne map file contains 0-based, closed intervals, hence +1 below.
					int index = Integer.parseInt(coord_parts[i]);
                    String contig = sequenceDictionary.getSequence(index).getSequenceName();
					int start = Integer.parseInt(coord_parts[i+1]);
					int stop = Integer.parseInt(coord_parts[i+2]);
					segments.add(genomeLocParser.createGenomeLoc(contig, start+1, stop+1));
				}
				
				addCustomContig(name, segments);
				
			}
			reader.close();
		} catch ( FileNotFoundException e) {
			throw new UserException.CouldNotReadInputFile(f, e);
		} catch (IOException e) {
			throw new UserException.CouldNotReadInputFile(f, e);
		}
	}

	/** Read genomic map from specified file in "new" format. Format: 
	 * name chr:start-stop,chr:start-stop,...,chr:start-stop
	 * where start, stop are 1 based, closed intervals  
	 * @param f
	 */
	public void read(GenomeLocParser genomeLocParser,File f) {
		
		try {
			BufferedReader reader = new BufferedReader( new FileReader(f) );

			String line = null;
			while( ( line = reader.readLine() ) != null ) {
				int p1 = 0;
				while ( p1 < line.length() && Character.isWhitespace(line.charAt(p1))) p1++;
				int p2 = p1;
				while ( p2 < line.length() && ! Character.isWhitespace(line.charAt(p2))) p2++;
				if ( p1 == p2 ) continue; // empty line

				String name = line.substring(p1, p2);
				
				List<GenomeLoc> segments = new ArrayList<GenomeLoc>( 5 );

				p1 = p2+1; // set p1 after first whitespace after the name
				while ( p1 < line.length() && Character.isWhitespace(line.charAt(p1))) p1++; // skip whitespaces
				p2 = p1;
				while ( p2 < line.length() && line.charAt(p2) != ',') p2++; // next comma or end-of-line

				while ( p2 != p1 ) {
					GenomeLoc newSegment = genomeLocParser.parseGenomeLoc(line.substring(p1, p2));
					if ( segments.size() > 0 &&
							segments.get(segments.size()-1).getStop()+1 == newSegment.getStart() &&
							segments.get(segments.size()-1).getContigIndex() == newSegment.getContigIndex())
						System.out.println("WARNING: strictly adjacent segments found in custom contig "+name);

					segments.add(newSegment);
				
					p1 = p2+1; // set p1 after the comma
					while ( p1 < line.length() && Character.isWhitespace(line.charAt(p1))) p1++; // skip whitespaces
					p2 = p1;
					while ( p2 < line.length() && line.charAt(p2) != ',') p2++; // next comma or end-of-line
				}
				if ( segments.size() == 0 ) throw new ReviewedStingException("Line "+line+" has no intervals specified");
				addCustomContig(name, segments);
			}
			reader.close();
		} catch ( FileNotFoundException e) {
			throw new UserException.CouldNotReadInputFile(f, e);
		} catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(f, e);
		}
	}

	public void write(File f) {
		try {
			BufferedWriter writer = new BufferedWriter( new FileWriter( f ));
			for ( String name : nameSet() ) {
				writer.append(name+" ");
				Iterator<GenomeLoc> iter = getContigMapping(name).iterator();
				if ( iter.hasNext() ) writer.append(iter.next().toString());
				while (iter.hasNext()) {
					writer.append(',');
					writer.append(iter.next().toString());
				}
				writer.append('\n');
			}
			writer.close();
		} catch (IOException e) {
			throw new UserException.CouldNotCreateOutputFile(f, e);
		}
	}
	
	/** Remaps a record (read) aligned to a custom contig back onto the master reference. 
	 * If the map does not have mapping information for
	 * the contig, an exception will be thrown. This method changes read's reference name, start position and 
	 * cigar, as well as the read's file header (must be provided). 
	 * 
	 * Some aligners (e.g. bwa) can return "alignments" spanning across contig boundaries. The last argument of this
	 * method controls the behavior in this case: if it is set to true, such alignments are ignored upon detection,
	 * and the method returns null. Otherwise, strict validation mode is used: if aligned read extends beyond the
	 * contig boundary, an exception is thrown.
	 *  
	 * @param r read, alignment information (contig, start position, cigar) will be modified by this method
	 * @param h SAM file header for the master reference the alignment is being mapped onto; will be substituted for the read's header.
	 * @return same read instance that was passed to this method, remapped
	 */
	public SAMRecord remapToMasterReference(SAMRecord r, SAMFileHeader h, boolean discardCrossContig) {
		if ( AlignmentUtils.isReadUnmapped(r) ) {
            // set to NO_... just in case: in principle, SAM format spec allows unmapped reads (with 'unmapped'
            // flag raised) to have reference contig and start position set to arbitrary values for sorting
            // purposes; after remapping, these values would make no sense or even cause a crash when reading
            // remapped bam
            r.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            r.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);

            // these are required because current santools jdk is over-validating and requiring MAPQ, CIGAR etc be
            // set to 0/null for unmapped reads. In principle, it should not matter.
            r.setMappingQuality(0);
            r.setCigar(new Cigar());
            r.setNotPrimaryAlignmentFlag(false);
            if ( r.getReadNegativeStrandFlag() ) {
                r.setReadBases(BaseUtils.simpleReverseComplement(r.getReadBases()));
                r.setBaseQualities(Utils.reverse(r.getBaseQualities()));
                r.setReadNegativeStrandFlag(false);
            }
            return r; // nothing to do if read is unmapped
        }
		
		int customStart = r.getAlignmentStart();

        // get mapping from read's contig onto a "global" contig (as a list of intervals on the latter):
		Collection<GenomeLoc> segments = getContigMapping(r.getReferenceName());
		if ( segments == null ) throw new UserException.MalformedBAM(r, "Can not remap a record: unknown custom contig name "+r.getReferenceName());

        // scroll the list of intervals until we find the interval that the alignment start falls into:
		Pair<? extends Iterator<GenomeLoc>, Integer> p = seekForward(segments,customStart);
		
		Iterator<GenomeLoc> iter = p.first;
		
		GenomeLoc gl = iter.next(); // initialization: get interval that contains customStart

        // p.second is 1-based position of alignment start relative to interval gl;
        // hence refPos is 1-based position of the alignment start on the master ref (since gl.getStart() is 1-based too)
		int refPos = (int)(p.second+gl.getStart()-1); 
		
		String oldRefName = r.getReferenceName();
		int oldStart = r.getAlignmentStart();
		int oldEnd = r.getAlignmentEnd();

		r.setAlignmentStart(refPos);		
		
		r.setHeader(h); // have to substitute here, or setReferenceIndex will not work correctly below
		
		r.setReferenceIndex(gl.getContigIndex());
		
		Cigar oldCigar = r.getCigar();
		Cigar newCigar = new Cigar();
		int N = oldCigar.numCigarElements() ;
		
		long currStop = gl.getStop();// end of the current segment of the custom contig on the master reference, 1-based inclusive
        int delayedGap = 0 ; // length of the 'N' gap between the segments (intervals) on the master ref, to be added only if followed by another cigar element		
		for ( int k = 0; k < N ; k++ ) {
			CigarElement ce = oldCigar.getCigarElement(k);
			int len = ce.getLength();      // length of the cigar element
			switch( ce.getOperator() ) {
			case S: // soft clip
			case H: // or hard clip - these are not included in getAlignmentStart, so pass them through
				if ( k != 0 && k != N-1 ) // paranoid
					throw new ReviewedStingException("Don't know what to do with S or N cigar element that is not at the either end of the cigar. Cigar: "+
							r.getCigarString());
			case I: // insertions are passed through as well
				newCigar.add(new CigarElement(len,ce.getOperator()));
				break;
			case D:
            case M:
            case EQ:
            case X:
///////////
                if ( delayedGap > 0 ) {
                    // we get here if previous M or D element ended exactly at the interval boundary; we need
                    // to add the stretch of N's only if that element turned out to be not the last one, so we do it now
                    newCigar.add(new CigarElement(delayedGap, CigarOperator.N));
                    delayedGap = 0;
                }
                while ( refPos + len - 1 > currStop ) { // current D or M cigar element extends beyond the end of current segment

                    // we have that many bases in the current cigar element till the end of the current segment:
                    int currLength = (int)(currStop-refPos+1);


                    // curr length can be exactly 0 if previous element ended exactly at the segment boundary:
                    // after that element was processed, refPos was set to currStop+1, so in this special case we need
                    // *first* to switch to next segment, *then* start adding bases from the current element.
                    if ( currLength > 0 ) {
                        newCigar.add(new CigarElement( currLength,ce.getOperator()) ); // record deletion/match till the end of the current segment
                        len -= currLength; // we still have 'len' bases remaining in the current cigar element
                    }

                    // NOTE: since we entered the loop, we were guaranteed that len > currLength, so now len > 0

                    // check if we have next segment to extend remaining matching bases to; if we don't, something's awfully wrong:
                    if ( ! iter.hasNext() ) {
                        String message = "Record "+r.getReadName()+" extends beyond its custom contig."+
                                        "\nRead aligns to: "+oldRefName+":"+oldStart+"-"+oldEnd+"; cigar="+
                                        r.getCigarString()+"; contig length="+contigLength(segments);
                        if ( discardCrossContig ) {
                //			System.out.println("WARNING: ALIGNMENT DISCARDED: "+message);
                            return null;
                        } else throw new UserException.MalformedBAM(r, message);
                    }

                    gl = iter.next(); // advance to next segment

                    refPos = (int)gl.getStart(); // we jump to the start of next segment on the master ref

                    if ( gl.getContigIndex() != r.getReferenceIndex() )
                        throw new UserException.MalformedBAM(r, "Contig "+oldRefName+
                        " has segments on different master contigs: currently unsupported");

                    if ( refPos < currStop + 1 )
                        throw new UserException.MalformedBAM(r, "Contig "+oldRefName+
                        " has segments that are out of order or strictly adjacent: currently unsupported");
                    if ( len > 0 && refPos > currStop + 1 ) {
                        // add "panning" N's w/respect to the master ref over the region between adjacent segments
                        // (and do not add anything if segments are strictly adjacent, i.e. refPos == currStop+1):
                        newCigar.add(new CigarElement((int)(refPos-currStop-1),CigarOperator.N));
                    } else {
                        // we jumped onto the next interval, but the current cigar element ended exactly
                        // at the end of the previous interval. We will need to end later, only if more M/D elements follow:
                        delayedGap = (int)(refPos-currStop-1);
                    }
                    currStop = gl.getStop();
                    // now we can continue with recording remaining matching bases over the current segment
                }
                // we get here when remaining matching bases fit completely inside the current segment:
                if ( len > 0 ) newCigar.add(new CigarElement(len,ce.getOperator()));
                refPos+=len;

                break;
////////////
			}
		}
		
		r.setCigar(newCigar);
		
		return r;
	}
	
	public int size() { return map.size(); }
	
	public Iterator<Map.Entry<String, Collection<GenomeLoc> > > iterator() { return map.entrySet().iterator(); }
	public Iterator<String> nameIterator() { return map.keySet().iterator(); }
	public Set<String> nameSet() { return map.keySet(); }
	
	/** Returns an iterator into the specified collection of segments that points right before the segment that contains
	 * specified position, and the offset of the position inside that segment. This helper method assumes that
	 * there is a "custom" contig built of intervals on the "master" reference; the first argument specifies
	 * the mapping (i.e. an ordered collection of master reference intervals the custom contig is built of), and the second argument
	 * is the 1-based position on that custom contig. Returned iterator is advanced towards the interval (element of the passed
	 * collection) that contains the specified position, namely a call to next() on the returned iterator will return that interval.
	 * Returned integer offset is the 1-based offset of the base at position <code>position</code> on the custom contig with respect
	 * to the start of the interval that base. If position is outside of the custom contig, runtime StingException will be thrown. 
	 * @param segments mapping of the custom contig onto the master reference
	 * @param position 1-based position on the custom contig
	 * @return
	 */
	private Pair<PushbackIterator<GenomeLoc>,Integer> seekForward(Collection<GenomeLoc> segments,int position) {
		
		if ( position < 1 ) throw new ReviewedStingException("Position "+position + " is outside of custom contig boundaries");
		
		PushbackIterator<GenomeLoc> iter = new PushbackIterator<GenomeLoc>(segments.iterator());
		
		while ( iter.hasNext() ) {
			GenomeLoc current = iter.next();
			long length = current.getStop() - current.getStart() + 1; // length of current segment
			if ( position <= length ) { // position is on the current segment
				iter.pushback(current);
				return new Pair<PushbackIterator<GenomeLoc>, Integer >( iter,position);
			}
			// no, position is beyond the current segment; subtract the length of current segment and step to next one
			position -= length;
		}
		// if we get here, position is to the right of the last segment; not good.
		throw new ReviewedStingException("Position "+position + " is outside of custom contig boundaries");
	}

	private long contigLength(Collection<GenomeLoc> segments) {
		long l = 0;
		for ( GenomeLoc g : segments ) l += (g.getStop() - g.getStart() + 1 );
		return l;
	}
	
	public static void main(String argv[]) {
		
//		SAMFileReader reader = new SAMFileReader(new java.io.File("/humgen/gsa-scr1/asivache/TCGA/Ovarian/C2K/0904/normal.bam"));
        SAMFileReader reader = new SAMFileReader(new java.io.File("X:/asivache/cDNA/new_pipeline/30BV1/test.1.sam"));

		
		SAMRecord r = new SAMRecord(reader.getFileHeader());
        GenomeLocParser genomeLocParser = new GenomeLocParser(reader.getFileHeader().getSequenceDictionary());

        r.setReferenceName("ENST00000378466");
        r.setAlignmentStart(1235);
        r.setCigarString("24M1D27M");

//		List<GenomeLoc> s = new ArrayList<GenomeLoc>();
//		s.add( GenomeLocParser.createGenomeLoc("chr1", 100, 199));	
//		s.add( GenomeLocParser.createGenomeLoc("chr1", 300, 499));	
//		s.add( GenomeLocParser.createGenomeLoc("chr1", 600, 799));	

		GenomicMap m = new GenomicMap(5);
		
//		m.readArachne(genomeLocParser,new File("/humgen/gsa-scr1/asivache/cDNA/Ensembl48.transcriptome.map"));
//        m.write(new File("/humgen/gsa-scr1/asivache/cDNA/new_pipeline/Ensembl48.new.transcriptome.map"));
        m.read(genomeLocParser,new File("W:/berger/cDNA_BAM/refs/Ensembl52.plus.Genome.map"));

        m.remapToMasterReference(r,reader.getFileHeader(),true);

//		if ( m.getContigMapping("ENST00000302418") == null ) System.out.println("ERROR! CONTIG IS MISSING!");

		int cnt = 0;
		
		System.out.println(m.size() + " contigs loaded");
        System.out.println("new alignment: "+r.format())  ;
/*
 		for ( String name : m.nameSet() ) {
 
			System.out.print(name);
			System.out.print(": ");
			for ( GenomeLoc g : m.getContigMapping(name)) {
				System.out.print(g.toString()+",  ");
			}
			System.out.println();
			cnt ++;
			if ( cnt > 10 ) break;
		}
*/		
//		m.addCustomContig("My", s);
/*		
		r.setReferenceName("My");
		r.setAlignmentStart(3);
		r.setCigarString("5S97M5D197M5H");
		
		m.remapToMasterReference(r);
		System.out.println(r.getReferenceName()+":"+r.getAlignmentStart()+" "+r.getCigarString());
*/		
		reader.close();
			
	}
	
}
