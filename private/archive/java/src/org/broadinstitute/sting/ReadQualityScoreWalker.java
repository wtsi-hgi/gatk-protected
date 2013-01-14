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

package org.broadinstitute.sting.walkers;

import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileWriter;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Oct 23, 2009
 *
 * This walker is designed to work as the second pass in a two-pass processing step.
 * It does a by-read traversal calculating a read quality score based on a number of factors:
 *  1.) Average neighborhood quality score for the length of the read (data generated by NeighborhoodQualityWalker)
 *  2.) Is this read's mate mapped to a different chromosome?
 *  3.) The mapping quality for this read.
 *  4.) Number of reference mismatches in this read.
 * This walker creates a new bam file in which each read is annotated by this read quality score
 *  in addition if the read quality score is below the given threshold, the read is flagged.
 *
 * This walker requires as input the file of (GenomeLoc QualityScore)'s generated by NeighborhoodQualityWalker.
 * This walker accepts as input a threshold in order to flag reads which are of unacceptable read quality.
 *
 * This walker is designed to be used in conjunction with NeighborhoodQualityWalker.
 */

public class
        ReadQualityScoreWalker extends ReadWalker<SAMRecord, SAMFileWriter> {
    @Output
    protected PrintStream out;
    @Argument(fullName = "inputQualityFile", shortName = "if", doc = "Input quality score file generated by NeighborhoodQualityWalker", required = true)
    protected String inputQualityFile = null;
    @Argument(fullName = "outputBamFile", shortName = "of", doc = "Write output to this BAM filename instead of STDOUT", required = false)
    protected SAMFileWriter outputBamFile = null;
    @Argument(fullName = "threshold", shortName = "th", doc="Flag reads whose read quality score is below this threshold", required = false)
    protected int qualityThreshold = 13;

    private BufferedReader inputReader = null;
    private static String line = null;

    public SAMRecord map( ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker ) {
        return read; // all the work is done in the reduce step for this walker
    }

    public SAMFileWriter reduceInit() {
        try {
            inputReader = new BufferedReader( new FileReader ( inputQualityFile ) );
        } catch ( FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(new File(inputQualityFile), e);
		} catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(new File(inputQualityFile), e);
		}
        return outputBamFile;
    }

    public SAMFileWriter reduce( SAMRecord read, SAMFileWriter output ) {

		int readQualityScore = 0;
		float meanNeighborhoodQuality = 0.0f;

        // The large block of code below is parsing through the input file and calculating the meanNeighborhoodQuality over the length of the read
        //  It does this by first skipping in the file to where the current read starts and marking that location
        //  Next it continues reading lines for the length of the read generating a sum of neighborhood quality
        //  When it reaches the end of the read it jumps back to the marker so that it can be used by the next read
        // BUGBUG: This assumes reads will be sorted by start location
        float sumNeighborhoodQuality = 0.0f;
        int numLines = 0;
        GenomeLoc readLoc = getToolkit().getGenomeLocParser().createGenomeLoc( read );
        if( readLoc.size() > 0 ) { // only calculate mean NQS if the read has a well formed GenomeLoc, if not NQS will be zero
            try {
                if( line == null ) {
                    line = inputReader.readLine();
                    if( line == null ) { throw new UserException.MalformedFile(new File(inputQualityFile), "Input file is empty" ); }
                }
                String[] halves = line.split( " ", 2 );
                GenomeLoc curLoc = getToolkit().getGenomeLocParser().parseGenomeLoc( halves[0] );
                while( curLoc.isBefore( readLoc ) ) { // Loop until the beginning of the read
                    line = inputReader.readLine();
                    if( line == null ) { throw new UserException.MalformedFile(new File(inputQualityFile), "Input file doesn't encompass all reads. Can't find beginning of read: " + readLoc ); }
                    halves = line.split( " ", 2 );
                    curLoc = getToolkit().getGenomeLocParser().parseGenomeLoc( halves[0] );
                }
                // now we have skipped ahead in the input file to where this read starts
                logger.debug( "Starting: " + curLoc + ", read: " + readLoc + "\t size: " + readLoc.size() );
                inputReader.mark( 30 * ( (int)readLoc.size() + 3 ) ); // BUGBUG: Is this a sufficient buffer size?
                String savedLine = line;

                while( !curLoc.isPast( readLoc ) ) {  // Loop until just past the end of the read
                    sumNeighborhoodQuality += Float.parseFloat( halves[1] );
                    numLines++;
                    line = inputReader.readLine();
                    if( line == null ) { throw new UserException.MalformedFile(new File(inputQualityFile), "Input file doesn't encompass all reads. Can't find end of read: " + readLoc ); }
                    halves = line.split( " ", 2 );
                    curLoc = getToolkit().getGenomeLocParser().parseGenomeLoc( halves[0] );
                }
                // now we have parsed the input file up to where the read ends
                // reset back to the mark in order to parse the next read in the next call to the reduce function
                inputReader.reset();
                line = savedLine;

            } catch ( FileNotFoundException e ) {
                throw new UserException.CouldNotReadInputFile(new File(inputQualityFile), e);
            } catch (IOException e ) {
                throw new UserException.CouldNotReadInputFile(new File(inputQualityFile), e);
            }

            meanNeighborhoodQuality = sumNeighborhoodQuality / ((float) numLines);
        }

		
        // Find out if this read's mate mapped to a different chromosome
        //boolean isGoodPair = ( read.getReadPairedFlag() ? read.getProperPairFlag() : true );
        boolean isGoodPair = ( !read.getReadPairedFlag() || read.getProperPairFlag() ); // optimized version of above line

        // Get the mapping quality for this read
        int mappingQuality = read.getMappingQuality();

        // Get the number of reference mismatches in this read
        assert read.getReadLength() > 0 : "Read length must be greater than zero.";
        float mismatchRate = 1.0f;
        if( read.getAttribute("NM") != null ) {
            mismatchRate = ((float) Integer.parseInt(read.getAttribute("NM").toString())) / ((float) read.getReadLength());
        }

		
		// Calculate the three additional metrics that go into a read quality score
		// BUGBUG: some analysis is needed to determine reasonable quality values and rates for the exponentials
		float scoreMate = ( isGoodPair ? 40.0f : 2.0f );
		float scoreMapping = 40.0f * (float) Math.exp( -0.02f * Math.max( 99.0f - mappingQuality, 0.0f ) );
									// exp decay with rate 0.02, scaled to Q=40 when mapping quality is 99
        float scoreMismatch = 40.0f * (float) Math.exp( -27.0f * mismatchRate );
									// exp decay with rate 27.0, scaled to Q=40 when the mismatch rate is 0% for this read
		
		// BUGBUG: some analysis is needed to determine reasonable weights for each metric
		readQualityScore = Math.round( 0.6f * meanNeighborhoodQuality + 0.1f * scoreMate + 0.05f * scoreMapping + 0.25f * scoreMismatch );
        if( readQualityScore == 0 ) { readQualityScore = 1; }
        assert readQualityScore > 0 : "Read quality score must be positive and nonzero.";

		// Add the read quality score to the read in the new bam file and flag it if quality is below the given threshold
        // BUGBUG: which attributes should be set here?
        read.setAttribute( "XR", readQualityScore );
        if( readQualityScore < qualityThreshold ) { 
            read.setAttribute( "ZR", 1 );
        }

        // verbose debug printing lines
        logger.debug( read.getReadName() + " " + readQualityScore );
        logger.debug( "neighborhood quality =\t" + meanNeighborhoodQuality );
        logger.debug( "mate mismatch? =\t" + isGoodPair + " --> " + scoreMate );
        logger.debug( "mapping quality =\t" + mappingQuality + " --> " + scoreMapping );
        logger.debug( "ref mismatch rate =\t" + mismatchRate + " --> " + scoreMismatch );

        // This printout useful for making histograms of scores in Matlab
        //out.println( readQualityScore + " " + meanNeighborhoodQuality + " " + scoreMate + " " + scoreMapping + " " + scoreMismatch );

        // Add the read to the output bam file or output to STDOUT
        if ( output != null ) {
            output.addAlignment( read );
        } else {
            out.println( read.format() );
        }

        return output;
    }

    public void onTraversalDone( SAMFileWriter reduceResult ) {
    }
    
}
