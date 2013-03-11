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

package org.broadinstitute.sting.gatk.walkers.readutils;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.picard.PicardException;
import net.sf.picard.fastq.*;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.sam.SamToFastq;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMValidationError;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Advanced;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
/**
 * Utility tool to blindly strip base adaptors. Main application is for FASTQ/unaligned BAM pre-processing where libraries
 * have very short inserts, and hence a substantial part of the sequencing data will have adaptor sequence present.
 * <p>
 * By design, tool will only work for Illumina-like library constructs, where the typical library architecture is:
 * [Adaptor 1]-[Genomic Insert]-[Adaptor 2 (index/barcode)]
 * <p>
 * It is assumed that when data is paired, one read will span the forward strand and one read will span the reverse strand.
 * Hence, when specifying adaptors they should be specified as both forward and reverse-complement to make sure they're removed in all cases.
 * By design, as well, "circular" constructions where a read can have an insert, then adaptor, then more genomic insert, are not supported.
 * When an adaptor is detected, all bases downstream from it (i.e. in the 3' direction) will be removed.
 * Adaptor detection is carried out by looking for overlaps between forward and reverse reads in a pair.
 * If a sufficiently high overlap is found, the insert size is computed and if insert size < read lengths adaptor bases are removed from reads.
 *
 * Advantages over ReadClipper:
 * - No previous knowledge of adaptors or library structure is necessary
 *
 * Advantages over 3rd party tools like SeqPrep:
 * - Can do BAM streaming instead of having to convert to fastq
 * - No need to merge reads - merging reads can have some advantages, but complicates downstream processing and loses information that can be used,
 *   e.g. in variant calling
 * <p>
 *
 * <h2>Input</h2>
 * <p>
 * The input read data in BAM format. Read data MUST be in query name ordering as produced, for example with Picard's FastqToBam
 *
 * <h2>Output</h2>
 * <p>
 * A merged BAM file with unaligned reads
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -T ReadAdaptorTrimmer \
 *   -I my_reads.bam \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -o trimmed_Reads.bam
 * </pre>
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_DATA, extraDocs = {CommandLineGATK.class} )
@PartitionBy(PartitionType.READ)
public class ReadAdaptorTrimmer extends ReadWalker<List<GATKSAMRecord>, SAMFileWriter> implements NanoSchedulable {
    @Output(doc="Write output to this BAM filename instead of STDOUT", required = false)
    SAMFileWriter out;

     /**
     * Only prints the first n reads of the file - for short testing
     */
     @Hidden
    @Argument(fullName = "number", shortName = "n", doc="Print the first n reads from the file, discarding the rest", required = false)
    int nReadsToPrint = -1;

    /**
     * Argument to control strictness of match between forward and reverse reads - by default, we require 15 matches between them to declare
     * an overlap.
     */
    @Advanced
    @Argument(fullName = "minMatches", shortName = "minMatches", doc="Minimum number of substring matches to detect pair overlaps", required = false)
    int minMatchesForOverlap = 15;


     /**
     * private class members
     */
    private GATKSAMRecord firstReadInPair;
    private MergeStats mergeStats = new MergeStats();

    static class MergeStats {
        long numReadsProcessed;
        long numReadsWithAdaptorTrimmed;
    }

   /**
     * The reads filter function.
     *
     * @param ref  the reference bases that correspond to our read, if a reference was provided
     * @param read the read itself, as a GATKSAMRecord
     * @return true if the read passes the filter, false if it doesn't
     */
    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
         // check if we've reached the output limit
        if ( nReadsToPrint == 0 ) {
            return false;          // n == 0 means we've printed all we needed.
        }
        else if (nReadsToPrint > 0) {
            nReadsToPrint--;       // n > 0 means there are still reads to be printed.
        }
        return true;
    }
    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     *
     * @return SAMFileWriter, set to the BAM output file if the command line option was set, null otherwise
     */
    public SAMFileWriter reduceInit() {
        return out;
    }

    public List<GATKSAMRecord> map( final ReferenceContext ref, final GATKSAMRecord readIn, final RefMetaDataTracker metaDataTracker ) {


        final List<GATKSAMRecord> readsToEmit = new ArrayList<GATKSAMRecord>();


        // cache first read in pair if flag set.
        if (readIn.getFirstOfPairFlag()) {
            firstReadInPair = GATKSAMRecord.emptyRead(readIn);
            firstReadInPair.setReadString(readIn.getReadString());
            firstReadInPair.setReadName(readIn.getReadName());
            firstReadInPair.setBaseQualities(readIn.getBaseQualities());
        }
        else {
            if (!readIn.getReadName().matches(firstReadInPair.getReadName()))
                throw new IllegalStateException("Second read in pair must follow first read in pair: data not ordered?");

            final int oldLength1 = firstReadInPair.getReadLength();
            final int oldLength2 = readIn.getReadLength();
            // try to strip any adaptor sequence in read pair
            final Integer result = mergeReads(firstReadInPair, readIn, minMatchesForOverlap, logger);

            logger.debug("mergeReads result = " + result);

            readsToEmit.add(firstReadInPair);
            readsToEmit.add(readIn);

            if (oldLength1 != firstReadInPair.getReadLength())
                mergeStats.numReadsWithAdaptorTrimmed++;
            if (oldLength2 != readIn.getReadLength())
                mergeStats.numReadsWithAdaptorTrimmed++;

         }


        mergeStats.numReadsProcessed++;
        return readsToEmit;

    }

    /**
     * given a read and a output location, reduce by emitting the read
     *
     * @param readsToEmit   the read itself
     * @param output the output source
     * @return the SAMFileWriter, so that the next reduce can emit to the same source
     */
    public SAMFileWriter reduce( final List<GATKSAMRecord> readsToEmit, final SAMFileWriter output ) {
        for (final GATKSAMRecord read : readsToEmit)
             output.addAlignment(read);

        return output;
    }

    @Override
    public void onTraversalDone(SAMFileWriter output) {

        logger.info("Finished merging:");
        logger.info("Number of processed reads:                     "+mergeStats.numReadsProcessed);
        logger.info("Number of reads with adaptor sequence trimmed: "+mergeStats.numReadsWithAdaptorTrimmed);
    }


    /**
     *
     * Workhorse routines...
     *
     */
        /**
         * Core routine that does most underlying work for walker. Takes two reads and looks for overlaps in them.
         * An overlap is defined as a contiguous chunk of N bases that matches reverse-complement between reads.
         * Currently, the only insert structure that it will look for overlaps is as follows:
         * CASE 1: Insert shorter than read length:
         * 3' XXXXXXXXXXXXXXXX 5'            (second read)
         * 5'      YYYYYYYYYYYYYYYY 3'       (first read)
         *         ***********
         *
         * In this case, if X and Y are complements at the 11 positions marked by *, routine will do the following
         * iff minMatchesForOverlap <= 11:
         *  a) Cleave adaptor from end of second read (leftmost dangling part in diagram above)
         *  b) Cleave adaptor from end of first read (rightmost part in diagram).
         *
         * CASE 2: Insert size >= read length:
         * 3'             XXXXXXXXXXXXXXXX 5'           (second read)
         * 5'      YYYYYYYYYYYYYYYY 3'                  (first read)
         *                *********                        (overlap)
         *
         * In this case, no trimming is done and reads are left unchanged
         * @param first                      (I/O) First read in pair - read contents (bases/quals) can be modified if adaptor is detected
         * @param second                     (I/O) Second read in pair - read contents (bases/quals) can be modified if adaptor is detected
         * @param minMatchesForOverlap       Reads need to match in these # of bases to be joined
         * @return                           Offset between second and first read.
         *                                   If there's no detectable offset, return Null
         */
    @Requires({"first != null","second != null","minMatchesForOverlap>0"})
    protected static Integer mergeReads(final GATKSAMRecord first,
                                              final GATKSAMRecord second,
                                              final int minMatchesForOverlap,
                                              final Logger logger) {

        final Integer insertSize = estimateInsertSize(first.getReadBases(), second.getReadBases(),
                minMatchesForOverlap, logger);

        if (insertSize == null)
            return insertSize;
        if (insertSize < first.getReadLength()) {
            // trim adaptor sequence from read
            first.setReadBases(Arrays.copyOfRange(first.getReadBases(),0,insertSize));
            first.setBaseQualities(Arrays.copyOfRange(first.getBaseQualities(),0,insertSize));
        }
        if (insertSize < second.getReadLength()) {
            // trim adaptor sequence from read
            second.setReadBases(Arrays.copyOfRange(second.getReadBases(),0,insertSize));
            second.setBaseQualities(Arrays.copyOfRange(second.getBaseQualities(),0,insertSize));
        }
        return insertSize;
    }

    /**
    * Brain-dead implementation of an aligner of two sequences, where it's assumed that there might be an overlap
    * from the first into the second. From this, an estimate of insert size is performed and returned
    * Assumes that reads come in reverse direction, so one of the base sequences needs to be reverse-complemented.]
    *
    * @param firstRead                           Bytes from first read
    * @param secondRead                          Bytes from second read (reverse direction)
    * @return                                  Estimated insert size based on offset between first and second read.
    *                                          If no overlap can be detected, return null
    */

    @Requires({"firstRead != null","secondRead != null","minMatches>0","firstRead.length == secondRead.length"})
    protected static Integer estimateInsertSize(final byte[] firstRead,
                                                                final byte[] secondRead,
                                                                final int minMatches,
                                                                final Logger logger) {
        final byte[] firstBases = firstRead;
        final byte[] secondBases = BaseUtils.simpleReverseComplement(secondRead);

        final Pair<Integer,Integer> overlaps = findOverlappingSequence(firstBases, secondBases);
        final int bestOffset = overlaps.first;
        final int maxScore = overlaps.second;
        if ( logger.isDebugEnabled()) {
            String sb="", s1 = new String(firstBases), s2 = new String(secondBases);
            for (int k=0; k < Math.abs(bestOffset); k++) sb+=" ";
            if (maxScore >= minMatches) {
                logger.debug(String.format("Match, Max Score = %d, best offset = %d\n",maxScore, bestOffset));
                if (bestOffset>0)
                    s2 = sb+s2;
                else
                    s1 = sb+s1;
            }
            else logger.debug("NoMatch:");
            logger.debug("R1:"+s1);
            logger.debug("R2:"+s2);


        }

        if (maxScore < minMatches)
            return null; // no overlap detected

        return bestOffset+secondRead.length;


    }


     /**
     * Tries to find overlapping sequence between two reads, and computes offset between them
      * For each possible offset, computes matching score, which is = MATCH_SCORE*Num_matches + MISMATCH_SCORE*num_mismatches
      * (like SW with infinite gap penalties).
     * @param first                              First read bytes
     * @param second                             Second read bytes
     * @return                                   Pair of integers (x,y). x = best offset between reads, y = corresponding score
     */
    @Requires({"first != null","second != null"})
    @Ensures("result != null")
    protected static Pair<Integer,Integer> findOverlappingSequence(final byte[] first,
                                                 final byte[] second) {
        final int MATCH_SCORE = 1;
        final int MISMATCH_SCORE = -1;
        // try every possible offset - O(N^2) algorithm

        // In case of following structure,
        //      111111111
        // 222222222
        // computed offset will be negative (=-5 in this case).
        // If however,
        //   111111111
        //      222222222
        // then offset will be positive (=3 in this case)
        int maxScore = 0, bestOffset =0;
        for (int offset = -second.length; offset < first.length; offset++) {
            int score = 0;
            // compute start index for each array
            int ind1 = (offset<0)?0:offset;
            int ind2 = (offset<0)?-offset:0;
            for (int k=0; k < Math.min(first.length, second.length) ; k++) {
                if (ind1 >= first.length)
                    break;
                if (ind2 >= second.length )
                    break;
                if (first[ind1] != 'N' && second[ind2] != 'N')  {
                    if (first[ind1] == second[ind2])
                        score += MATCH_SCORE;
                    else
                        score += MISMATCH_SCORE;
                }
                ind1++;
                ind2++;
            }
            if (score > maxScore) {
                maxScore = score;
                bestOffset = offset;
            }
        }
        return new Pair<Integer, Integer>(bestOffset,maxScore);
    }

}
