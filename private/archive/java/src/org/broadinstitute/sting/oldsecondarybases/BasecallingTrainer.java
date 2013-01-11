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

package org.broadinstitute.sting.secondarybase;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.reference.ReferenceSequence;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * BasecallingTrainingSet holds a set of raw read sequences, their raw intensities, and quality scores.
 *
 * @author Kiran Garimella
 */
public class BasecallingTrainer {
    private File bustardDir;
    private int lane;
    private int trainingLimit;

    private ArrayList<RawRead> trainingData;

    /**
     * Constructor for BasecallingTrainingSet.
     *
     * @param bustardDir     the Bustard directory for the sample
     * @param lane           the lane for the sample
     * @param trainingLimit  the number of training reads to accept
     */
    public BasecallingTrainer(File bustardDir, int lane, int trainingLimit) {
        this.bustardDir = bustardDir;
        this.lane = lane;
        this.trainingLimit = trainingLimit;
    }

    /**
     * Get the training data array list.
     *
     * @return  the arraylist of raw training reads
     */
    public ArrayList<RawRead> getTrainingData() {
        return this.trainingData;
    }

    /**
     * Set the training data array list.
     *
     * @param trainingData  the arraylist of raw training reads
     */
    public void setTrainingData(ArrayList<RawRead> trainingData) {
        this.trainingData = trainingData;
    }

    /**
     * Take the first N reads that have no ambiguous bases, an average quality score greater
     * than or equal to 15, and are not largely homopolymers and add them to the training set.
     */
    public void loadFirstNReasonableReadsTrainingSet() {
        this.trainingData = new ArrayList<RawRead>(trainingLimit);

        IlluminaParser iparser = new IlluminaParser(bustardDir, lane);

        RawRead rawread;
        int numreads = 0;

        while (numreads < trainingLimit && iparser.next()) {
            rawread = iparser.getRawRead();

            int numAmbiguous = 0;
            byte[] sequence = rawread.getSequence();
            
            for ( byte byteBase : sequence ) {
                if (BaseUtils.simpleBaseToBaseIndex((char) byteBase) == -1) {
                    numAmbiguous++;
                }
            }

            if (numAmbiguous == 0 && getAverageQualityScore(rawread) >= 15 && BaseUtils.mostFrequentBaseFraction(rawread.getSequence()) < 0.4) {
                trainingData.add(rawread);
                numreads++;
            }
        }
    }

    /**
     * Take the first N reads that have no ambiguous bases, an average quality score greater
     * than or equal to 15, and are not largely homopolymers and add them to the training set.
     *
     * @param bustardDir     the bustard directory
     * @param lane           the lane number
     * @param trainingLimit  how many reads should we use to train?
     */
    public static void loadNReasonableReadsTrainingSet(File bustardDir, int lane, int trainingLimit) {
        ArrayList<RawRead> trainingData = new ArrayList<RawRead>(trainingLimit);

        IlluminaParser iparser = new IlluminaParser(bustardDir, lane);

        RawRead rawread;
        int numreads = 0;

        while (numreads < trainingLimit && iparser.next()) {
            rawread = iparser.getRawRead();

            int numAmbiguous = 0;
            byte[] sequence = rawread.getSequence();

            for ( byte byteBase : sequence ) {
                if (BaseUtils.simpleBaseToBaseIndex((char) byteBase) == -1) {
                    numAmbiguous++;
                }
            }

            if (numAmbiguous == 0 && getAverageQualityScore(rawread) >= 15 && BaseUtils.mostFrequentBaseFraction(rawread.getSequence()) < 0.4) {
                trainingData.add(rawread);
                numreads++;
            }
        }
    }

    /**
     * Return the average quality score of a raw read.
     *
     * @param read  the raw read
     * @return  the average quality score
     */
    private static double getAverageQualityScore(RawRead read) {
        double averageQual = 0;

        for ( byte qual : read.getQuals() ) {
            averageQual += qual;
        }

        return averageQual / ((double) read.getReadLength());
    }

    /**
     * Load a training set from perfect reads in an already-aligned bam file.
     *
     * @param samIn      the SAM/BAM file to load the reads from
     * @param reference  the reference to which the reads should be compared
     */
    public void loadPreAlignedTrainingSet(File samIn, File reference) {
        Vector< HashMap<String, SAMRecord> > trainingReads = getPerfectAlignmentsByTile(samIn, reference);

        trainingData = correlateReadsAndIntensities(trainingReads);
    }

    /**
     * Find perfect reads and group them by tile.
     *
     * @param samIn      the SAM/BAM file to load the raeds from
     * @param reference  the reference to which the reads should be compared
     * @return a vector of perfect reads, grouped by tile
     */
    private Vector<HashMap<String, SAMRecord>> getPerfectAlignmentsByTile(File samIn, File reference) {
        ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(reference);
        String currentContig = "none";
        byte[] refbases = null;

        SAMFileReader sf = new SAMFileReader(samIn);
        SAMRecord sr;
        CloseableIterator<SAMRecord> sfit = sf.iterator();

        Vector< HashMap<String, SAMRecord> > trainingReads = new Vector< HashMap<String, SAMRecord> >(101);
        int numTrainingReads = 0;

        while (numTrainingReads < trainingLimit && (sr = sfit.next()) != null) {
            if (sr.getCigar().numCigarElements() == 1) {
                int offset = sr.getAlignmentStart();

                if (!currentContig.matches(sr.getReferenceName())) {
                    ReferenceSequence refSeq = ref.nextSequence();
                    while( !refSeq.getName().equals(sr.getReferenceName()) )
                        refSeq = ref.nextSequence();
                    currentContig = sr.getReferenceName();
                    refbases = refSeq.getBases();
                }

                int mismatches = 0;

                String refString = "";
                for (int i = offset, j = 0; i < offset + sr.getReadBases().length; i++, j++) {
                    refString += (char) refbases[i - 1];

                    mismatches += (BaseUtils.simpleBaseToBaseIndex((char) refbases[i - 1]) !=
                                   BaseUtils.simpleBaseToBaseIndex((char) sr.getReadBases()[j]))
                                   ? 1 : 0;
                }

                if (mismatches == 0) {
                    Pattern p = Pattern.compile(":(\\d):(\\d+):(\\d+):(\\d+)#");
                    Matcher m = p.matcher(sr.getReadName());

                    if (m.find()) {
                        int tile = Integer.valueOf(m.group(2));
                        String readKey = String.format("%s:%s:%s:%s", m.group(1), m.group(2), m.group(3), m.group(4));

                        if (tile > trainingReads.size()) {
                            trainingReads.setSize(tile + 1);
                        }

                        if (trainingReads.get(tile) == null) {
                            trainingReads.set(tile, new HashMap<String, SAMRecord>());
                        }

                        trainingReads.get(tile).put(readKey, sr);
                        numTrainingReads++;
                    } else {
                        throw new StingException("Pattern '" + p.pattern() + "' does not match against read name '" + sr.getReadName() + "'");
                    }
                }
            }
        }

        return trainingReads;
    }

    /**
     * Correlate the perfect reads with their intensities (at least, theoretically).  This doesn't work right now...
     *
     * @param trainingReads  the set of training reads, hashed by tile
     * @return the final training set with intensities added in
     */
    private ArrayList<RawRead> correlateReadsAndIntensities(Vector<HashMap<String, SAMRecord>> trainingReads) {
        ArrayList<RawRead> newTrainingData = new ArrayList<RawRead>(trainingLimit);

        IlluminaParser iparser = new IlluminaParser(bustardDir, lane);

        int totalReadCount = 0;

        for (int tile = 1; tile < trainingReads.size(); tile++) {
            iparser.seekToTile(tile);

            int tileReadCount = 0;

            RawRead iread;
            while (trainingReads.get(tile) != null && tileReadCount < trainingReads.get(tile).size() && iparser.next()) {
                iread = iparser.getRawRead();
                String readKey = iread.getReadKey();

                if (trainingReads.get(tile).containsKey(readKey)) {
                    System.out.printf("\tTile %d: found %d of %d (%4.4f in tile, %4.4f total)            \r",
                                      tile,
                                      tileReadCount,
                                      trainingReads.get(tile).size(),
                                      ((double) tileReadCount)/((double) trainingReads.get(tile).size()),
                                      ((double) totalReadCount)/((double) trainingLimit));

                    byte[] quals = new byte[iread.getReadLength()];
                    for (int qualIndex = 0; qualIndex < quals.length; qualIndex++) {
                        quals[qualIndex] = 40;
                    }

                    iread.setQuals(quals);
                    newTrainingData.add(iread);

                    tileReadCount++;
                    totalReadCount++;
                }
            }
        }

        iparser.close();

        return newTrainingData;
    }
}
