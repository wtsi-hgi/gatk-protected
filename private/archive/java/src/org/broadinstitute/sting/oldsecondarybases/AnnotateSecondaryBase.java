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

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.variant.utils.BaseUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * AnnotateSecondaryBase computes the second best base for every base in an Illumina lane.
 * First, a statistical model is fit to a subset of the raw Illumina intensities (i.e. those
 * generated by Illumina's "Firecrest" package).  Then, every read's set of raw intensities
 * is evaluated against this model to determine the base probability distribution of a given
 * base observation.
 *
 * Approximately 95% of the time, this method and Illumina's basecalling package, "Bustard",
 * agree on the identity of the best base.  In these cases, we simply annotate our estimate
 * of the second-best base.  In cases where this method and Bustard disagree, we annotate
 * the secondary base as this method's primary base.
 *
 * @author Kiran Garimella
 */

/*
 An example invocation:
   java -Xmx2048m -Djava.io.tmpdir=/broad/hptmp/ -jar /path/to/AnnotateSecondaryBase.jar \
    -D /seq/solexaproc2/SL-XAX/analyzed/090217_SL-XAX_0003_FC30R47AAXX/Data/C1-152_Firecrest1.3.2_25-02-2009_prodinfo/Bustard1.3.2_25-02-2009_prodinfo/ \
    -L 5 \
    -B 30R47AAXX090217
    -CR '0-75,76-151' \
    -SO ~/test.sam \
    -SI /seq/picard/30R47AAXX/C1-152_2009-02-17_2009-04-02/5/Solexa-10265/30R47AAXX.5.aligned.bam \
 */
    
public class AnnotateSecondaryBase extends CommandLineProgram {
    public static AnnotateSecondaryBase Instance = null;

    @Argument(fullName="bustard_dir", shortName="D", doc="Illumina Bustard directory") public File BUSTARD_DIR;
    @Argument(fullName="lane", shortName="L", doc="Illumina flowcell lane") public int LANE;
    @Argument(fullName="run_barcode", shortName="B", doc="Run barcode (embedded as part of the read name; i.e. 30R47AAXX090217)") public String RUN_BARCODE;
    @Argument(fullName="cycle_ranges", shortName="CR", doc="Cycle ranges for single-end or paired reads (i.e. '0-50,56-106') (0-based, inclusive)") public String CYCLE_RANGES;
    @Argument(fullName="sam_out", shortName="SO", doc="Output path for sam file") public File SAM_OUT;
    @Argument(fullName="sam_in", shortName="SI", doc="The file to use for training and annotation", required=false) public File SAM_IN;
    @Argument(fullName="training_limit", shortName="T", doc="Number of reads to use for parameter initialization", required=false) public int TRAINING_LIMIT = 100000;
    @Argument(fullName="calling_limit", shortName="C", doc="Number of reads to basecall", required=false) public int CALLING_LIMIT = Integer.MAX_VALUE;
    @Argument(fullName="unaligned_sam", shortName="US", doc="Unaligned sam file, so we can skip making it", required=false) public File USAM;
    @Argument(fullName="aligned_sam", shortName="AS", doc="Aligned, queryname-sorted sam file, so we can skip resorting it", required=false) public File ASAM;

    public static void main(String[] argv) {
        Instance = new AnnotateSecondaryBase();
        start(Instance, argv);
    }

    protected int execute() {
        ArrayList<Pair<Integer, Integer>> cycleRanges = getCycleRanges(CYCLE_RANGES);
        File unalignedSam;

        if (USAM == null || !USAM.exists()) {
            BasecallingTrainer trainer = new BasecallingTrainer(BUSTARD_DIR, LANE, TRAINING_LIMIT);

            // Iterate through raw Firecrest data and store the first N reasonable reads up to TRAINING_LIMIT
            System.out.println("Loading training set from the first " + TRAINING_LIMIT + " reasonable reads in the raw data...");
            trainer.loadFirstNReasonableReadsTrainingSet();

            // Iterate through the stored training data and add the info to the BasecallingReadModel
            System.out.println("Applying training set...");
            BasecallingReadModel model = new BasecallingReadModel(trainer.getTrainingData());

            // Call bases and write results
            System.out.println("Calling bases...");

            SAMFileHeader sfh = new SAMFileHeader();
            sfh.setSortOrder(SAMFileHeader.SortOrder.queryname);

            unalignedSam = (canAnnotate(SAM_IN)) ? getTempSAMFile("unaligned") : SAM_OUT;
            SAMFileWriter sfw = new SAMFileWriterFactory().makeSAMOrBAMWriter(sfh, false, unalignedSam);

            IlluminaParser iparser = new IlluminaParser(BUSTARD_DIR, LANE);
    
            BasecallingStats bstats = new BasecallingStats();

            while (bstats.getReadsTotal() < CALLING_LIMIT && iparser.next()) {
                RawRead rr = iparser.getRawRead();
                FourProbRead fpr = model.call(rr);

                for (int cycleRangeIndex = 0; cycleRangeIndex < cycleRanges.size(); cycleRangeIndex++) {
                    Pair<Integer, Integer> cycleRange = cycleRanges.get(cycleRangeIndex);

                    RawRead rrEnd = iparser.getSubset(cycleRange.getFirst(), cycleRange.getSecond());
                    FourProbRead fprEnd = fpr.getSubset(cycleRange.getFirst(), cycleRange.getSecond());

                    sfw.addAlignment(constructSAMRecord(rrEnd, fprEnd, sfh, RUN_BARCODE, cycleRanges.size() == 2, cycleRangeIndex == 1));

                    if (cycleRangeIndex == 0) {
                        bstats.update(rrEnd, fprEnd);
                        bstats.notifyOnInterval(5000);
                    }
                }
            }

            bstats.notifyNow();

            iparser.close();
            sfw.close();
        } else {
            unalignedSam = USAM;
        }

        if (canAnnotate(SAM_IN)) {
            // If we're in annotation mode, annotate the aligned BAM file with the SQ tag
            System.out.println("Annotating aligned SAM file...");

            File alignedSam;
            if (ASAM == null || !ASAM.exists()) {
                System.out.println("  sorting aligned SAM file by read name...");
                alignedSam = getTempSAMFile("aligned");
                sortBAMByReadName(SAM_IN, alignedSam);
            } else {
                alignedSam = ASAM;
            }

            System.out.println("  merging unaligned and aligned SAM files...");
            File mergedSam = SAM_OUT;
            mergeUnalignedAndAlignedBams(unalignedSam, alignedSam, mergedSam);
        }

        System.out.println("Done.");

        return 0;
    }

    /**
     * Return a tempfile.  This is a laziness method so that I don't have to litter my code with try/catch blocks for IOExceptions.
     *
     * @param prefix  the prefix for the temp file
     * @return  the temp file
     */
    private File getTempSAMFile(String prefix) {
        try {
            File tempFile = File.createTempFile(prefix, ".sam", SAM_OUT.getParentFile());
            //tempFile.deleteOnExit();

            // Ensure that the volumes we're about to use are ready.
            PathUtils.refreshVolume(tempFile);
            PathUtils.refreshVolume(new File(System.getProperty("java.io.tmpdir")));

            return tempFile;
        } catch (IOException e) {
            throw new StingException("Unable to create tempfile in directory " + SAM_OUT.getParent());
        }
    }

    /**
     * Parse the cycle_ranges string that defines the cycle where a read starts and stops.
     * Comma-separated ranges are interpreted to be the first and second end of a pair.
     *
     * @param cycleRangesString  the 0-based, inclusive, comma-separated ranges (i.e. '0-50,51-100')
     * @return an ArrayList of cycle ranges
     */
    private ArrayList< Pair<Integer, Integer> > getCycleRanges(String cycleRangesString) {
        ArrayList< Pair<Integer, Integer> > cycleRanges = new ArrayList< Pair<Integer, Integer> >();

        String[] pieces = cycleRangesString.split(",");
        
        Pattern p = Pattern.compile("(\\d+)-(\\d+)");

        for (String piece : pieces) {
            Matcher m = p.matcher(piece);

            if (m.find()) {
                Integer cycleStart = new Integer(m.group(1));
                Integer cycleStop = new Integer(m.group(2));

                cycleRanges.add(new Pair<Integer, Integer>(cycleStart, cycleStop));
            }
        }

        if (cycleRanges.size() == 0) {
            throw new StingException("At least one cycle range must be specified.");
        }

        if (cycleRanges.size() > 2) {
            throw new StingException(cycleRanges.size() + " specified, but we're unable to handle more than 2.");
        }

        return cycleRanges;
    }

    /**
     * Simple test to determine whether we're in aligned bam annotation mode or not.
     *
     * @param samfile  the aligned sam file
     * @return true if the file exists, false otherwise
     */
    private boolean canAnnotate(File samfile) {
        return (samfile != null && samfile.exists());
    }

    /**
     * Construct a SAMRecord object with the specified information.  The secondary bases
     * will be annotated suchthat they will not conflict with the primary base.
     *
     * @param rr                 the raw Illumina read
     * @param fpr                the four-base distributions for every base in the read
     * @param sfh                the SAM header
     * @param runBarcode         the run barcode of the lane (used to prefix the reads)
     * @param isPaired           is this a paired-end lane?
     * @param isSecondEndOfPair  is this the second end of the pair?
     *
     * @return a fully-constructed SAM record
     */
    private SAMRecord constructSAMRecord(RawRead rr, FourProbRead fpr, SAMFileHeader sfh, String runBarcode, boolean isPaired, boolean isSecondEndOfPair) {
        SAMRecord sr = new SAMRecord(sfh);

        sr.setReadName(runBarcode + ":" + rr.getReadKey() + "#0");
        sr.setReadUmappedFlag(true);
        sr.setReadString(rr.getSequenceAsString());
        sr.setBaseQualities(rr.getQuals());

        sr.setReadPairedFlag(isPaired);
        if (isPaired) {
            sr.setMateUnmappedFlag(true);
            sr.setFirstOfPairFlag(!isSecondEndOfPair);
        }
        
        sr.setAttribute("SQ", fpr.getSQTag(rr));

        return sr;
    }

    /**
     * Resorts a SAM file to queryname order.
     *
     * @param samFile        the input SAM file
     * @param sortedSamFile  the sorted SAM output file
     */
    private void sortBAMByReadName(File samFile, File sortedSamFile) {
        SAMFileReader samIn = new SAMFileReader(samFile);

        SAMFileHeader sfh = samIn.getFileHeader();
        sfh.setSortOrder(SAMFileHeader.SortOrder.queryname);

        SAMFileWriter samOut = new SAMFileWriterFactory().makeSAMOrBAMWriter(sfh, false, sortedSamFile);

        for (SAMRecord sr : samIn) {
            samOut.addAlignment(sr);
        }

        samIn.close();
        samOut.close();
    }

    /**
     * Merges two SAM files that have been sorted in queryname order
     * 
     * @param queryNameSortedUnalignedSam  the sorted unaligned SAM file
     * @param queryNameSortedAlignedSam    the sorted aligned SAM file
     * @param mergedSam the output file where the merged results should be stored
     */
    private void mergeUnalignedAndAlignedBams(File queryNameSortedUnalignedSam, File queryNameSortedAlignedSam, File mergedSam) {
        SAMFileReader usam = new SAMFileReader(queryNameSortedUnalignedSam);
        SAMFileReader asam = new SAMFileReader(queryNameSortedAlignedSam);

        SAMFileHeader sfh = asam.getFileHeader();
        sfh.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        SAMFileWriter samOut = new SAMFileWriterFactory().makeSAMOrBAMWriter(sfh, false, mergedSam);

        CloseableIterator<SAMRecord> usamIt = usam.iterator();
        CloseableIterator<SAMRecord> asamIt = asam.iterator();

        SAMRecord usr = usamIt.next();
        SAMRecord asr = asamIt.next();

        int annotatedRecords = 0;

        do {
            // Pull a record from the unaligned file and advance the aligned file until we find the matching record.  We
            // don't have to advance the unaligned file until we find our record because we assume every record we generate
            // will be in the aligned file (which also contains unaligned reads).
            //
            // If Picard ever stops storing the unaligned reads, this logic will need to be rewritten.
            System.out.println(asr.getReadString());
            System.out.println(BaseUtils.simpleReverseComplement(asr.getReadString()));
            System.out.println();
            
            if (usr.getReadName().equals(asr.getReadName()) && usr.getFirstOfPairFlag() == asr.getFirstOfPairFlag()) {
                byte[] sqtag = (byte[]) usr.getAttribute("SQ");
                String usrread = usr.getReadString();
                String asrread = asr.getReadString();

                System.out.println(asrread);

                if (asr.getReadNegativeStrandFlag()) {
                    sqtag = QualityUtils.reverseComplementCompressedQualityArray(sqtag);
                    asrread = BaseUtils.simpleReverseComplement(asrread);

                    System.out.println(asrread);
                }

                if (usrread != null && asrread != null && !usrread.equals(asrread)) {
                    throw new StingException(
                        String.format("Purportedly identical unaligned and aligned reads have different read sequences.  Perhaps this lane was reanalyzed by the Illumina software but not the production pipeline?\n '%s:%b:%s'\n '%s:%b:%s'",
                                      usr.getReadName(), usr.getFirstOfPairFlag(), usrread,
                                      asr.getReadName(), asr.getFirstOfPairFlag(), asrread));
                }

                asr.setAttribute("SQ", sqtag);
                annotatedRecords++;

                System.out.println("Annotated " + annotatedRecords + " records.");

                usr = usamIt.next();
            } else {
                asr = asamIt.next();
            }

            samOut.addAlignment(asr);
        } while (usamIt.hasNext() && asamIt.hasNext());

        usam.close();
        asam.close();
        samOut.close();
    }

    /**
     * For debugging purposes.  Spits out relevant information for two SAMRecords.
     * 
     * @param sra  first SAMRecord
     * @param srb  second SAMRecord
     */
    private void printRecords(SAMRecord sra, SAMRecord srb) {
        System.out.println("a: " + sra.getReadName() + " " + sra.getFirstOfPairFlag());
        System.out.println("b: " + srb.getReadName() + " " + srb.getFirstOfPairFlag());
        System.out.println();

    }
}
