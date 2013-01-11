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

package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.*;
import net.sf.samtools.util.SequenceUtil;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.recalibration.BaseRecalibration;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.MissingResourceException;
import java.util.ResourceBundle;
import java.util.regex.Pattern;

/**
 * Second pass of the base quality score recalibration -- Uses the table generated by CountCovariates to update the base quality scores of the input bam file using a sequential table calculation making the base quality scores more accurately reflect the actual quality of the bases as measured by reference mismatch rate.
 *
 * <p>
 * This walker is designed to work as the second pass in a two-pass processing step, doing a by-read traversal. For each
 * base in each read this walker calculates various user-specified covariates (such as read group, reported quality score,
 * cycle, and dinuc). Using these values as a key in a large hashmap the walker calculates an empirical base quality score
 * and overwrites the quality score currently in the read. This walker then outputs a new bam file with these updated (recalibrated) reads.
 *
 * <p>
 * See the GATK wiki for a tutorial and example recalibration accuracy plots.
 * http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
 *
 * <h2>Input</h2>
 * <p>
 * The input read data whose base quality scores need to be recalibrated.
 * <p>
 * The recalibration table file in CSV format that was generated by the CountCovariates walker.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A bam file in which the quality scores in each read have been recalibrated.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -I my_reads.bam \
 *   -T TableRecalibration \
 *   -o my_reads.recal.bam \
 *   -recalFile my_reads.recal_data.csv
 * </pre>
 */

@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = BAQ.ApplicationTime.ON_OUTPUT)
@WalkerName("TableRecalibration")
@Requires({DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES})
// This walker requires -I input.bam, it also requires -R reference.fasta
public class TableRecalibrationWalker extends ReadWalker<SAMRecord, SAMFileWriter> {

    public static final String PROGRAM_RECORD_NAME = "GATK TableRecalibration";

    /////////////////////////////
    // Shared Arguments
    /////////////////////////////
    @ArgumentCollection
    private RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    /**
     * After the header, data records occur one per line until the end of the file. The first several items on a line are the
     * values of the individual covariates and will change depending on which covariates were specified at runtime. The last
     * three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches,
     * and the raw empirical quality score calculated by phred-scaling the mismatch rate.
     */
    @Input(fullName = "recal_file", shortName = "recalFile", required = true, doc = "Filename for the input covariates table recalibration .csv file")
    public File RECAL_FILE = null;
    /**
     * A new bam file in which the quality scores in each read have been recalibrated. The alignment of the reads is left untouched.
     */
    @Output(doc = "The output recalibrated BAM file", required = true)
    private StingSAMFileWriter OUTPUT_BAM = null;

    /**
     * TableRacalibration accepts a --preserve_qscores_less_than / -pQ <Q> flag that instructs TableRecalibration to not modify
     * quality scores less than <Q> but rather just write them out unmodified in the recalibrated BAM file. This is useful
     * because Solexa writes Q2 and Q3 bases when the machine has really gone wrong. This would be fine in and of itself,
     * but when you select a subset of these reads based on their ability to align to the reference and their dinucleotide effect,
     * your Q2 and Q3 bins can be elevated to Q8 or Q10, leading to issues downstream. With the default value of 5, all Q0-Q4 bases
     * are unmodified during recalibration, so they don't get inappropriately evaluated.
     */
    @Argument(fullName = "preserve_qscores_less_than", shortName = "pQ", doc = "Bases with quality scores less than this threshold won't be recalibrated. In general it's unsafe to change qualities scores below < 5, since base callers use these values to indicate random or bad bases", required = false)
    private int PRESERVE_QSCORES_LESS_THAN = 5;

    /**
     * By default TableRecalibration applies a Yates' correction to account for overfitting when it calculates the empirical
     * quality score, in particular, ( # mismatches + 1 ) / ( # observations + 1 ). TableRecalibration accepts a --smoothing / -sm <int>
     * argument which sets how many unobserved counts to add to every bin. Use --smoothing 0 to turn off all smoothing or, for example,
     * --smoothing 15 for a large amount of smoothing.
     */
    @Argument(fullName = "smoothing", shortName = "sm", required = false, doc = "Number of imaginary counts to add to each bin in order to smooth out bins with few data points")
    private int SMOOTHING = 1;

    /**
     * Combinations of covariates in which there are zero mismatches technically have infinite quality. We get around this situation
     * by capping at the specified value. We've found that Q40 is too low when using a more completely database of known variation like dbSNP build 132 or later.
     */
    @Argument(fullName = "max_quality_score", shortName = "maxQ", required = false, doc = "The integer value at which to cap the quality scores")
    private int MAX_QUALITY_SCORE = 50;

    /**
     * By default TableRecalibration emits the OQ field -- so you can go back and look at the original quality scores, rerun
     * the system using the OQ flags, etc, on the output BAM files; to turn off emission of the OQ field use this flag.
     */
    @Argument(fullName = "doNotWriteOriginalQuals", shortName = "noOQs", required = false, doc = "If true, we will not write the original quality (OQ) tag for each read")
    private boolean DO_NOT_WRITE_OQ = false;

    /////////////////////////////
    // Debugging-only Arguments
    /////////////////////////////
    @Hidden
    @Argument(fullName = "no_pg_tag", shortName = "noPG", required = false, doc = "Don't output the usual PG tag in the recalibrated bam file header. FOR DEBUGGING PURPOSES ONLY. This option is required in order to pass integration tests.")
    private boolean NO_PG_TAG = false;
    @Hidden
    @Argument(fullName = "fail_with_no_eof_marker", shortName = "requireEOF", required = false, doc = "If no EOF marker is present in the covariates file, exit the program with an exception.")
    private boolean REQUIRE_EOF = false;
    @Hidden
    @Argument(fullName = "skipUQUpdate", shortName = "skipUQUpdate", required = false, doc = "If true, we will skip the UQ updating step for each read, speeding up the calculations")
    private boolean skipUQUpdate = false;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private RecalDataManager dataManager; // Holds the data HashMap, mostly used by TableRecalibrationWalker to create collapsed data hashmaps
    private final ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>(); // List of covariates to be used in this calculation
    public static final Pattern COMMENT_PATTERN = Pattern.compile("^#.*");
    public static final Pattern OLD_RECALIBRATOR_HEADER = Pattern.compile("^rg,.*");
    public static final Pattern COVARIATE_PATTERN = Pattern.compile("^ReadGroup,QualityScore,.*");
    public static final String EOF_MARKER = "EOF";
    private long numReadsWithMalformedColorSpace = 0;

    /////////////////////////////
    //  Optimization
    /////////////////////////////
    private NestedHashMap qualityScoreByFullCovariateKey = new NestedHashMap(); // Caches the result of performSequentialQualityCalculation(..) for all sets of covariate values.

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Read in the recalibration table input file.
     * Parse the list of covariate classes used during CovariateCounterWalker.
     * Parse the CSV data and populate the hashmap.
     */
    public void initialize() {

        if (RAC.FORCE_PLATFORM != null) {
            RAC.DEFAULT_PLATFORM = RAC.FORCE_PLATFORM;
        }

        // Get a list of all available covariates
        final List<Class<? extends Covariate>> classes = new PluginManager<Covariate>(Covariate.class).getPlugins();

        int lineNumber = 0;
        boolean foundAllCovariates = false;

        // Read in the data from the csv file and populate the data map and covariates list
        logger.info("Reading in the data from input csv file...");

        boolean sawEOF = false;
        try {
            for (String line : new XReadLines(RECAL_FILE)) {
                lineNumber++;
                if (EOF_MARKER.equals(line)) {
                    sawEOF = true;
                }
                else if (COMMENT_PATTERN.matcher(line).matches() || OLD_RECALIBRATOR_HEADER.matcher(line).matches()) {
                    ; // Skip over the comment lines, (which start with '#')
                }
                // Read in the covariates that were used from the input file
                else if (COVARIATE_PATTERN.matcher(line).matches()) { // The line string is either specifying a covariate or is giving csv data
                    if (foundAllCovariates) {
                        throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration file. Found covariate names intermingled with data in file: " + RECAL_FILE);
                    }
                    else { // Found the covariate list in input file, loop through all of them and instantiate them
                        String[] vals = line.split(",");
                        for (int iii = 0; iii < vals.length - 3; iii++) { // There are n-3 covariates. The last three items are nObservations, nMismatch, and Qempirical
                            boolean foundClass = false;
                            for (Class<?> covClass : classes) {
                                if ((vals[iii] + "Covariate").equalsIgnoreCase(covClass.getSimpleName())) {
                                    foundClass = true;
                                    try {
                                        Covariate covariate = (Covariate) covClass.newInstance();
                                        requestedCovariates.add(covariate);
                                    } catch (Exception e) {
                                        throw new DynamicClassResolutionException(covClass, e);
                                    }

                                }
                            }

                            if (!foundClass) {
                                throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration file. The requested covariate type (" + (vals[iii] + "Covariate") + ") isn't a valid covariate option.");
                            }
                        }
                    }

                }
                else { // Found a line of data
                    if (!foundAllCovariates) {
                        foundAllCovariates = true;

                        // At this point all the covariates should have been found and initialized
                        if (requestedCovariates.size() < 2) {
                            throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration csv file. Covariate names can't be found in file: " + RECAL_FILE);
                        }

                        final boolean createCollapsedTables = true;

                        // Initialize any covariate member variables using the shared argument collection
                        for (Covariate cov : requestedCovariates) {
                            cov.initialize(RAC);
                        }
                        // Initialize the data hashMaps
                        dataManager = new RecalDataManager(createCollapsedTables, requestedCovariates.size());

                    }
                    addCSVData(RECAL_FILE, line); // Parse the line and add the data to the HashMap
                }
            }

        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(RECAL_FILE, "Can not find input file", e);
        } catch (NumberFormatException e) {
            throw new UserException.MalformedFile(RECAL_FILE, "Error parsing recalibration data at line " + lineNumber + ". Perhaps your table was generated by an older version of CovariateCounterWalker.");
        }
        logger.info("...done!");

        if (!sawEOF) {
            final String errorMessage = "No EOF marker was present in the recal covariates table; this could mean that the file is corrupted or was generated with an old version of the CountCovariates tool.";
            if (REQUIRE_EOF)
                throw new UserException.MalformedFile(RECAL_FILE, errorMessage);
            logger.warn(errorMessage);
        }

        logger.info("The covariates being used here: ");
        for (Covariate cov : requestedCovariates) {
            logger.info("\t" + cov.getClass().getSimpleName());
        }

        if (dataManager == null) {
            throw new UserException.MalformedFile(RECAL_FILE, "Can't initialize the data manager. Perhaps the recal csv file contains no data?");
        }

        // Create the tables of empirical quality scores that will be used in the sequential calculation
        logger.info("Generating tables of empirical qualities for use in sequential calculation...");
        dataManager.generateEmpiricalQualities(SMOOTHING, MAX_QUALITY_SCORE);
        logger.info("...done!");

        // Take the header of the input SAM file and tweak it by adding in a new programRecord with the version number and list of covariates that were used
        final SAMFileHeader header = getToolkit().getSAMFileHeader().clone();
        if (!NO_PG_TAG) {
            final SAMProgramRecord programRecord = new SAMProgramRecord(PROGRAM_RECORD_NAME);
            final ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
            try {
                final String version = headerInfo.getString("org.broadinstitute.sting.gatk.version");
                programRecord.setProgramVersion(version);
            } catch (MissingResourceException e) {
            }

            StringBuffer sb = new StringBuffer();
            sb.append(getToolkit().createApproximateCommandLineArgumentString(getToolkit(), this));
            sb.append(" Covariates=[");
            for (Covariate cov : requestedCovariates) {
                sb.append(cov.getClass().getSimpleName());
                sb.append(", ");
            }
            sb.setCharAt(sb.length() - 2, ']');
            sb.setCharAt(sb.length() - 1, ' ');
            programRecord.setCommandLine(sb.toString());

            List<SAMProgramRecord> oldRecords = header.getProgramRecords();
            List<SAMProgramRecord> newRecords = new ArrayList<SAMProgramRecord>(oldRecords.size() + 1);
            for (SAMProgramRecord record : oldRecords) {
                if (!record.getId().startsWith(PROGRAM_RECORD_NAME))
                    newRecords.add(record);
            }
            newRecords.add(programRecord);
            header.setProgramRecords(newRecords);

            // Write out the new header
            OUTPUT_BAM.writeHeader(header);
        }
    }

    /**
     * For each covariate read in a value and parse it. Associate those values with the data itself (num observation and num mismatches)
     *
     * @param line A line of CSV data read from the recalibration table data file
     */
    private void addCSVData(final File file, final String line) {
        final String[] vals = line.split(",");

        // Check if the data line is malformed, for example if the read group string contains a comma then it won't be parsed correctly
        if (vals.length != requestedCovariates.size() + 3) { // +3 because of nObservations, nMismatch, and Qempirical
            throw new UserException.MalformedFile(file, "Malformed input recalibration file. Found data line with too many fields: " + line +
                    " --Perhaps the read group string contains a comma and isn't being parsed correctly.");
        }

        final Object[] key = new Object[requestedCovariates.size()];
        Covariate cov;
        int iii;
        for (iii = 0; iii < requestedCovariates.size(); iii++) {
            cov = requestedCovariates.get(iii);
            key[iii] = cov.getValue(vals[iii]);
        }

        // Create a new datum using the number of observations, number of mismatches, and reported quality score
        final RecalDatum datum = new RecalDatum(Long.parseLong(vals[iii]), Long.parseLong(vals[iii + 1]), Double.parseDouble(vals[1]), 0.0);
        // Add that datum to all the collapsed tables which will be used in the sequential calculation
        dataManager.addToAllTables(key, datum, PRESERVE_QSCORES_LESS_THAN);
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * For each base in the read calculate a new recalibrated quality score and replace the quality scores in the read
     *
     * @param refBases References bases over the length of the read
     * @param read     The read to be recalibrated
     * @return The read with quality scores replaced
     */
    public SAMRecord map(ReferenceContext refBases, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {

        if (read.getReadLength() == 0) { // Some reads have '*' as the SEQ field and samtools returns length zero. We don't touch these reads.
            return read;
        }

        RecalDataManager.parseSAMRecord(read, RAC);

        byte[] originalQuals = read.getBaseQualities();
        final byte[] recalQuals = originalQuals.clone();

        final String platform = read.getReadGroup().getPlatform();
        if (platform.toUpperCase().contains("SOLID") && !(RAC.SOLID_RECAL_MODE == RecalDataManager.SOLID_RECAL_MODE.DO_NOTHING)) {
            if (!(RAC.SOLID_NOCALL_STRATEGY == RecalDataManager.SOLID_NOCALL_STRATEGY.THROW_EXCEPTION)) {
                final boolean badColor = RecalDataManager.checkNoCallColorSpace(read);
                if (badColor) {
                    numReadsWithMalformedColorSpace++;
                    if (RAC.SOLID_NOCALL_STRATEGY == RecalDataManager.SOLID_NOCALL_STRATEGY.LEAVE_READ_UNRECALIBRATED) {
                        return read; // can't recalibrate a SOLiD read with no calls in the color space, and the user wants to skip over them
                    }
                    else if (RAC.SOLID_NOCALL_STRATEGY == RecalDataManager.SOLID_NOCALL_STRATEGY.PURGE_READ) {
                        read.setReadFailsVendorQualityCheckFlag(true);
                        return read;
                    }
                }
            }
            originalQuals = RecalDataManager.calcColorSpace(read, originalQuals, RAC.SOLID_RECAL_MODE, refBases == null ? null : refBases.getBases());
        }

        //compute all covariate values for this read
        final Comparable[][] covariateValues_offset_x_covar = RecalDataManager.computeCovariates(read, requestedCovariates);

        // For each base in the read
        for (int offset = 0; offset < read.getReadLength(); offset++) {

            final Object[] fullCovariateKey = covariateValues_offset_x_covar[offset];

            Byte qualityScore = (Byte) qualityScoreByFullCovariateKey.get(fullCovariateKey);
            if (qualityScore == null) {
                qualityScore = performSequentialQualityCalculation(fullCovariateKey);
                qualityScoreByFullCovariateKey.put(qualityScore, fullCovariateKey);
            }

            recalQuals[offset] = qualityScore;
        }

        preserveQScores(originalQuals, recalQuals); // Overwrite the work done if original quality score is too low

        read.setBaseQualities(recalQuals); // Overwrite old qualities with new recalibrated qualities
        if (!DO_NOT_WRITE_OQ && read.getAttribute(RecalDataManager.ORIGINAL_QUAL_ATTRIBUTE_TAG) == null) { // Save the old qualities if the tag isn't already taken in the read
            try {
                read.setAttribute(RecalDataManager.ORIGINAL_QUAL_ATTRIBUTE_TAG, SAMUtils.phredToFastq(originalQuals));
            } catch (IllegalArgumentException e) {
                throw  new UserException.MalformedBAM(read, "illegal base quality encountered; " + e.getMessage());
            }
        }

        if (!skipUQUpdate && refBases != null && read.getAttribute(SAMTag.UQ.name()) != null) {
            read.setAttribute(SAMTag.UQ.name(), SequenceUtil.sumQualitiesOfMismatches(read, refBases.getBases(), read.getAlignmentStart() - 1, false));
        }

        if (RAC.SOLID_RECAL_MODE == RecalDataManager.SOLID_RECAL_MODE.SET_Q_ZERO_BASE_N && refBases != null && read.getAttribute(SAMTag.NM.name()) != null) {
            read.setAttribute(SAMTag.NM.name(), SequenceUtil.calculateSamNmTag(read, refBases.getBases(), read.getAlignmentStart() - 1, false));
        }

        return read;
    }

    /**
     * Implements a serial recalibration of the reads using the combinational table.
     * First, we perform a positional recalibration, and then a subsequent dinuc correction.
     *
     * Given the full recalibration table, we perform the following preprocessing steps:
     *
     * - calculate the global quality score shift across all data [DeltaQ]
     * - calculate for each of cycle and dinuc the shift of the quality scores relative to the global shift
     * -- i.e., DeltaQ(dinuc) = Sum(pos) Sum(Qual) Qempirical(pos, qual, dinuc) - Qreported(pos, qual, dinuc) / Npos * Nqual
     * - The final shift equation is:
     *
     * Qrecal = Qreported + DeltaQ + DeltaQ(pos) + DeltaQ(dinuc) + DeltaQ( ... any other covariate ... )
     *
     * @param key The list of Comparables that were calculated from the covariates
     * @return A recalibrated quality score as a byte
     */
    private byte performSequentialQualityCalculation(final Object... key) {

        final byte qualFromRead = (byte) Integer.parseInt(key[1].toString());
        final Object[] readGroupCollapsedKey = new Object[1];
        final Object[] qualityScoreCollapsedKey = new Object[2];
        final Object[] covariateCollapsedKey = new Object[3];

        // The global quality shift (over the read group only)
        readGroupCollapsedKey[0] = key[0];
        final RecalDatum globalRecalDatum = ((RecalDatum) dataManager.getCollapsedTable(0).get(readGroupCollapsedKey));
        double globalDeltaQ = 0.0;
        if (globalRecalDatum != null) {
            final double globalDeltaQEmpirical = globalRecalDatum.getEmpiricalQuality();
            final double aggregrateQReported = globalRecalDatum.getEstimatedQReported();
            globalDeltaQ = globalDeltaQEmpirical - aggregrateQReported;
        }

        // The shift in quality between reported and empirical
        qualityScoreCollapsedKey[0] = key[0];
        qualityScoreCollapsedKey[1] = key[1];
        final RecalDatum qReportedRecalDatum = ((RecalDatum) dataManager.getCollapsedTable(1).get(qualityScoreCollapsedKey));
        double deltaQReported = 0.0;
        if (qReportedRecalDatum != null) {
            final double deltaQReportedEmpirical = qReportedRecalDatum.getEmpiricalQuality();
            deltaQReported = deltaQReportedEmpirical - qualFromRead - globalDeltaQ;
        }

        // The shift in quality due to each covariate by itself in turn
        double deltaQCovariates = 0.0;
        double deltaQCovariateEmpirical;
        covariateCollapsedKey[0] = key[0];
        covariateCollapsedKey[1] = key[1];
        for (int iii = 2; iii < key.length; iii++) {
            covariateCollapsedKey[2] = key[iii]; // The given covariate
            final RecalDatum covariateRecalDatum = ((RecalDatum) dataManager.getCollapsedTable(iii).get(covariateCollapsedKey));
            if (covariateRecalDatum != null) {
                deltaQCovariateEmpirical = covariateRecalDatum.getEmpiricalQuality();
                deltaQCovariates += (deltaQCovariateEmpirical - qualFromRead - (globalDeltaQ + deltaQReported));
            }
        }

        final double newQuality = qualFromRead + globalDeltaQ + deltaQReported + deltaQCovariates;
        return QualityUtils.boundQual((int) Math.round(newQuality), (byte) MAX_QUALITY_SCORE);

        // Verbose printouts used to validate with old recalibrator
        //if(key.contains(null)) {
        //    System.out.println( key  + String.format(" => %d + %.2f + %.2f + %.2f + %.2f = %d",
        //                 qualFromRead, globalDeltaQ, deltaQReported, deltaQPos, deltaQDinuc, newQualityByte));
        //}
        //else {
        //    System.out.println( String.format("%s %s %s %s => %d + %.2f + %.2f + %.2f + %.2f = %d",
        //                 key.get(0).toString(), key.get(3).toString(), key.get(2).toString(), key.get(1).toString(), qualFromRead, globalDeltaQ, deltaQReported, deltaQPos, deltaQDinuc, newQualityByte) );
        //}

        //return newQualityByte;
    }

    /**
     * Loop over the list of qualities and overwrite the newly recalibrated score to be the original score if it was less than some threshold
     *
     * @param originalQuals The list of original base quality scores
     * @param recalQuals    A list of the new recalibrated quality scores
     */
    private void preserveQScores(final byte[] originalQuals, final byte[] recalQuals) {
        for (int iii = 0; iii < recalQuals.length; iii++) {
            if (originalQuals[iii] < PRESERVE_QSCORES_LESS_THAN) {
                recalQuals[iii] = originalQuals[iii];
            }
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Start the reduce with a handle to the output bam file
     *
     * @return A FileWriter pointing to a new bam file
     */
    public SAMFileWriter reduceInit() {
        return OUTPUT_BAM;
    }

    /**
     * Output each read to disk
     *
     * @param read   The read to output
     * @param output The FileWriter to write the read to
     * @return The FileWriter
     */
    public SAMFileWriter reduce(SAMRecord read, SAMFileWriter output) {
        if (output != null) {
            output.addAlignment(read);
        }
        return output;
    }

    /**
     * Do nothing
     *
     * @param output The SAMFileWriter that outputs the bam file
     */
    public void onTraversalDone(SAMFileWriter output) {
        if (numReadsWithMalformedColorSpace != 0) {
            if (RAC.SOLID_NOCALL_STRATEGY == RecalDataManager.SOLID_NOCALL_STRATEGY.LEAVE_READ_UNRECALIBRATED) {
                Utils.warnUser("Discovered " + numReadsWithMalformedColorSpace + " SOLiD reads with no calls in the color space. Unfortunately these reads cannot be recalibrated with this recalibration algorithm " +
                        "because we use reference mismatch rate as the only indication of a base's true quality. These reads have had reference bases inserted as a way of correcting " +
                        "for color space misalignments and there is now no way of knowing how often it mismatches the reference and therefore no way to recalibrate the quality score. " +
                        "These reads remain in the output bam file but haven't been corrected for reference bias. !!! USE AT YOUR OWN RISK !!!");
            }
            else if (RAC.SOLID_NOCALL_STRATEGY == RecalDataManager.SOLID_NOCALL_STRATEGY.PURGE_READ) {
                Utils.warnUser("Discovered " + numReadsWithMalformedColorSpace + " SOLiD reads with no calls in the color space. Unfortunately these reads cannot be recalibrated with this recalibration algorithm " +
                        "because we use reference mismatch rate as the only indication of a base's true quality. These reads have had reference bases inserted as a way of correcting " +
                        "for color space misalignments and there is now no way of knowing how often it mismatches the reference and therefore no way to recalibrate the quality score. " +
                        "These reads were completely removed from the output bam file.");

            }
        }
    }
}
