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

package org.broadinstitute.sting.analyzecovariates;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.gatk.walkers.recalibration.Covariate;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.sting.utils.R.RScriptExecutor;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.io.Resource;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * Call R scripts to plot residual error versus the various covariates.
 *
 * <p>
 * After counting covariates in either the initial BAM File or again in the recalibrated BAM File, an analysis tool is available which
 * reads the .csv file and outputs several PDF (and .dat) files for each read group in the given BAM. These PDF files graphically
 * show the various metrics and characteristics of the reported quality scores (often in relation to the empirical qualities).
 * In order to show that any biases in the reported quality scores have been generally fixed through recalibration one should run
 * CountCovariates again on a bam file produced by TableRecalibration. In this way users can compare the analysis plots generated
 * by pre-recalibration and post-recalibration .csv files. Our usual chain of commands that we use to generate plots of residual
 * error is: CountCovariates, TableRecalibrate, samtools index on the recalibrated bam file, CountCovariates again on the recalibrated
 * bam file, and then AnalyzeCovariates on both the before and after recal_data.csv files to see the improvement in recalibration.
 *
 * <p>
 * The color coding along with the RMSE is included in the plots to give some indication of the number of observations that went into
 * each of the quality score estimates. It is defined as follows for N, the number of observations:
 *
 * <ul>
 * <li>light blue means N < 1,000</li>
 * <li>cornflower blue means 1,000 <= N < 10,000</li>
 * <li>dark blue means N >= 10,000</li>
 * <li>The pink dots indicate points whose quality scores are special codes used by the aligner and which are mathematically
 * meaningless and so aren't included in any of the numerical calculations.</li>
 * </ul>
 *
 * <p>
 * NOTE: Rscript needs to be in your environment PATH (this is the scripting version of R, not the interactive version).
 * See <a target="r-project" href="http://www.r-project.org">http://www.r-project.org</a> for more info on how to download and install R.
 *
 * <p>
 * See the GATK wiki for a tutorial and example recalibration accuracy plots.
 * <a target="gatkwiki" href="http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration"
 * >http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration</a>
 *
 * <h2>Input</h2>
 * <p>
 * The recalibration table file in CSV format that was generated by the CountCovariates walker.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx4g -jar AnalyzeCovariates.jar \
 *   -recalFile /path/to/recal.table.csv  \
 *   -outputDir /path/to/output_dir/  \
 *   -ignoreQ 5
 * </pre>
 *
 */

@DocumentedGATKFeature(
        groupName = "AnalyzeCovariates",
        summary = "Package to plot residual accuracy versus error covariates for the base quality score recalibrator")
public class AnalyzeCovariates extends CommandLineProgram {
    final private static Logger logger = Logger.getLogger(AnalyzeCovariates.class);

    private static final String PLOT_RESDIUAL_ERROR_QUALITY_SCORE_COVARIATE = "plot_residualError_QualityScoreCovariate.R";
    private static final String PLOT_RESDIUAL_ERROR_OTHER_COVARIATE = "plot_residualError_OtherCovariate.R";
    private static final String PLOT_INDEL_QUALITY_RSCRIPT = "plot_indelQuality.R";

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    /**
     * After the header, data records occur one per line until the end of the file. The first several items on a line are the
     * values of the individual covariates and will change depending on which covariates were specified at runtime. The last
     * three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches,
     * and the raw empirical quality score calculated by phred-scaling the mismatch rate.
     */
    @Input(fullName = "recal_file", shortName = "recalFile", doc = "The input recal csv file to analyze", required = false)
    private String RECAL_FILE = "output.recal_data.csv";
    @Argument(fullName = "output_dir", shortName = "outputDir", doc = "The directory in which to output all the plots and intermediate data files", required = false)
    private File OUTPUT_DIR = new File("analyzeCovariates");
    @Argument(fullName = "ignoreQ", shortName = "ignoreQ", doc = "Ignore bases with reported quality less than this number.", required = false)
    private int IGNORE_QSCORES_LESS_THAN = 5;
    @Argument(fullName = "numRG", shortName = "numRG", doc = "Only process N read groups. Default value: -1 (process all read groups)", required = false)
    private int NUM_READ_GROUPS_TO_PROCESS = -1; // -1 means process all read groups

    /**
     * Combinations of covariates in which there are zero mismatches technically have infinite quality. We get around this situation
     * by capping at the specified value. We've found that Q40 is too low when using a more completely database of known variation like dbSNP build 132 or later.
     */
    @Argument(fullName="max_quality_score", shortName="maxQ", required = false, doc="The integer value at which to cap the quality scores, default is 50")
    private int MAX_QUALITY_SCORE = 50;

    /**
     * This argument is useful for comparing before/after plots and you want the axes to match each other.
     */
    @Argument(fullName="max_histogram_value", shortName="maxHist", required = false, doc="If supplied, this value will be the max value of the histogram plots")
    private int MAX_HISTOGRAM_VALUE = 0;

    @Hidden
    @Argument(fullName="do_indel_quality", shortName="indels", required = false, doc="If supplied, do indel quality plotting")
    private boolean DO_INDEL_QUALITY = false;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private AnalysisDataManager dataManager; // Holds the data HashMap, mostly used by TableRecalibrationWalker to create collapsed data hashmaps
    private ArrayList<Covariate> requestedCovariates; // List of covariates to be used in this calculation
    private final Pattern COMMENT_PATTERN = Pattern.compile("^#.*");
    private final Pattern OLD_RECALIBRATOR_HEADER = Pattern.compile("^rg,.*");
    private final Pattern COVARIATE_PATTERN = Pattern.compile("^ReadGroup,QualityScore,.*");
    protected static final String EOF_MARKER = "EOF";

    protected int execute() {

        // create the output directory where all the data tables and plots will go
        if (!OUTPUT_DIR.exists() && !OUTPUT_DIR.mkdirs())
            throw new UserException.BadArgumentValue("--output_dir/-outDir", "Unable to create output directory: " + OUTPUT_DIR);

        if (!RScriptExecutor.RSCRIPT_EXISTS)
            Utils.warnUser(logger, "Rscript not found in environment path. Plots will not be generated.");

        // initialize all the data from the csv file and allocate the list of covariates
        logger.info("Reading in input csv file...");
        initializeData();
        logger.info("...Done!");

        // output data tables for Rscript to read in
        logger.info("Writing out intermediate tables for R...");
        writeDataTables();
        logger.info("...Done!");

        // perform the analysis using Rscript and output the plots
        logger.info("Calling analysis R scripts and writing out figures...");
        callRScripts();
        logger.info("...Done!");

        return 0;
    }

    private void initializeData() {

        // Get a list of all available covariates
        Collection<Class<? extends Covariate>> classes = new PluginManager<Covariate>(Covariate.class).getPlugins();

        int lineNumber = 0;
        boolean foundAllCovariates = false;

        // Read in the covariates that were used from the input file
        requestedCovariates = new ArrayList<Covariate>();

        try {
            for ( final String line : new XReadLines(new File( RECAL_FILE )) ) {
                lineNumber++;
                if( COMMENT_PATTERN.matcher(line).matches() || OLD_RECALIBRATOR_HEADER.matcher(line).matches() || line.equals(EOF_MARKER) )  {
                    ; // Skip over the comment lines, (which start with '#')
                }
                else if( COVARIATE_PATTERN.matcher(line).matches() ) { // The line string is either specifying a covariate or is giving csv data
                    if( foundAllCovariates ) {
                        throw new RuntimeException( "Malformed input recalibration file. Found covariate names intermingled with data in file: " + RECAL_FILE );
                    } else { // Found the covariate list in input file, loop through all of them and instantiate them
                        String[] vals = line.split(",");
                        for( int iii = 0; iii < vals.length - 3; iii++ ) { // There are n-3 covariates. The last three items are nObservations, nMismatch, and Qempirical
                            boolean foundClass = false;
                            for( Class<?> covClass : classes ) {
                                if( (vals[iii] + "Covariate").equalsIgnoreCase( covClass.getSimpleName() ) ) {
                                    foundClass = true;
                                    try {
                                        Covariate covariate = (Covariate)covClass.newInstance();
                                        requestedCovariates.add( covariate );
                                    } catch (Exception e) {
                                        throw new DynamicClassResolutionException(covClass, e);
                                    }
                                }
                            }

                            if( !foundClass ) {
                                throw new RuntimeException( "Malformed input recalibration file. The requested covariate type (" + (vals[iii] + "Covariate") + ") isn't a valid covariate option." );
                            }
                        }

                    }

                } else { // Found a line of data
                    if( !foundAllCovariates ) {

                        foundAllCovariates = true;

                        // At this point all the covariates should have been found and initialized
                        if( requestedCovariates.size() < 2 ) {
                            throw new RuntimeException( "Malformed input recalibration file. Covariate names can't be found in file: " + RECAL_FILE );
                        }

                        // Initialize any covariate member variables using the shared argument collection
                        for( Covariate cov : requestedCovariates ) {
                            cov.initialize( new RecalibrationArgumentCollection() );
                        }

                        // Initialize the data hashMaps
                        dataManager = new AnalysisDataManager( requestedCovariates.size() );

                    }
                    addCSVData(line); // Parse the line and add the data to the HashMap
                }
            }

        } catch ( FileNotFoundException e ) {
            throw new RuntimeException("Can not find input file: " + RECAL_FILE);
        } catch ( NumberFormatException e ) {
            throw new RuntimeException("Error parsing recalibration data at line " + lineNumber + ". Perhaps your table was generated by an older version of CovariateCounterWalker.");
        }
    }

    private void addCSVData(String line) {
        String[] vals = line.split(",");

        // Check if the data line is malformed, for example if the read group string contains a comma then it won't be parsed correctly
        if( vals.length != requestedCovariates.size() + 3 ) { // +3 because of nObservations, nMismatch, and Qempirical
            throw new RuntimeException("Malformed input recalibration file. Found data line with too many fields: " + line +
                    " --Perhaps the read group string contains a comma and isn't being parsed correctly.");
        }

        Object[] key = new Object[requestedCovariates.size()];
        Covariate cov;
        int iii;
        for( iii = 0; iii < requestedCovariates.size(); iii++ ) {
            cov = requestedCovariates.get( iii );
            key[iii] = cov.getValue( vals[iii] );
        }
        // Create a new datum using the number of observations, number of mismatches, and reported quality score
        final RecalDatum datum = new RecalDatum( Long.parseLong( vals[iii] ), Long.parseLong( vals[iii + 1] ), Double.parseDouble( vals[1] ), 0.0 );
        // Add that datum to all the collapsed tables which will be used in the sequential calculation
        dataManager.addToAllTables( key, datum, IGNORE_QSCORES_LESS_THAN );
    }

    private void writeDataTables() {

        int numReadGroups = 0;

        // for each read group
        for( final Object readGroupKey : dataManager.getCollapsedTable(0).data.keySet() ) {

            if( NUM_READ_GROUPS_TO_PROCESS == -1 || ++numReadGroups <= NUM_READ_GROUPS_TO_PROCESS ) {
                final String readGroup = readGroupKey.toString();
                final RecalDatum readGroupDatum = (RecalDatum) dataManager.getCollapsedTable(0).data.get(readGroupKey);
                logger.info(String.format(
                        "Writing out data tables for read group: %s\twith %s observations\tand aggregate residual error = %.3f",
                        readGroup, readGroupDatum.getNumObservations(),
                        readGroupDatum.empiricalQualDouble(0, MAX_QUALITY_SCORE) - readGroupDatum.getEstimatedQReported()));

                // for each covariate
                for( int iii = 1; iii < requestedCovariates.size(); iii++ ) {
                    Covariate cov = requestedCovariates.get(iii);

                    // Create a PrintStream
                    File outputFile = new File(OUTPUT_DIR, readGroup + "." + cov.getClass().getSimpleName()+ ".dat");
                    PrintStream output;
                    try {
                        output = new PrintStream(FileUtils.openOutputStream(outputFile));
                    } catch (IOException e) {
                        throw new UserException.CouldNotCreateOutputFile(outputFile, e);
                    }

                    try {
                        // Output the header
                        output.println("Covariate\tQreported\tQempirical\tnMismatches\tnBases");

                        for( final Object covariateKey : ((Map)dataManager.getCollapsedTable(iii).data.get(readGroupKey)).keySet() ) {
                            output.print( covariateKey.toString() + "\t" );                                                     // Covariate
                            final RecalDatum thisDatum = (RecalDatum)((Map)dataManager.getCollapsedTable(iii).data.get(readGroupKey)).get(covariateKey);
                            output.print( String.format("%.3f", thisDatum.getEstimatedQReported()) + "\t" );                    // Qreported
                            output.print( String.format("%.3f", thisDatum.empiricalQualDouble(0, MAX_QUALITY_SCORE)) + "\t" );  // Qempirical
                            output.print( thisDatum.getNumMismatches() + "\t" );                                                // nMismatches
                            output.println( thisDatum.getNumObservations() );                                                   // nBases
                        }
                    } finally {
                        // Close the PrintStream
                        IOUtils.closeQuietly(output);
                    }
                }
            } else {
                break;
            }

        }
    }

    private void callRScripts() {
        int numReadGroups = 0;

        // for each read group
        for( Object readGroupKey : dataManager.getCollapsedTable(0).data.keySet() ) {
            if(++numReadGroups <= NUM_READ_GROUPS_TO_PROCESS || NUM_READ_GROUPS_TO_PROCESS == -1) {

                String readGroup = readGroupKey.toString();
                logger.info("Analyzing read group: " + readGroup);

                // for each covariate
                for( int iii = 1; iii < requestedCovariates.size(); iii++ ) {
                    final Covariate cov = requestedCovariates.get(iii);
                    final File outputFile = new File(OUTPUT_DIR, readGroup + "." + cov.getClass().getSimpleName()+ ".dat");
                    if (DO_INDEL_QUALITY) {
                        RScriptExecutor executor = new RScriptExecutor();
                        executor.addScript(new Resource(PLOT_INDEL_QUALITY_RSCRIPT, AnalyzeCovariates.class));
                        // The second argument is the name of the covariate in order to make the plots look nice
                        executor.addArgs(outputFile, cov.getClass().getSimpleName().split("Covariate")[0]);
                        executor.exec();
                    } else {
                        if( iii == 1 ) {
                            // Analyze reported quality
                            RScriptExecutor executor = new RScriptExecutor();
                            executor.addScript(new Resource(PLOT_RESDIUAL_ERROR_QUALITY_SCORE_COVARIATE, AnalyzeCovariates.class));
                            // The second argument is the Q scores that should be turned pink in the plot because they were ignored
                            executor.addArgs(outputFile, IGNORE_QSCORES_LESS_THAN, MAX_QUALITY_SCORE, MAX_HISTOGRAM_VALUE);
                            executor.exec();
                        } else { // Analyze all other covariates
                            RScriptExecutor executor = new RScriptExecutor();
                            executor.addScript(new Resource(PLOT_RESDIUAL_ERROR_OTHER_COVARIATE, AnalyzeCovariates.class));
                            // The second argument is the name of the covariate in order to make the plots look nice
                            executor.addArgs(outputFile, cov.getClass().getSimpleName().split("Covariate")[0]);
                            executor.exec();
                        }
                    }
                }
            } else { // at the maximum number of read groups so break out
                break;
            }
        }
    }

    public static void main(String args[]) {
        try {
            AnalyzeCovariates clp = new AnalyzeCovariates();
            start(clp, args);
            System.exit(CommandLineProgram.result);
        } catch (Exception e) {
            exitSystemWithError(e);
        }
    }
}
