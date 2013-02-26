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

package org.broadinstitute.sting.pipeline;

import edu.mit.broad.picard.genotype.concordance.DbSnpMatchMetrics;
import edu.mit.broad.picard.genotype.fingerprint.v2.FingerprintingSummaryMetrics;
import net.sf.picard.analysis.AlignmentSummaryMetrics;
import net.sf.picard.analysis.InsertSizeMetrics;
import net.sf.picard.analysis.directed.HsMetrics;
import net.sf.picard.metrics.MetricBase;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.util.TabbedTextFileWithHeaderParser;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.util.CloseableIterator;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.utils.R.RUtils;
import org.broadinstitute.sting.utils.classloader.JVMUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * Aggregates Picard metrics with variant eval GATK reports into a single tab delimited file easily imported into R.
 */
public class AggregateMetrics extends CommandLineProgram {
    @Input(doc = "One or more TSV files", shortName = "tsv")
    public List<File> tsv_file;

    @Input(doc = "One or more eval files", shortName = "eval")
    public List<File> eval_file;

    @Output(doc = "Output file", shortName = "out", required = false)
    public File output = null;

    @Override
    protected int execute() throws Exception {
        if (tsv_file.size() != eval_file.size())
            throw new UserException.BadArgumentValue("tsv_file/eval_file",
                    String.format("Different number of tsvs and eval files specified: %d tsvs, %d evals", tsv_file.size(), eval_file.size()));

        PrintStream outputStream;
        if (output != null)
            outputStream = new PrintStream(FileUtils.openOutputStream(output));
        else
            outputStream = System.out;

        try {
            // print out headers
            List<String> headers = new ArrayList<String>();
            headers.addAll(Arrays.asList("project_name", "squid_project", "squid_sample", "aggregation_version"));
            headers.addAll(Arrays.asList(COUNT_VARIANTS_COLUMNS));
            headers.addAll(Arrays.asList(TITV_VARIANT_EVALUATOR_COLUMNS));
            headers.addAll(getFullMetricsFields());
            headers.add("sequencing_dates");

            printColumns(outputStream, headers);

            int fileCount = tsv_file.size();
            for (int i = 0; i < fileCount; i++) {
                File tsvFilename = tsv_file.get(i);
                File evalFilename = eval_file.get(i);

                logger.info(String.format("processing %s %s", tsvFilename, evalFilename));

                if (!tsvFilename.exists()) {
                    logger.warn(String.format("tsv file %s does not exist", tsvFilename));
                    continue;
                }

                if (!evalFilename.exists()) {
                    logger.warn(String.format("eval file %s does not exist", evalFilename));
                    continue;
                }

                String projectName = FilenameUtils.removeExtension(tsvFilename.getName());

                GATKReport report = new GATKReport(evalFilename);
                GATKReportTable cvTable = report.getTable("CountVariants");
                GATKReportTable titvTable = report.getTable("TiTvVariantEvaluator");

                for (PicardSample picardSample : PicardAggregationUtils.parseSamples(tsvFilename, false)) {
                    String squidProject = picardSample.getProject();
                    String squidSample = picardSample.getSample();
                    int aggregationVersion = picardSample.getVersion();

                    logger.info(String.format("  processing squid project = %s, sample = %s, version = %d", squidProject, squidSample, aggregationVersion));

                    // NOTE: This data is duplicated into the row for every functional class / novelty combination.
                    List<Object> fullMetrics = getFullMetrics(squidProject, squidSample, aggregationVersion);
                    List<Date> sequencingDates = getBamRunDates(squidProject, squidSample, aggregationVersion);

                    String sequencingDatesString = RUtils.toDateList(sequencingDates);

                    for (String functionalClass : FUNCTIONAL_CLASSES) {
                        for (String novelty : NOVELTIES) {
                            List<Object> columns = new ArrayList<Object>();

                            columns.addAll(Arrays.asList(projectName, squidProject, squidSample, aggregationVersion));

                            Object[] key = new Object[]{"dbsnp", "eval", functionalClass, novelty, squidSample};

                            int cvKey = cvTable.findRowByData(key);
                            if (cvKey == -1) {
                                logger.warn(String.format("  Could not find CountVariants key %s in %s", Arrays.asList(key), evalFilename));
                                continue;
                            }

                            int titvKey = titvTable.findRowByData(key);
                            if (titvKey == -1) {
                                logger.warn(String.format("  Could not find TiTvVariantEvaluator key %s in %s", Arrays.asList(key), evalFilename));
                                continue;
                            }

                            for (String column : COUNT_VARIANTS_COLUMNS)
                                columns.add(cvTable.get(cvKey, column));

                            for (String column : TITV_VARIANT_EVALUATOR_COLUMNS)
                                columns.add(titvTable.get(titvKey, column));

                            columns.addAll(fullMetrics);
                            columns.add(sequencingDatesString);

                            // replace any null in the columns with the text "NA"
                            Collections.replaceAll(columns, null, "NA");

                            printColumns(outputStream, columns);
                        }
                    }
                }
            }

        } finally {
            if (output != null)
                IOUtils.closeQuietly(outputStream);
        }
        return 0;
    }

    public static void main(String[] args) {
        try {
            AggregateMetrics instance = new AggregateMetrics();
            start(instance, args);
            System.exit(CommandLineProgram.result);
        } catch (UserException e) {
            exitSystemWithUserError(e);
        } catch (Exception e) {
            exitSystemWithError(e);
        }
    }

    private static final Logger logger = Logger.getLogger(AggregateMetrics.class);

    private static final String[] FUNCTIONAL_CLASSES = {"all", "missense", "nonsense", "silent"};
    private static final String[] NOVELTIES = {"all", "known", "novel"};

    private static final String[] COUNT_VARIANTS_COLUMNS = {
            "CompRod", "EvalRod", "FunctionalClass", "Novelty", "nProcessedLoci", "nCalledLoci", "nRefLoci", "nVariantLoci", "variantRate", "variantRatePerBp",
            "nSNPs", "nMNPs", "nInsertions", "nDeletions", "nComplex", "nNoCalls", "nHets", "nHomRef", "nHomVar", "nSingletons", "nHomDerived", "heterozygosity",
            "heterozygosityPerBp", "hetHomRatio", "indelRate", "indelRatePerBp", "deletionInsertionRatio"};

    private static final String[] TITV_VARIANT_EVALUATOR_COLUMNS = {
            "nTi", "nTv", "tiTvRatio", "nTiInComp", "nTvInComp", "TiTvRatioStandard", "nTiDerived", "nTvDerived", "tiTvDerivedRatio"};

    private static final MetricsType[] SAMPLE_SUMMARY_METRICS_TYPES = {
            new MetricsType(HsMetrics.class, "hybrid_selection_metrics", null, null),
            new MetricsType(AlignmentSummaryMetrics.class, "alignment_summary_metrics", "CATEGORY=PAIR", null),
            new MetricsType(InsertSizeMetrics.class, "insert_size_metrics", "PAIR_ORIENTATION=FR", "FR"),
            new MetricsType(InsertSizeMetrics.class, "insert_size_metrics", "PAIR_ORIENTATION=RF", "RF"),
            new MetricsType(InsertSizeMetrics.class, "insert_size_metrics", "PAIR_ORIENTATION=TANDEM", "TANDEM"),
            new MetricsType(DbSnpMatchMetrics.class, "dbsnp_matches", null, null)};

    private static List<String> getFullMetricsFields() {
        ArrayList<String> fields = new ArrayList<String>();

        // Add the initiative name from analysis_files.txt
        fields.add("INITIATIVE");
        // Add in the fingerprint lods.
        fields.add("FINGERPRINT_LODS");
        fields.add("HAPLOTYPES_CONFIDENTLY_MATCHING");

        for (MetricsType metricsType : SAMPLE_SUMMARY_METRICS_TYPES) {
            List<String> metricsFields = MetricsUtils.getFieldNames(metricsType.getType());
            if (metricsType.getSuffix() != null) {
                int numFields = metricsFields.size();
                for (int i = 0; i < numFields; i++) {
                    metricsFields.set(i, metricsFields.get(i) + "_" + metricsType.getSuffix());
                }
            }
            fields.addAll(metricsFields);
        }
        return fields;
    }

    private static List<Date> getBamRunDates(String project, String sample, int aggregationVersion) {
        SAMFileReader reader = new SAMFileReader(new File(PicardAggregationUtils.getSampleBam(project, sample, aggregationVersion)));
        ArrayList<Date> dates = new ArrayList<Date>();
        try {
            reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            SAMFileHeader header = reader.getFileHeader();
            for (SAMReadGroupRecord record : header.getReadGroups()) {
                if (!dates.contains(record.getRunDate()))
                    dates.add(record.getRunDate());
            }
        } finally {
            IOUtils.closeQuietly(reader);
        }
        return dates;
    }

    @SuppressWarnings("unchecked")
    private static List<Object> getFullMetrics(String squidProject, String squidSample, int aggregationVersion) throws IOException {
        ArrayList<Object> data = new ArrayList<Object>();

        String initiative = "NA";
        String fingerprintLodsString = "NA";
        String haplotypesConfidentlyMatchingString = "NA";

        // Load in the initiative from analysis_files.txt. I believe this data is lane-level, so we grab only the first row for the initiative data.
        File analysisFile = new File(PicardAggregationUtils.getSampleDir(squidProject, squidSample, aggregationVersion), "analysis_files.txt");
        CloseableIterator<TabbedTextFileWithHeaderParser.Row> itor = new TabbedTextFileWithHeaderParser(analysisFile).iterator();
        try {
            if (itor.hasNext())
                initiative = itor.next().getField("INITIATIVE");
        } finally {
            itor.close();
        }

        File fingerprintMetricsFile = new File(PicardAggregationUtils.getSampleFile(squidProject, squidSample, aggregationVersion, "fingerprinting_summary_metrics"));
        if (!fingerprintMetricsFile.exists()) {
            logger.warn(String.format("metrics file %s does not exist", fingerprintMetricsFile));
        } else {
            List<FingerprintingSummaryMetrics> fingerprintingSummaryMetrics = (List<FingerprintingSummaryMetrics>) MetricsFile.readBeans(fingerprintMetricsFile);
            List<Double> fingerprintLods = new ArrayList<Double>();
            List<Integer> haplotypesConfidentlyMatching = new ArrayList<Integer>();
            for (FingerprintingSummaryMetrics metric : fingerprintingSummaryMetrics) {
                fingerprintLods.add(metric.LOD_EXPECTED_SAMPLE);
                haplotypesConfidentlyMatching.add(metric.HAPLOTYPES_CONFIDENTLY_MATCHING);
            }
            fingerprintLodsString = RUtils.toNumberList(fingerprintLods);
            haplotypesConfidentlyMatchingString = RUtils.toNumberList(haplotypesConfidentlyMatching);
        }

        data.add(initiative);
        data.add(fingerprintLodsString);
        data.add(haplotypesConfidentlyMatchingString);

        for (MetricsType metricsType : SAMPLE_SUMMARY_METRICS_TYPES) {
            String metricsExtension = metricsType.getExtension();
            String metricsFilter = metricsType.getFilter();

            File metricsFile = new File(PicardAggregationUtils.getSampleFile(squidProject, squidSample, aggregationVersion, metricsExtension));
            MetricBase metrics;
            if (!metricsFile.exists()) {
                logger.warn(String.format("metrics file %s does not exist", metricsFile));
                metrics = null;
            } else {
                List<? extends MetricBase> filteredMetrics = MetricsFile.readBeans(metricsFile);
                if (metricsFilter != null) {
                    String[] tokens = metricsFilter.split("=", 2);
                    filteredMetrics = MetricsUtils.filterMetrics(filteredMetrics, tokens[0], tokens[1]);
                }
                metrics = MetricsUtils.getSampleMetric(filteredMetrics);
            }

            List<String> fields = MetricsUtils.getFieldNames(metricsType.getType());
            for (String field : fields) {
                if (metrics != null) {
                    // TODO: Optimize expensive field lookup.
                    data.add(JVMUtils.getFieldValue(JVMUtils.findField(metrics.getClass(), field), metrics));
                } else {
                    data.add("NA");
                }
            }
        }

        return data;
    }

    private static void printColumns(PrintStream stream, List<?> columns) {
        stream.println(StringUtils.join(columns, "\t"));
    }
}
