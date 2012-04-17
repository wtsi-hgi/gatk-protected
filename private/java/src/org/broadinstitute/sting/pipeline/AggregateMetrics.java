/*
 * Copyright (c) 2012, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
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

                            Object cvKey = cvTable.findPrimaryKeyByData(key);
                            if (cvKey == null) {
                                logger.warn(String.format("  Could not find CountVariants key %s in %s", Arrays.asList(key), evalFilename));
                                continue;
                            }

                            Object titvKey = titvTable.findPrimaryKeyByData(key);
                            if (titvKey == null) {
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
