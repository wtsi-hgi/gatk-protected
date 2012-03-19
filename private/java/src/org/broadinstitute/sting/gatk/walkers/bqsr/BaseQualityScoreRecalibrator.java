/*
 * Copyright (c) 2010 The Broad Institute
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.recalibration.QualQuantizer;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

/**
 * First pass of the base quality score recalibration -- Generates recalibration table based on various user-specified covariates (such as reported quality score, cycle, and dinucleotide).
 *
 * <p>
 * This walker is designed to work as the first pass in a two-pass processing step. It does a by-locus traversal operating
 * only at sites that are not in dbSNP. We assume that all reference mismatches we see are therefore errors and indicative
 * of poor base quality. This walker generates tables based on various user-specified covariates (such as read group,
 * reported quality score, cycle, and dinucleotide). Since there is a large amount of data one can then calculate an empirical
 * probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations.
 * The output file is a CSV list of (the several covariate values, num observations, num mismatches, empirical quality score).
 * <p>
 * Note: ReadGroupCovariate and QualityScoreCovariate are required covariates and will be added for the user regardless of whether or not they were specified.
 *
 * <p>
 * See the GATK wiki for a tutorial and example recalibration accuracy plots.
 * http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
 *
 * <h2>Input</h2>
 * <p>
 * The input read data whose base quality scores need to be assessed.
 * <p>
 * A database of known polymorphic sites to skip over.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A GATK Report file with many tables:
 * <ol>
 *     <li>The list of arguments</li>
 *     <li>The quantized qualities table</li>
 *     <li>The recalibration table by read group</li>
 *     <li>The recalibration table by quality score</li>
 *     <li>The recalibration table for all the optional covariates</li>
 * </ol>
 *
 * The GATK Report is intended to be easy to read by humans or computers. Check out the documentation of the GATKReport to learn how to manipulate this table.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -T BaseQualityScoreRecalibrator \
 *   -I my_reads.bam \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
 *   -knownSites another/optional/setOfSitesToMask.vcf \
 *   -o recal_data.grp
 * </pre>
 */

@BAQMode(ApplicationTime = BAQ.ApplicationTime.FORBIDDEN)
@By(DataSource.READS)

@ReadFilters({MappingQualityZeroFilter.class, MappingQualityUnavailableFilter.class})                                   // only look at covered loci, not every loci of the reference file
@Requires({DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES})                                         // filter out all reads with zero or unavailable mapping quality
@PartitionBy(PartitionType.LOCUS)                                                                                       // this walker requires both -I input.bam and -R reference.fasta

public class BaseQualityScoreRecalibrator extends LocusWalker<BaseQualityScoreRecalibrator.CountedData, BaseQualityScoreRecalibrator.CountedData> implements TreeReducible<BaseQualityScoreRecalibrator.CountedData> {
    @ArgumentCollection
    private RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();                                // all the command line arguments for BQSR and it's covariates

    /**
     * Enables debug information in the GATK Report output (example: number of observations, number of mismatches, etc.)
     */
    @Argument(fullName = "debug_mode", shortName = "debug", doc = "enable debug information output", required = false)
    private boolean DEBUG_MODE = false;

    private Map<BQSRKeyManager, Map<BitSet, RecalDatum>> tablesAndKeysMap;                                              // a map mapping each recalibration table to its corresponding key manager
    private List<Byte> qualQuantizationMap;                                                                             // histogram containing the map for qual quantization (calculated after recalibration is done)

    protected final ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>();                              // list to hold the all the covariate objects that were requested (required + standard + experimental)
    protected ArrayList<Covariate> requiredCovariates;                                                                  // list to hold the all the required covaraite objects that were requested
    protected ArrayList<Covariate> optionalCovariates;                                                                  // list to hold the all the optional covariate objects that were requested

    private static final String SKIP_RECORD_ATTRIBUTE = "SKIP";                                                         // used to label reads that should be skipped.
    private static final String SEEN_ATTRIBUTE = "SEEN";                                                                // used to label reads as processed.
    private static final String COVARS_ATTRIBUTE = "COVARS";                                                            // used to store covariates array as a temporary attribute inside GATKSAMRecord.\
    private static final int SMOOTHING_CONSTANT = 1;                                                                    // used when calculating empirical qualities to avoid division by zero


    private static final String NO_DBSNP_EXCEPTION = "This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.";

    /**
     * Parse the -cov arguments and create a list of covariates to be used here
     * Based on the covariates' estimates for initial capacity allocate the data hashmap
     */
    public void initialize() {

        if (RAC.FORCE_PLATFORM != null)
            RAC.DEFAULT_PLATFORM = RAC.FORCE_PLATFORM;


        if (RAC.knownSites.isEmpty() && !RAC.RUN_WITHOUT_DBSNP)                                                         // Warn the user if no dbSNP file or other variant mask was specified
            throw new UserException.CommandLineException(NO_DBSNP_EXCEPTION);

        if (RAC.LIST_ONLY) {
            RecalDataManager.listAvailableCovariates(logger);
            System.exit(0);
        }

        Pair<ArrayList<Covariate>, ArrayList<Covariate>> covariates = RecalDataManager.initializeCovariates(RAC);       // initialize the required and optional covariates
        requiredCovariates = covariates.getFirst();                                                                     // assign them to our global variables
        optionalCovariates = covariates.getSecond();                                                                    // assign them to our global variables
        requestedCovariates.addAll(requiredCovariates);                                                                 // add all required covariates to the list of requested covariates
        requestedCovariates.addAll(optionalCovariates);                                                                 // add all optional covariates to the list of requested covariates

        logger.info("The covariates being used here: ");
        for (Covariate cov : requestedCovariates) {                                                                     // list all the covariates being used
            logger.info("\t" + cov.getClass().getSimpleName());
            cov.initialize(RAC);                                                                                        // initialize any covariate member variables using the shared argument collection
        }

        tablesAndKeysMap = RecalDataManager.initializeTables(requiredCovariates, optionalCovariates);                   // initialize the recalibration tables and their relative key managers
    }


    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     *
     * @param tracker the reference metadata tracker
     * @param ref     the reference context
     * @param context the alignment context
     * @return returns 1, but this value isn't used in the reduce step
     */
    public CountedData map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        // Only use data from non-dbsnp sites
        // Assume every mismatch at a non-dbsnp site is indicative of poor quality
        CountedData counter = new CountedData();

        // Only analyze sites not present in the provided known sites
        if (tracker.getValues(RAC.knownSites).size() == 0) {
            for (final PileupElement p : context.getBasePileup()) {
                final GATKSAMRecord read = p.getRead();
                final int offset = p.getOffset();

                if (read.containsTemporaryAttribute(SKIP_RECORD_ATTRIBUTE) || read.getBaseQualities()[offset] == 0)
                    continue;                                                                                           // This read has been marked to be skipped or has quality 0.

                if (!read.containsTemporaryAttribute(SEEN_ATTRIBUTE)) {
                    read.setTemporaryAttribute(SEEN_ATTRIBUTE, true);
                    RecalDataManager.parseSAMRecord(read, RAC);

                    if (!(RAC.SOLID_NOCALL_STRATEGY == RecalDataManager.SOLID_NOCALL_STRATEGY.THROW_EXCEPTION) && RecalDataManager.checkNoCallColorSpace(read)) {
                        read.setTemporaryAttribute(SKIP_RECORD_ATTRIBUTE, true);
                        continue;                                                                                       // Skip over reads with no calls in the color space if the user requested it
                    }

                    RecalDataManager.parseColorSpace(read);
                    read.setTemporaryAttribute(COVARS_ATTRIBUTE, RecalDataManager.computeCovariates(read, requestedCovariates));
                }

                final byte[] bases = read.getReadBases();
                final byte refBase = ref.getBase();

                // SOLID bams have inserted the reference base into the read
                // if the color space in inconsistent with the read base so
                // skip it
                if (!ReadUtils.isSOLiDRead(read) || RAC.SOLID_RECAL_MODE == RecalDataManager.SOLID_RECAL_MODE.DO_NOTHING || !RecalDataManager.isInconsistentColorSpace(read, offset))
                    updateDataForPileupElement(counter, p, refBase);                                                    // This base finally passed all the checks for a good base, so add it to the big data hashmap

                else {
                    if (refBase == bases[offset])
                        counter.solidInsertedReferenceBases++;                                                          // calculate SOLID reference insertion rate
                    else
                        counter.otherColorSpaceInconsistency++;
                }
            }
            counter.countedSites++;
        }
        else
            counter.skippedSites++;

        return counter;
    }

    /**
     * Initialize the reduce step by creating a PrintStream from the filename specified as an argument to the walker.
     *
     * @return returns A PrintStream created from the -recalFile filename argument specified to the walker
     */
    public CountedData reduceInit() {
        return new CountedData();
    }

    /**
     * The Reduce method doesn't do anything for this walker.
     *
     * @param mapped Result of the map. This value is immediately ignored.
     * @param sum    The summing CountedData used to output the CSV data
     * @return returns The sum used to output the CSV data
     */
    public CountedData reduce(CountedData mapped, CountedData sum) {
        return sum.add(mapped);
    }

    public CountedData treeReduce(CountedData sum1, CountedData sum2) {
        return sum1.add(sum2);
    }

    /**
     * Write out the full data hashmap to disk in CSV format
     *
     * @param sum The CountedData to write out to RECAL_FILE
     */
    public void onTraversalDone(CountedData sum) {
        logger.info("Calculating empirical quality scores...");
        calculateEmpiricalQuals();
        logger.info("Calculating quantized quality scores...");
        quantizeQualityScores();
        logger.info("Writing GATK Report...");
        generateReport();
        logger.info("...done!");
    }


    /**
     * Major workhorse routine for this walker.
     * Loop through the list of requested covariates and pick out the value from the read, offset, and reference
     * Using the list of covariate values as a key, pick out the RecalDatum and increment,
     * adding one to the number of observations and potentially one to the number of mismatches for all three
     * categories (mismatches, insertions and deletions).
     *
     * @param counter       Data structure which holds the counted bases
     * @param pileupElement The pileup element to update
     * @param refBase       The reference base at this locus
     */
    private void updateDataForPileupElement(final CountedData counter, final PileupElement pileupElement, final byte refBase) {
        final int offset = pileupElement.getOffset();
        final ReadCovariates readCovariates = covariateKeySetFrom(pileupElement.getRead());

        final RecalDatum mismatchesDatum = createDatumObject(pileupElement.getQual(), !BaseUtils.basesAreEqual(pileupElement.getBase(), refBase));
        final RecalDatum insertionsDatum = createDatumObject(pileupElement.getQual(), (pileupElement.getRead().getReadNegativeStrandFlag()) ? pileupElement.isAfterInsertion() : pileupElement.isBeforeInsertion());
        final RecalDatum deletionsDatum  = createDatumObject(pileupElement.getQual(), (pileupElement.getRead().getReadNegativeStrandFlag()) ? pileupElement.isAfterDeletion() : pileupElement.isBeforeDeletion());

        for (Map.Entry<BQSRKeyManager, Map<BitSet, RecalDatum>> entry : tablesAndKeysMap.entrySet()) {
            BQSRKeyManager keyManager = entry.getKey();
            Map<BitSet, RecalDatum> table = entry.getValue();

            List<BitSet> mismatchesKeys = keyManager.bitSetsFromAllKeys(readCovariates.getMismatchesKeySet(offset), EventType.BASE_SUBSTITUTION);
            List<BitSet> insertionsKeys = keyManager.bitSetsFromAllKeys(readCovariates.getInsertionsKeySet(offset), EventType.BASE_INSERTION);
            List<BitSet> deletionsKeys  = keyManager.bitSetsFromAllKeys(readCovariates.getDeletionsKeySet(offset), EventType.BASE_DELETION);


            for (BitSet key : mismatchesKeys)
                updateCovariateWithKeySet(table, key, mismatchesDatum);                                                 // the three arrays WON'T always have the same length

            for (BitSet key : insertionsKeys)
                updateCovariateWithKeySet(table, key, insertionsDatum);                                                 // negative strand reads should be check if the previous base is an insertion. Positive strand reads check the next base.

            for (BitSet key : deletionsKeys)
                updateCovariateWithKeySet(table, key, deletionsDatum);                                                  // negative strand reads should be check if the previous base is a deletion. Positive strand reads check the next base.
        }
        counter.countedBases++;
    }

    /**
     * creates a datum object with one observation and one or zero error 
     * 
     * @param reportedQual  the quality score reported by the instrument for this base
     * @param isError       whether or not the observation is an error
     * @return a new RecalDatum object with the observation and the error
     */
    private RecalDatum createDatumObject(byte reportedQual, boolean isError) {
        return new RecalDatum(1, isError ? 1:0, reportedQual, 0.0);                                                     // initialized with zeros, will be incremented at end of method.
    }

    /**
     * Generic functionality to add to the number of observations and mismatches given a covariate key set
     *
     * @param recalTable the recalibration table to use
     * @param hashKey the key to the hash map in bitset representation aggregating all the covariate keys and the event type
     * @param datum the RecalDatum object with the observation/error information 
     */
    private void updateCovariateWithKeySet(Map<BitSet, RecalDatum> recalTable, BitSet hashKey, RecalDatum datum) {
        RecalDatum previousDatum = recalTable.get(hashKey);                                                             // using the list of covariate values as a key, pick out the RecalDatum from the data HashMap
        if (previousDatum == null)                                                                                      // key doesn't exist yet in the map so make a new bucket and add it
            recalTable.put(hashKey, datum);                                                                             // initialized with zeros, will be incremented at end of method.
        else
            previousDatum.increment(datum);                                                                             // add one to the number of observations and potentially one to the number of mismatches
    }


    /**
     * Get the covariate key set from a read
     *
     * @param read the read
     * @return the covariate keysets for this read
     */
    private ReadCovariates covariateKeySetFrom(GATKSAMRecord read) {
        return (ReadCovariates) read.getTemporaryAttribute(COVARS_ATTRIBUTE);
    }

    /**
     * Calculates the empirical qualities in all recalibration tables
     */
    private void calculateEmpiricalQuals() {
        for (Map<BitSet, RecalDatum> table : tablesAndKeysMap.values()) {
            for (RecalDatum datum : table.values()) {
                datum.calcCombinedEmpiricalQuality(SMOOTHING_CONSTANT, QualityUtils.MAX_QUAL_SCORE);
                datum.calcEstimatedReportedQuality();
            }
        }
    }

    /**
     * go through the quality score table and use the # observations and the empirical quality score 
     * to build a quality score histogram for quantization. Then use the QuantizeQual algorithm to 
     * generate a quantization map (recalibrated_qual -> quantized_qual) 
     */
    private void quantizeQualityScores() {
        final Long [] qualHistogram = new Long[QualityUtils.MAX_QUAL_SCORE+1];                                          // create a histogram with the empirical quality distribution
        for (int i = 0; i < qualHistogram.length; i++)
            qualHistogram[i] = 0L;

        Map<BitSet, RecalDatum> qualTable = null;                                                                       // look for the quality score table                                                                      
        for (Map.Entry<BQSRKeyManager, Map<BitSet, RecalDatum>> entry : tablesAndKeysMap.entrySet()) {
            BQSRKeyManager keyManager = entry.getKey();
            if (keyManager.getRequiredCovariates().size() == 2)                                                         // it should be the only one with 2 required covaraites
                qualTable = entry.getValue();
        }
        
        if (qualTable == null)
            throw new ReviewedStingException("Could not find QualityScore table.");
        
        for (RecalDatum datum : qualTable.values()) {
            int empiricalQual = (int) Math.round(datum.getEmpiricalQuality());                                          // convert the empirical quality to an integer ( it is already capped by MAX_QUAL )
            long nObservations = datum.numObservations;
            qualHistogram[empiricalQual] += nObservations;                                                              // add the number of observations for every key
        }
        QualQuantizer quantizer = new QualQuantizer(Arrays.asList(qualHistogram), RAC.QUANTIZING_LEVELS, QualityUtils.MIN_USABLE_Q_SCORE);
        qualQuantizationMap = quantizer.getOriginalToQuantizedMap();
    }

    private void generateReport() {
        GATKReport report = new GATKReport();

        GATKReportTable argumentsTable = new GATKReportTable("Arguments", "Recalibration argument collection values used in this run");
        argumentsTable.addPrimaryKey("Argument");
        argumentsTable.addColumn("Value", "null");
        argumentsTable.set("covariate", "Value", (RAC.COVARIATES == null) ? "null" : Utils.join(",", RAC.COVARIATES));
        argumentsTable.set("standard_covs", "Value", RAC.USE_STANDARD_COVARIATES);
        argumentsTable.set("run_without_dbsnp", "Value", RAC.RUN_WITHOUT_DBSNP);
        argumentsTable.set("solid_recal_mode", "Value", RAC.SOLID_RECAL_MODE);
        argumentsTable.set("solid_nocall_strategy", "Value", RAC.SOLID_NOCALL_STRATEGY);
        argumentsTable.set("mismatches_context_size", "Value", RAC.MISMATCHES_CONTEXT_SIZE);
        argumentsTable.set("insertions_context_size", "Value", RAC.INSERTIONS_CONTEXT_SIZE);
        argumentsTable.set("deletions_context_size", "Value", RAC.DELETIONS_CONTEXT_SIZE);
        argumentsTable.set("mismatches_default_quality", "Value", RAC.MISMATCHES_DEFAULT_QUALITY);
        argumentsTable.set("insertions_default_quality", "Value", RAC.INSERTIONS_DEFAULT_QUALITY);
        argumentsTable.set("low_quality_tail", "Value", RAC.LOW_QUAL_TAIL);
        argumentsTable.set("default_platform", "Value", RAC.DEFAULT_PLATFORM);
        argumentsTable.set("force_platform", "Value", RAC.FORCE_PLATFORM);
        argumentsTable.set("quantizing_levels", "Value", RAC.QUANTIZING_LEVELS);
        report.addTable(argumentsTable);
        
        GATKReportTable quantizedTable = new GATKReportTable("Quantized", "Quality quantization map");
        quantizedTable.addPrimaryKey("QualityScore");
        quantizedTable.addColumn("QuantizedScore", (byte) 0);
        int i = 0;
        for (byte qual : qualQuantizationMap)
            quantizedTable.set(i++, "QuantizedScore", qual);
        report.addTable(quantizedTable);

        int tableIndex = 0;
        for (Map.Entry<BQSRKeyManager, Map<BitSet, RecalDatum>> entry : tablesAndKeysMap.entrySet()) {
            BQSRKeyManager keyManager = entry.getKey();
            Map<BitSet, RecalDatum> recalTable = entry.getValue();

            GATKReportTable reportTable = new GATKReportTable("RecalTable" + tableIndex++, "");
            final Pair<String, String> covariateValue     = new Pair<String, String>("CovariateValue", "%s");
            final Pair<String, String> covariateName      = new Pair<String, String>("CovariateName", "%s");
            final Pair<String, String> eventType          = new Pair<String, String>("EventType", "%s");
            final Pair<String, String> empiricalQuality   = new Pair<String, String>("EmpiricalQuality", "%.2f");
            final Pair<String, String> estimatedQReported = new Pair<String, String>("EstimatedQReported", "%.2f");

            long primaryKey = 0L;

            List<Covariate> requiredList = keyManager.getRequiredCovariates();                                          // ask the key manager what required covariates were used in this recal table 
            List<Covariate> optionalList = keyManager.getOptionalCovariates();                                          // ask the key manager what optional covariates were used in this recal table

            ArrayList<Pair<String, String>> columnNames = new ArrayList<Pair<String, String>>();                        // initialize the array to hold the column names

            for (Covariate covariate : requiredList) {
                String name = covariate.getClass().getSimpleName().split("Covariate")[0];                               // get the covariate names and put them in order
                columnNames.add(new Pair<String,String>(name, "%s"));                                                   // save the required covariate name so we can reference it in the future
            }

            if (optionalList.size() > 0) {
                columnNames.add(covariateValue);
                columnNames.add(covariateName);
            }
            
            columnNames.add(eventType);                                                                                 // the order of these column names is important here
            columnNames.add(empiricalQuality);
            columnNames.add(estimatedQReported);

            if (DEBUG_MODE) {
                final Pair<String, String> nObservations   = new Pair<String, String>("Observations", "%d");
                final Pair<String, String> nErrors = new Pair<String, String>("Errors", "%d");
                columnNames.add(nObservations);
                columnNames.add(nErrors);
            }


            reportTable.addPrimaryKey("PrimaryKey", false);                                                             // every table must have a primary key (hidden)
            for (Pair<String, String> columnName : columnNames)
                reportTable.addColumn(columnName.getFirst(), true, columnName.getSecond());                             // every table must have the event type

            for (Map.Entry<BitSet, RecalDatum> recalTableEntry : recalTable.entrySet()) {                               // create a map with column name => key value for all covariate keys
                BitSet bitSetKey = recalTableEntry.getKey();
                Map<String, Object> columnData = new HashMap<String, Object>(columnNames.size());
                Iterator<Pair<String, String>> iterator = columnNames.iterator();
                for (Object key : keyManager.keySetFrom(bitSetKey)) {                                                   
                    String columnName = iterator.next().getFirst();
                    columnData.put(columnName, key);                                                                    
                }                
                RecalDatum datum = recalTableEntry.getValue();
                columnData.put(iterator.next().getFirst(), datum.getEmpiricalQuality());                                // iterator.next() gives the column name for Empirical Quality
                columnData.put(iterator.next().getFirst(), Math.round(datum.getEstimatedQReported()));                  // iterator.next() gives the column name for EstimatedQReported

                if (DEBUG_MODE) {
                    columnData.put(iterator.next().getFirst(), datum.numObservations);
                    columnData.put(iterator.next().getFirst(), datum.numMismatches);
                }

                for (Map.Entry<String, Object> dataEntry : columnData.entrySet()) {
                    String columnName = dataEntry.getKey();
                    Object value = dataEntry.getValue();
                    reportTable.set(primaryKey, columnName, value.toString());
                }
                primaryKey++;
            }
            report.addTable(reportTable);
        }

        report.print(RAC.RECAL_FILE);
    }

    public class CountedData {
        private long countedSites = 0;                                                                                  // Number of loci used in the calculations, used for reporting in the output file
        private long countedBases = 0;                                                                                  // Number of bases used in the calculations, used for reporting in the output file
        private long skippedSites = 0;                                                                                  // Number of loci skipped because it was a dbSNP site, used for reporting in the output file
        private long solidInsertedReferenceBases = 0;                                                                   // Number of bases where we believe SOLID has inserted the reference because the color space is inconsistent with the read base
        private long otherColorSpaceInconsistency = 0;                                                                  // Number of bases where the color space is inconsistent with the read but the reference wasn't inserted.

        /**
         * Adds the values of other to this, returning this
         *
         * @param other another object
         * @return this object with the other object incremented
         */
        public CountedData add(CountedData other) {
            countedSites += other.countedSites;
            countedBases += other.countedBases;
            skippedSites += other.skippedSites;
            solidInsertedReferenceBases += other.solidInsertedReferenceBases;
            otherColorSpaceInconsistency += other.otherColorSpaceInconsistency;
            return this;
        }
    }


}

