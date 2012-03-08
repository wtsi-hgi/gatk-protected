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

import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.io.PrintStream;
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
 * A recalibration table file in CSV format that is used by the TableRecalibration walker.
 * It is a comma-separated text file relating the desired covariates to the number of such bases and their rate of mismatch in the genome, and its implied empirical quality score.
 *
 * The first 20 lines of such a file is shown below.
 * * The file begins with a series of comment lines describing:
 * ** The number of counted loci
 * ** The number of counted bases
 * ** The number of skipped loci and the fraction skipped, due to presence in dbSNP or bad reference bases
 *
 * * After the comments appears a header line indicating which covariates were used as well as the ordering of elements in the subsequent records.
 *
 * * After the header, data records occur one per line until the end of the file. The first several items on a line are the values of the individual covariates and will change
 * depending on which covariates were specified at runtime. The last three items are the data- that is, number of observations for this combination of covariates, number of
 * reference mismatches, and the raw empirical quality score calculated by phred-scaling the mismatch rate.
 *
 * <pre>
 * # Counted Sites    19451059
 * # Counted Bases    56582018
 * # Skipped Sites    82666
 * # Fraction Skipped 1 / 235 bp
 * ReadGroup,QualityScore,Cycle,Context,EventType,nObservations,nMismatches,Qempirical
 * SRR006446,11,65,CA,M,9,1,10
 * SRR006446,11,48,TA,M,10,0,40
 * SRR006446,11,67,AAGCTGAC,I,27,0,40
 * SRR006446,11,61,GAAAGCTG,D,11,1,10
 * SRR006446,12,34,CAGAAAGC,D,47,1,17
 * SRR006446,12,30,GAAAAATC,I,52,1,17
 * SRR006446,12,36,AA,M,352,1,25
 * SRR006446,12,17,TA,M,182,11,12
 * SRR006446,11,48,TG,M,2,0,40
 * SRR006446,11,67,AGTCCTGC,D,1,0,40
 * SRR006446,12,34,CGGTCCTG,I,9,0,40
 * SRR006446,12,30,GGGTCCTG,I,43,0,40
 * ERR001876,4,31,AGAAAAAT,D,1,0,40
 * ERR001876,4,31,AT,M,2,2,1
 * ERR001876,4,31,CA,M,1,0,40
 * </pre>
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
 *   -knownSites another/optional/setOfSitesToMask.vcf \
 *   -I my_reads.bam \
 *   -T BaseQualityScoreRecalibrator \
 *   -cov ReadGroupCovariate \
 *   -cov QualityScoreCovariate \
 *   -cov CycleCovariate \
 *   -cov ContextCovariate \
 *   -o my_reads.recal_data.csv
 * </pre>
 */

@BAQMode(ApplicationTime = BAQ.ApplicationTime.FORBIDDEN)
@By(DataSource.READS)
// Only look at covered loci, not every loci of the reference file
@ReadFilters({MappingQualityZeroFilter.class, MappingQualityUnavailableFilter.class})
// Filter out all reads with zero or unavailable mapping quality
@Requires({DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES})
// This walker requires both -I input.bam and -R reference.fasta
@PartitionBy(PartitionType.LOCUS)
public class BaseQualityScoreRecalibrator extends LocusWalker<BaseQualityScoreRecalibrator.CountedData, BaseQualityScoreRecalibrator.CountedData> implements TreeReducible<BaseQualityScoreRecalibrator.CountedData> {
    @ArgumentCollection
    private RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();                        // All the command line arguments for BQSR and it's covariates

    private BQSRKeyManager keyManager;                                                                          // The key manager to handle the BitSet key representation conversions
    private final HashMap<BitSet, RecalDatumOptimized> allCovariatesTable = new HashMap<BitSet, RecalDatumOptimized>(Short.MAX_VALUE);      // The big table that holds ALL the data for the covariates

    protected final ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>();                      // A list to hold the all the covariate objects that were requested (required + standard + experimental)
    protected final ArrayList<Covariate> requiredCovariates  = new ArrayList<Covariate>();                      // A list to hold the all the required covaraite objects that were requested
    protected final ArrayList<Covariate> optionalCovariates  = new ArrayList<Covariate>();                      // A list to hold the all the optional covariate objects that were requested 

    protected final String SKIP_RECORD_ATTRIBUTE = "SKIP";                                                      // used to label reads that should be skipped.
    protected final String SEEN_ATTRIBUTE = "SEEN";                                                             // used to label reads as processed.

    /**
     * Parse the -cov arguments and create a list of covariates to be used here
     * Based on the covariates' estimates for initial capacity allocate the data hashmap
     */
    public void initialize() {

        if (RAC.FORCE_PLATFORM != null) {
            RAC.DEFAULT_PLATFORM = RAC.FORCE_PLATFORM;
        }

        // Get a list of all available covariates
        final List<Class<? extends Covariate>> covariateClasses = new PluginManager<Covariate>(Covariate.class).getPlugins();
        final List<Class<? extends RequiredCovariate>> requiredClasses = new PluginManager<RequiredCovariate>(RequiredCovariate.class).getPlugins();
        final List<Class<? extends StandardCovariate>> standardClasses = new PluginManager<StandardCovariate>(StandardCovariate.class).getPlugins();

        // Print and exit if that's what was requested
        if (RAC.LIST_ONLY) {
            logger.info("Available covariates:");
            for (Class<?> covClass : covariateClasses) 
                logger.info(covClass.getSimpleName());            
            logger.info("");
            System.exit(0);                                                                                     // Early exit here because user requested it
        }
        

        // Warn the user if no dbSNP file or other variant mask was specified
        if (RAC.knownSites.isEmpty() && !RAC.RUN_WITHOUT_DBSNP) {
            throw new UserException.CommandLineException("This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.");
        }

        addRequiredCovariatesToList(requiredCovariates, requiredClasses);                                       // add the required covariates 
        if (RAC.USE_STANDARD_COVARIATES)                                                  
            addStandardCovariatesToList(optionalCovariates, standardClasses);                                   // add the standard covariates if -standard was specified by the user    
                                                                                                        
        if (RAC.COVARIATES != null) {                                                                           // parse the -cov arguments that were provided, skipping over the ones already specified
            for (String requestedCovariateString : RAC.COVARIATES) {
                boolean foundClass = false;
                for (Class<? extends Covariate> covClass : covariateClasses) {
                    if (requestedCovariateString.equalsIgnoreCase(covClass.getSimpleName())) {                  // -cov argument matches the class name for an implementing class
                        foundClass = true;
                        if (!requiredClasses.contains(covClass) && (!RAC.USE_STANDARD_COVARIATES || !standardClasses.contains(covClass))) {
                            try {
                                final Covariate covariate = covClass.newInstance();                             // Now that we've found a matching class, try to instantiate it
                                optionalCovariates.add(covariate);
                            } catch (Exception e) {
                                throw new DynamicClassResolutionException(covClass, e);
                            }
                        }
                    }
                }

                if (!foundClass) {
                    throw new UserException.CommandLineException("The requested covariate type (" + requestedCovariateString + ") isn't a valid covariate option. Use --list to see possible covariates.");
                }
            }
        }
        
        requestedCovariates.addAll(requiredCovariates);                                                         // add all the required covariates to the full list of requested covariates
        requestedCovariates.addAll(optionalCovariates);                                                         // add all the optional covariates to the full list of requested covariates
        
        logger.info("The covariates being used here: ");
        for (Covariate cov : requestedCovariates) {
            logger.info("\t" + cov.getClass().getSimpleName());
            cov.initialize(RAC);                                                                                // Initialize any covariate member variables using the shared argument collection
        }

        keyManager = new BQSRKeyManager(requiredCovariates, optionalCovariates);                                // Now that we know how many covariates are going to be in this run, initialize the Key Manager
    }

    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     *
     * @param tracker The reference metadata tracker
     * @param ref     The reference context
     * @param context The alignment context
     * @return Returns 1, but this value isn't used in the reduce step
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

                if (read.containsTemporaryAttribute(SKIP_RECORD_ATTRIBUTE))
                    continue;                                           // This read has been marked to be skipped.

                if (read.getBaseQualities()[offset] == 0)
                    continue;                                           // Skip this base if quality is zero

                if (p.isDeletion())
                    continue;                                           // We look at deletions from their left aligned base, not at the base itself

                if (!read.containsTemporaryAttribute(SEEN_ATTRIBUTE)) {
                    read.setTemporaryAttribute(SEEN_ATTRIBUTE, true);
                    RecalDataManager.parseSAMRecord(read, RAC);

                    if (!(RAC.SOLID_NOCALL_STRATEGY == RecalDataManager.SOLID_NOCALL_STRATEGY.THROW_EXCEPTION) && RecalDataManager.checkNoCallColorSpace(read)) {
                        read.setTemporaryAttribute(SKIP_RECORD_ATTRIBUTE, true);
                        continue;                                       // Skip over reads with no calls in the color space if the user requested it
                    }

                    RecalDataManager.parseColorSpace(read);
                    RecalDataManager.computeCovariates(read, requestedCovariates);
                }

                final byte[] bases = read.getReadBases();
                final byte refBase = ref.getBase();

                // SOLID bams have inserted the reference base into the read
                // if the color space in inconsistent with the read base so
                // skip it
                if (!ReadUtils.isSOLiDRead(read) || RAC.SOLID_RECAL_MODE == RecalDataManager.SOLID_RECAL_MODE.DO_NOTHING || !RecalDataManager.isInconsistentColorSpace(read, offset))
                    updateDataForPileupElement(counter, p, refBase);    // This base finally passed all the checks for a good base, so add it to the big data hashmap

                else {
                    if (refBase == bases[offset])
                        counter.solidInsertedReferenceBases++;          // calculate SOLID reference insertion rate
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
        logger.info("Writing raw recalibration data...");
        outputToCSV(sum);
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
    private void updateDataForPileupElement(CountedData counter, PileupElement pileupElement, final byte refBase) {
        final int offset = pileupElement.getOffset();
        final ReadCovariates readCovariates = RecalDataManager.covariateKeySetFrom(pileupElement.getRead());

        List<BitSet> mismatchesKeys = keyManager.bitSetsFromAllKeys(readCovariates.getMismatchesKeySet(offset), EventType.BASE_SUBSTITUTION);
        List<BitSet> insertionsKeys = keyManager.bitSetsFromAllKeys(readCovariates.getInsertionsKeySet(offset), EventType.BASE_INSERTION);
        List<BitSet> deletionsKeys  = keyManager.bitSetsFromAllKeys(readCovariates.getDeletionsKeySet(offset), EventType.BASE_DELETION);

        // the three arrays WON'T always have the same length
        for (BitSet key : mismatchesKeys)
            updateCovariateWithKeySet(key, !BaseUtils.basesAreEqual(pileupElement.getBase(), refBase));

        // negative strand reads should be check if the previous base is an insertion. Positive strand reads check the next base.
        for (BitSet key : insertionsKeys)
            updateCovariateWithKeySet(key, (pileupElement.getRead().getReadNegativeStrandFlag()) ? pileupElement.isAfterInsertion() : pileupElement.isBeforeInsertion());

        // negative strand reads should be check if the previous base is a deletion. Positive strand reads check the next base.
        for (BitSet key : deletionsKeys)
            updateCovariateWithKeySet(key, (pileupElement.getRead().getReadNegativeStrandFlag()) ? pileupElement.isAfterDeletion() : pileupElement.isBeforeDeletion());

        counter.countedBases++;
    }

    /**
     * Generic functionality to add to the number of observations and mismatches given a covariate key set
     *
     * @param hashKey the key to the hash map in bitset representation aggregating all the covariate keys and the event type
     * @param isError whether or not this base is an error (reference mismatch or precedes insertion or deletion)
     */
    private void updateCovariateWithKeySet(BitSet hashKey, boolean isError) {
        RecalDatumOptimized datum = allCovariatesTable.get(hashKey);            // Using the list of covariate values as a key, pick out the RecalDatum from the data HashMap
        if (datum == null) {                                                    // key doesn't exist yet in the map so make a new bucket and add it
            datum = new RecalDatumOptimized();
            allCovariatesTable.put(hashKey, datum);                             // initialized with zeros, will be incremented at end of method.
        }
        datum.increment(1, isError ? 1 : 0);                                    // Add one to the number of observations and potentially one to the number of mismatches
    }

    /**
     * For each entry (key-value pair) in the data hashmap output the Covariate's values as well as the RecalDatum's data in CSV format
     *
     * @param sum the counted data object to use as a reference for the number of bases/sites counted and skipped
     */
    private void outputToCSV(final CountedData sum) {
        PrintStream out = RAC.RECAL_FILE;
        out.printf("# Counted Sites    %d%n", sum.countedSites);
        out.printf("# Counted Bases    %d%n", sum.countedBases);
        out.printf("# Skipped Sites    %d%n", sum.skippedSites);
        out.printf("# Fraction Skipped 1 / %.0f bp%n", (double) sum.countedSites / sum.skippedSites);

        if (sum.solidInsertedReferenceBases != 0) {
            out.printf("# Fraction SOLiD inserted reference 1 / %.0f bases%n", (double) sum.countedBases / sum.solidInsertedReferenceBases);
            out.printf("# Fraction other color space inconsistencies 1 / %.0f bases%n", (double) sum.countedBases / sum.otherColorSpaceInconsistency);
        }
        
        ArrayList<String> requiredCovariateNames = new ArrayList<String>(requiredCovariates.size());
        ArrayList<String> optionalCovariateNames = new ArrayList<String>(optionalCovariates.size());

        for (Covariate covariate : requiredCovariates)
            requiredCovariateNames.add(covariate.getClass().getSimpleName().split("Covariate")[0]);
        
        for (Covariate covariate : optionalCovariates)
            optionalCovariateNames.add(covariate.getClass().getSimpleName().split("Covariate")[0]);

        out.println("# Required Covariates (in order): " + Utils.join(",", requiredCovariateNames));        // Print the required covariates used in order
        out.println("# Optional Covariates (in order): " + Utils.join(",", optionalCovariateNames));        // Print the optional covariates (standard + experimental) used in order
        out.println("# Recalibration Data  (in order): CovariateID,EventType,nObservations,nMismatches,Qempirical");

        for (Map.Entry<BitSet, RecalDatumOptimized> entry : allCovariatesTable.entrySet()) {
            for (Object key : keyManager.keySetFrom(entry.getKey()))
                out.print(key + ",");                                                                       // Print the keys (covariate values) plus the event type
            out.println(entry.getValue());                                                                  // Print the recalibration data
        }
        out.println("EOF");                                                                                 // Print an EOF marker
    }


    private void addRequiredCovariatesToList(List<Covariate> dest, List<Class<? extends RequiredCovariate>> classes) {
        if (classes.size() != 2) 
            throw new ReviewedStingException("The number of required covariate has changed, this is a hard change in the code and needs to be inspected");
        
        dest.add(new ReadGroupCovariate());             // enforce the order with RG first and QS next.
        dest.add(new QualityScoreCovariate());
    }

    private void addStandardCovariatesToList(List<Covariate> dest, List<Class<? extends StandardCovariate>> classes) {
        for (Class<?> covClass : classes) {
            try {
                final Covariate covariate = (Covariate) covClass.newInstance();
                dest.add(covariate);
            } catch (Exception e) {
                throw new DynamicClassResolutionException(covClass, e);
            }
        }
    }

    
    public class CountedData {
        private long countedSites = 0;                  // Number of loci used in the calculations, used for reporting in the output file
        private long countedBases = 0;                  // Number of bases used in the calculations, used for reporting in the output file
        private long skippedSites = 0;                  // Number of loci skipped because it was a dbSNP site, used for reporting in the output file
        private long solidInsertedReferenceBases = 0;   // Number of bases where we believe SOLID has inserted the reference because the color space is inconsistent with the read base
        private long otherColorSpaceInconsistency = 0;  // Number of bases where the color space is inconsistent with the read but the reference wasn't inserted.

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

