package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 8/3/11
 * Time: 2:24 PM
 */
public class SlidingWindow {

    // Sliding Window data
    protected LinkedList<GATKSAMRecord> readsInWindow;
    protected LinkedList<HeaderElement> windowHeader;
    protected int contextSizeMismatches;
    protected int contextSizeIndels;
    protected int contextSize;                                      // the largest context size (between mismatches and indels)
    protected int startLocation;
    protected int stopLocation;
    protected String contig;
    protected int contigIndex;
    protected SAMFileHeader header;
    protected GATKSAMReadGroupRecord readGroupAttribute;
    protected int downsampleCoverage;

    // Running consensus data
    protected SyntheticRead runningConsensus;
    protected int consensusCounter;
    protected String consensusReadName;

    // Filtered Data Consensus data
    protected SyntheticRead filteredDataConsensus;
    protected int filteredDataConsensusCounter;
    protected String filteredDataReadName;


    // Additional parameters
    protected double MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT;    // proportion has to be greater than this value to trigger variant region due to mismatches
    protected double MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT;  // proportion has to be greater than this value to trigger variant region due to deletions
    protected int MIN_BASE_QUAL_TO_COUNT;                           // qual has to be greater than or equal to this value
    protected int MIN_MAPPING_QUALITY;

    protected ReduceReadsWalker.DownsampleStrategy downsampleStrategy;
    
    public int getStopLocation() {
        return stopLocation;
    }

    public String getContig() {
        return contig;
    }

    public int getContigIndex() {
        return contigIndex;
    }


    public SlidingWindow(String contig, int contigIndex, int contextSizeMismatches, int contextSizeIndels, SAMFileHeader header, GATKSAMReadGroupRecord readGroupAttribute, int windowNumber, final double minAltProportionToTriggerVariant, final double minIndelProportionToTriggerVariant, int minBaseQual, int minMappingQuality, int downsampleCoverage, final ReduceReadsWalker.DownsampleStrategy downsampleStrategy) {
        this.startLocation = -1;
        this.stopLocation = -1;
        this.contextSizeMismatches = contextSizeMismatches;
        this.contextSizeIndels = contextSizeIndels;
        this.contextSize = Math.max(contextSizeIndels, contextSizeMismatches);
        this.downsampleCoverage = downsampleCoverage;

        this.MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT = minAltProportionToTriggerVariant;
        this.MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT = minIndelProportionToTriggerVariant;
        this.MIN_BASE_QUAL_TO_COUNT = minBaseQual;
        this.MIN_MAPPING_QUALITY = minMappingQuality;

        this.windowHeader = new LinkedList<HeaderElement>();
        this.readsInWindow = new LinkedList<GATKSAMRecord>();

        this.contig = contig;
        this.contigIndex = contigIndex;
        this.header = header;
        this.readGroupAttribute = readGroupAttribute;

        this.consensusCounter = 0;
        this.consensusReadName = "Consensus-" + windowNumber + "-";

        this.filteredDataConsensusCounter = 0;
        this.filteredDataReadName = "Filtered-" + windowNumber + "-";

        this.runningConsensus = null;
        this.filteredDataConsensus = null;
        
        this.downsampleStrategy = downsampleStrategy;
    }

    /**
     * Add a read to the sliding window and slides the window accordingly.
     * <p/>
     * Reads are assumed to be in order, therefore, when a read is added the sliding window can
     * assume that no more reads will affect read.getUnclippedStart() - contextSizeMismatches. The window
     * slides forward to that position and returns all reads that may have been finalized in the
     * sliding process.
     *
     * @param read
     * @return a list of reads that have been finished by sliding the window.
     */
    public List<GATKSAMRecord> addRead(GATKSAMRecord read) {
        // If this is the first read in the window, update startLocation
        if (startLocation < 0)
            startLocation = read.getAlignmentStart();

        List<GATKSAMRecord> finalizedReads = slideIfPossible(read.getUnclippedStart());

        // update the window header counts
        updateHeaderCounts(read);

        // add read to sliding reads
        readsInWindow.add(read);

        return finalizedReads;
    }

    /**
     * Determines if the window can be slid given the new incoming read.
     * <p/>
     * We check from the start of the window to the (unclipped) start of the new incoming read if there
     * is any variant.
     * <p/>
     * If there are no variant sites, we slide the left side of the window to the unclipped start of the
     * new incoming read minus the context size (we can only guarantee that the reads before this point
     * will not be in any variant region ever)
     * <p/>
     * If there are variant sites, we check if it's time to close the variant region.
     *
     * @param incomingReadUnclippedStart the incoming read's start position. Must be the unclipped start!
     * @return all reads that have fallen to the left of the sliding window after the slide
     */
    protected List<GATKSAMRecord> slideIfPossible(int incomingReadUnclippedStart) {

        List<GATKSAMRecord> finalizedReads = new LinkedList<GATKSAMRecord>();

        // No point doing any calculation if the read's unclipped start is behind
        // the start of the window
        if (incomingReadUnclippedStart - contextSize > startLocation) {
            int positionShift = incomingReadUnclippedStart - startLocation;
            boolean[] variantSite = markSites(startLocation + positionShift);

            // this is the limit of what we can close/send to consensus (non-inclusive)
            int breakpoint = Math.max(positionShift - contextSize, 0);
            int start = 0;
            int i = 0;

            while (i < breakpoint) {
                while (i < breakpoint && !variantSite[i]) i++;
                finalizedReads.addAll(addToSyntheticReads(start, i));

                start = i;
                while (i < breakpoint && variantSite[i]) i++;
                if (i < breakpoint || (i - 1 > 0 && variantSite[i - 1] && !variantSite[i])) {
                    finalizedReads.addAll(closeVariantRegion(start, i - 1));
                    start = i;
                }
            }

            slide(start);
        }

        return finalizedReads;
    }

    /**
     * returns an array marked with variant and non-variant regions (it uses
     * markVariantRegions to make the marks)
     *
     * @param stop check the window from start to stop (not-inclusive)
     * @return a boolean array with 'true' marking variant regions and false marking consensus sites
     */
    protected boolean[] markSites(int stop) {

        boolean[] markedSites = new boolean[stop - startLocation + contextSize + 1];

        Iterator<HeaderElement> headerElementIterator = windowHeader.iterator();
        for (int i = startLocation; i < stop; i++) {
            if (headerElementIterator.hasNext()) {
                HeaderElement headerElement = headerElementIterator.next();
                int context = 0;

                if (headerElement.isVariantFromDeletions() || headerElement.isVariantFromInsertions())
                    context = Math.max(context, contextSizeIndels);

                // if context has been set to indels, only do mismatches if contextSize is bigger
                if ((context > 0 && contextSizeMismatches > contextSizeIndels) || headerElement.isVariantFromMismatches())
                    context = Math.max(context, contextSizeMismatches);

                if (context > 0)
                    markVariantRegion(markedSites, i - startLocation, context);

            } else
                break;
        }
        return markedSites;
    }

    /**
     * Marks the sites around the variant site (as true)
     *
     * @param markedSites         the boolean array to bear the marks
     * @param variantSiteLocation the location where a variant site was found
     * @param context             the number of bases to mark on each side of the variant site
     */
    protected void markVariantRegion(boolean[] markedSites, int variantSiteLocation, int context) {
        int from = (variantSiteLocation < context) ? 0 : variantSiteLocation - context;
        int to = (variantSiteLocation + context + 1 > markedSites.length) ? markedSites.length : variantSiteLocation + context + 1;
        for (int i = from; i < to; i++)
            markedSites[i] = true;
    }

    /**
     * Adds bases to the running consensus or filtered data accordingly
     * <p/>
     * If adding a sequence with gaps, it will finalize multiple consensus reads and keep the last running consensus
     *
     * @param start the first header index to add to consensus
     * @param end   the first header index NOT TO add to consensus
     * @return a list of consensus reads generated by this call. Empty list if no consensus was generated.
     */
    protected List<GATKSAMRecord> addToSyntheticReads(int start, int end) {
        LinkedList<GATKSAMRecord> reads = new LinkedList<GATKSAMRecord>();
        if (start < end) {

            ListIterator<HeaderElement> headerElementIterator = windowHeader.listIterator(start);

            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException(String.format("Requested to add to synthetic reads a region that contains no header element at index: %d  - %d / %d", start, windowHeader.size(), end));

            HeaderElement headerElement = headerElementIterator.next();

            if (headerElement.hasConsensusData()) {
                reads.addAll(finalizeAndAdd(ConsensusType.FILTERED));

                int endOfConsensus = findNextNonConsensusElement(start, end);
                addToRunningConsensus(start, endOfConsensus);

                if (endOfConsensus <= start)
                    throw new ReviewedStingException(String.format("next start is <= current start: (%d <= %d)", endOfConsensus, start));

                reads.addAll(addToSyntheticReads(endOfConsensus, end));
            } else if (headerElement.hasFilteredData()) {
                reads.addAll(finalizeAndAdd(ConsensusType.CONSENSUS));

                int endOfFilteredData = findNextNonFilteredDataElement(start, end);
                addToFilteredData(start, endOfFilteredData);

                if (endOfFilteredData <= start)
                    throw new ReviewedStingException(String.format("next start is <= current start: (%d <= %d)", endOfFilteredData, start));

                reads.addAll(addToSyntheticReads(endOfFilteredData, end));
            } else if (headerElement.isEmpty()) {
                reads.addAll(finalizeAndAdd(ConsensusType.BOTH));

                int endOfEmptyData = findNextNonEmptyElement(start, end);

                if (endOfEmptyData <= start)
                    throw new ReviewedStingException(String.format("next start is <= current start: (%d <= %d)", endOfEmptyData, start));

                reads.addAll(addToSyntheticReads(endOfEmptyData, end));
            } else
                throw new ReviewedStingException(String.format("Header Element %d is neither Consensus, Data or Empty. Something is wrong.", start));

        }

        return reads;
    }

    /**
     * Finalizes one or more synthetic reads.
     *
     * @param type the synthetic reads you want to close
     * @return the GATKSAMRecords generated by finalizing the synthetic reads
     */
    private List<GATKSAMRecord> finalizeAndAdd(ConsensusType type) {
        GATKSAMRecord read = null;
        List<GATKSAMRecord> list = new LinkedList<GATKSAMRecord>();

        switch (type) {
            case CONSENSUS:
                read = finalizeRunningConsensus();
                break;
            case FILTERED:
                read = finalizeFilteredDataConsensus();
                break;
            case BOTH:
                read = finalizeRunningConsensus();
                if (read != null) list.add(read);
                read = finalizeFilteredDataConsensus();
        }
        if (read != null) list.add(read);

        return list;
    }

    /**
     * Looks for the next position without consensus data
     *
     * @param start beginning of the filtered region
     * @param upTo  limit to search for another consensus element
     * @return next position with consensus data or empty
     */
    private int findNextNonConsensusElement(int start, int upTo) {
        Iterator<HeaderElement> headerElementIterator = windowHeader.listIterator(start);
        int index = start;
        while (index < upTo) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("There are no more header elements in this window");

            HeaderElement headerElement = headerElementIterator.next();
            if (!headerElement.hasConsensusData())
                break;
            index++;
        }
        return index;
    }

    /**
     * Looks for the next position without filtered data
     *
     * @param start beginning of the region
     * @param upTo  limit to search for
     * @return next position with no filtered data
     */
    private int findNextNonFilteredDataElement(int start, int upTo) {
        Iterator<HeaderElement> headerElementIterator = windowHeader.listIterator(start);
        int index = start;
        while (index < upTo) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("There are no more header elements in this window");

            HeaderElement headerElement = headerElementIterator.next();
            if (!headerElement.hasFilteredData() || headerElement.hasConsensusData())
                break;
            index++;
        }
        return index;
    }

    /**
     * Looks for the next non-empty header element
     *
     * @param start beginning of the region
     * @param upTo  limit to search for
     * @return next position with non-empty element
     */
    private int findNextNonEmptyElement(int start, int upTo) {
        ListIterator<HeaderElement> headerElementIterator = windowHeader.listIterator(start);
        int index = start;
        while (index < upTo) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("There are no more header elements in this window");

            HeaderElement headerElement = headerElementIterator.next();
            if (!headerElement.isEmpty())
                break;
            index++;
        }
        return index;
    }


    /**
     * Adds bases to the filtered data synthetic read.
     * <p/>
     * Different from the addToConsensus method, this method assumes a contiguous sequence of filteredData
     * bases.
     *
     * @param start the first header index to add to consensus
     * @param end   the first header index NOT TO add to consensus
     * @return a list of consensus reads generated by this call. Empty list if no consensus was generated.
     */
    private void addToFilteredData(int start, int end) {
        if (filteredDataConsensus == null)
            filteredDataConsensus = new SyntheticRead(header, readGroupAttribute, contig, contigIndex, filteredDataReadName + filteredDataConsensusCounter++, windowHeader.get(start).location, GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG);

        ListIterator<HeaderElement> headerElementIterator = windowHeader.listIterator(start);
        for (int index = start; index < end; index++) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("Requested to create a filtered data synthetic read from " + start + " to " + end + " but " + index + " does not exist");

            HeaderElement headerElement = headerElementIterator.next();
            if (headerElement.hasConsensusData())
                throw new ReviewedStingException("Found consensus data inside region to add to filtered data.");

            if (!headerElement.hasFilteredData())
                throw new ReviewedStingException("No filtered data in " + index);

            BaseIndex base = headerElement.filteredBaseCounts.baseIndexWithMostCounts();
            byte count = (byte) Math.min(headerElement.filteredBaseCounts.countOfMostCommonBase(), Byte.MAX_VALUE);
            byte qual = headerElement.filteredBaseCounts.averageQualsOfMostCommonBase();
            filteredDataConsensus.add(base, count, qual, headerElement.getRMS());
        }
    }

    /**
     * Adds bases to the filtered data synthetic read.
     * <p/>
     * Different from the addToConsensus method, this method assumes a contiguous sequence of filteredData
     * bases.
     *
     * @param start the first header index to add to consensus
     * @param end   the first header index NOT TO add to consensus
     * @return a list of consensus reads generated by this call. Empty list if no consensus was generated.
     */
    private void addToRunningConsensus(int start, int end) {
        if (runningConsensus == null)
            runningConsensus = new SyntheticRead(header, readGroupAttribute, contig, contigIndex, consensusReadName + consensusCounter++, windowHeader.get(start).location, GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG);

        Iterator<HeaderElement> headerElementIterator = windowHeader.listIterator(start);
        for (int index = start; index < end; index++) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("Requested to create a running consensus synthetic read from " + start + " to " + end + " but " + index + " does not exist");

            HeaderElement headerElement = headerElementIterator.next();
            if (!headerElement.hasConsensusData())
                throw new ReviewedStingException("No CONSENSUS data in " + index);

            BaseIndex base = headerElement.consensusBaseCounts.baseIndexWithMostCounts();
            byte count = (byte) Math.min(headerElement.consensusBaseCounts.countOfMostCommonBase(), Byte.MAX_VALUE);
            byte qual = headerElement.consensusBaseCounts.averageQualsOfMostCommonBase();
            runningConsensus.add(base, count, qual, headerElement.getRMS());
        }
    }


    /**
     * Finalizes a variant region, any adjacent synthetic reads.
     *
     * @param start the first position in the variant region (inclusive)
     * @param stop  the last position of the variant region (inclusive)
     * @return all reads contained in the variant region plus any adjacent synthetic reads
     */
    @Requires("start <= stop")
    protected List<GATKSAMRecord> closeVariantRegion(int start, int stop) {
        List<GATKSAMRecord> syntheticReads = new LinkedList<GATKSAMRecord>();
        List<GATKSAMRecord> allReads = new LinkedList<GATKSAMRecord>();

        syntheticReads.addAll(finalizeAndAdd(ConsensusType.BOTH));  // Finalize consensus if there is one to finalize before the Variant Region

        int refStart = windowHeader.get(start).location;            // Clipping operations are reference based, not read based
        int refStop = windowHeader.get(stop).location;

        for (GATKSAMRecord read : readsInWindow) {
            GATKSAMRecord trimmedRead = ReadClipper.hardClipToRegion(read, refStart, refStop);
            if (trimmedRead.getReadLength() > 0) {
                allReads.add(trimmedRead);
            }
        }

        List<GATKSAMRecord> result = (downsampleCoverage > 0) ? downsampleVariantRegion(allReads, refStart, refStop) : allReads;
        result.addAll(syntheticReads);
        return result;                                              // finalized reads will be downsampled if necessary
    }


    /**
     * Downsamples a variant region to the downsample coverage of the sliding window.
     *
     * It will use the downsampling strategy defined by the SlidingWindow
     *
     * @param allReads the reads to select from (all reads that cover the window)
     * @param refStart start of the window (inclusive)
     * @param refStop  end of the window (inclusive)
     * @return a list of reads selected by the downsampler to cover the window to at least the desired coverage
     */
    protected List<GATKSAMRecord> downsampleVariantRegion(final List<GATKSAMRecord> allReads, final int refStart, final int refStop) {
        LinkedList<GATKSAMRecord> readList = new LinkedList<GATKSAMRecord>();

        Pair<HashMap<Integer, HashSet<GATKSAMRecord>>, HashMap<GATKSAMRecord, Boolean[]>> mappings = ReadUtils.getBothReadToLociMappings(allReads, refStart, refStop);
        int [] coverageDistribution = ReadUtils.getCoverageDistributionOfReads(allReads, refStart, refStop); // the full coverage distribution array for the variant region with all the reads
        switch (downsampleStrategy) {
            case Normal:
                readList.addAll(downsampleVariantRegionNormally(mappings, refStart, coverageDistribution));   // todo -- maybe return the set to avoid going through the list every time?
                break;
            case Adaptive:
                readList.addAll(downsampleVariantRegionAdaptively(mappings, refStart, coverageDistribution));   // todo -- maybe return the set to avoid going through the list every time?
                break;
            default:
                throw new UserException.BadArgumentValue("dm" + downsampleStrategy.toString(), "Invalid value for downsample strategy");
        }

        return readList;
    }

    /**
     * Downsampling using the Normally Distributed strategy (this function is called by downsampleVariantRegion)
     *
     * This function will select reads at random from the given pool of reads starting from the middle loci
     * (likely where the variation is) guaranteeing the minimum coverage to be >= downsampleCoverage and
     * expands to the adjacent loci increasing their coverage as needed (by selecting more reads at random)
     * until all loci in the window are covered to the downsampleCoverage level.
     *
     * @param mappings             the two maps from read=>locus and loci=>read
     * @param refStart             the alignment position at the start of the array so we can use locusIndex appropriately
     * @param coverageDistribution the distribution of coverage before downsampling
     * @return all reads that pass the downsampling filtering process
     */
    protected Set<GATKSAMRecord> downsampleVariantRegionNormally(final Pair<HashMap<Integer, HashSet<GATKSAMRecord>>, HashMap<GATKSAMRecord, Boolean[]>> mappings, int refStart, int [] coverageDistribution ) {
        int [] downsampledCoverageDistribution = new int [coverageDistribution.length];              // the downsampled distribution array with only the selected reads
        HashSet<GATKSAMRecord> downsampledReads = new HashSet<GATKSAMRecord>();

        int middle = downsampledCoverageDistribution.length / 2;                                     // start covering with randomly selected reads from the middle of the window
        for (int i=middle; i < downsampledCoverageDistribution.length; i++) {                        // cover the entire window forward and backward
            downsampleLocus(mappings, i, refStart, downsampleCoverage, downsampledCoverageDistribution, downsampledReads);
            if (middle - i >= 0)
                downsampleLocus(mappings, middle - i, refStart, downsampleCoverage, downsampledCoverageDistribution, downsampledReads);
        }
        return downsampledReads;
    }

    protected Set<GATKSAMRecord> downsampleVariantRegionAdaptively(final Pair<HashMap<Integer, HashSet<GATKSAMRecord>>, HashMap<GATKSAMRecord, Boolean[]>> mappings, int refStart, int [] coverageDistribution ) {
        int [] downsampledCoverageDistribution = new int [coverageDistribution.length];
        HashSet<GATKSAMRecord> downsampledReads = new HashSet<GATKSAMRecord>();
        
        int minCoverage = MathUtils.arrayMin(coverageDistribution);    // find the base with the least coverage in the region
        int transform = Math.max(minCoverage-downsampleCoverage, 0);   // define our transformation to be subtraction that takes the minimum coverage point and brings it to the downsampleCoverage level

        for (int i=0; i < coverageDistribution.length; i++) {
            int goalCoverage = coverageDistribution[i] - transform;    // find our goal coverage for this locus and downsample it
            downsampleLocus(mappings, i, refStart, goalCoverage, downsampledCoverageDistribution, downsampledReads);
        }
        return downsampledReads;
    }

    /**
     * Internal function to downsample a given locus given all the necessary parameters.
     *
     * Note: This function WILL change (update) downsampledCoverageDistribution and downsampledReads accordingly.
     *
     * @param mappings                        the two maps from read=>locus and loci=>read
     * @param locusIndex                      the index of the loci in the array
     * @param refStart                        the alignment position at the start of the array so we can use locusIndex appropriately
     * @param coverageGoal                    the number of coverage we want to achieve for this site
     * @param downsampledCoverageDistribution the current coverage distribution with the selected reads so far -- will be updated in this function
     * @param downsampledReads                the list of reads selected by the downsampler so far -- will be updated in this function
     */
    protected void downsampleLocus(final Pair<HashMap<Integer, HashSet<GATKSAMRecord>>, HashMap<GATKSAMRecord, Boolean[]>> mappings, final int locusIndex, final int refStart, int coverageGoal, int[] downsampledCoverageDistribution, HashSet<GATKSAMRecord> downsampledReads) {
        HashMap<Integer, HashSet<GATKSAMRecord>> locusToReadMap = mappings.getFirst();                            // a map from every locus in the window to the reads that cover it
        HashMap<GATKSAMRecord, Boolean[]> readToLocusMap = mappings.getSecond();                                  // a map from every read in the window to the loci it covers
        HashSet readsOnThisLocus = locusToReadMap.get(refStart+locusIndex);                                       // Get the reads that are represented in this loci
        readsOnThisLocus.removeAll(downsampledReads);                                                             // Remove all reads that have already been chosen by previous loci
        int numberOfReadsToAdd = coverageGoal - downsampledCoverageDistribution[locusIndex];                      // we need to add this many reads to achieve the minimum downsampleCoverage
        if (numberOfReadsToAdd > 0) {                                                                             // no need to add reads if we're already covered or over covered.
            GATKSAMRecord [] readArray = convertArray(readsOnThisLocus.toArray());                                // convert to array so we can get a random subset
            GATKSAMRecord [] selectedReads = convertArray(MathUtils.randomSubset(readArray, numberOfReadsToAdd)); // get a random subset of the reads with the exact (or less if not available) number of reads we need to hit the downsample coverage
            downsampledReads.addAll(Arrays.asList(selectedReads));                                                // add the selected reads to the final set of reads

            for (GATKSAMRecord selectedRead : selectedReads) {                                                    // update the coverage distribution with the newly added reads
                Boolean [] lociAffected = readToLocusMap.get(selectedRead);                                       // get the boolean array to know which loci the read affected

                for (int j=0; j<downsampledCoverageDistribution.length; j++)
                    downsampledCoverageDistribution[j] += lociAffected[j] ? 1 : 0;                                // if it affects, increase coverage by one (no reduced reads in this pileup)
            }
        }
    }


    /**
     * Slides the window to the new position
     *
     * @param shift the new starting location of the window
     */
    protected void slide(int shift) {
        if (shift > 0) {
            LinkedList<GATKSAMRecord> newSlidingReads = new LinkedList<GATKSAMRecord>();

            // drop base baseCounts out of the window
            for (int i = 0; i < shift; i++)
                windowHeader.remove();
            startLocation += shift;

            // clip reads to new window
            int refLeftPosition = windowHeader.peekFirst().location;  // todo -- should I just get startLocation ?  (works if there are no gaps)
            for (GATKSAMRecord read : readsInWindow) {
                GATKSAMRecord clippedRead = clipStart(read, refLeftPosition);
                if (clippedRead != null)
                    newSlidingReads.add(clippedRead);
            }
            readsInWindow = newSlidingReads;
        }
    }

    /**
     * Properly closes a Sliding Window, finalizing all consensus and variant
     * regions that still exist regardless of being able to fulfill the
     * context size requirement in the end.
     *
     * @return All reads generated
     */
    public List<GATKSAMRecord> close() {
        // mark variant regions
        List<GATKSAMRecord> finalizedReads = new LinkedList<GATKSAMRecord>();

        if (!windowHeader.isEmpty()) {
            boolean[] variantSite = markSites(stopLocation + 1);

            // close everything (+1 to include the last site) -- consensus or variant region
            int sitesToClose = stopLocation - startLocation + 1;
            int start = 0;
            int i = 0;

            while (i < sitesToClose) {
                while (i < sitesToClose && !variantSite[i]) i++;
                finalizedReads.addAll(addToSyntheticReads(start, i));
                start = i;

                // close all variant regions regardless of having enough for context size
                // on the last one
                while (i < sitesToClose && variantSite[i]) i++;
                if (start <= i - 1)
                    finalizedReads.addAll(closeVariantRegion(start, i - 1));
                start = i;
            }

            // if it ended in running consensus, finish it up
            finalizedReads.addAll(finalizeAndAdd(ConsensusType.BOTH));

        }
        return finalizedReads;
    }

    /**
     * generates the SAM record for the running consensus read and resets it (to null)
     *
     * @return the read contained in the running consensus
     */
    protected GATKSAMRecord finalizeRunningConsensus() {
        GATKSAMRecord finalizedRead = null;
        if (runningConsensus != null) {
            if (runningConsensus.size() > 0)
                finalizedRead = runningConsensus.close();
            else
                consensusCounter--;

            runningConsensus = null;
        }
        return finalizedRead;
    }

    /**
     * generates the SAM record for the filtered data consensus and resets it (to null)
     *
     * @return the read contained in the running consensus
     */
    protected GATKSAMRecord finalizeFilteredDataConsensus() {
        GATKSAMRecord finalizedRead = null;
        if (filteredDataConsensus != null) {
            if (filteredDataConsensus.size() > 0)
                finalizedRead = filteredDataConsensus.close();
            else
                filteredDataConsensusCounter--;

            filteredDataConsensus = null;
        }
        return finalizedRead;
    }


    /**
     * Updates the sliding window's header counts with the incoming read bases, insertions
     * and deletions.
     *
     * @param read the incoming read to be added to the sliding window
     */
    @Requires("read.getAlignmentStart() >= startLocation")
    protected void updateHeaderCounts(GATKSAMRecord read) {
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();
        Cigar cigar = read.getCigar();

        int readBaseIndex = 0;
        int locationIndex = read.getAlignmentStart() - startLocation;

        // Do we need to add extra elements before the start of the header?
        // -- this may happen if the previous read was clipped and this alignment starts before the beginning of the window
        if (locationIndex < 0) {
            for (int i = 1; i <= -locationIndex; i++)
                windowHeader.addFirst(new HeaderElement(startLocation - i));

            // update start location accordingly
            startLocation = read.getAlignmentStart();
            locationIndex = 0;
        }

        // Do we need to add extra elements to the header?
        if (stopLocation < read.getAlignmentEnd()) {
            int elementsToAdd = (stopLocation < 0) ? read.getAlignmentEnd() - read.getAlignmentStart() + 1 : read.getAlignmentEnd() - stopLocation;
            while (elementsToAdd-- > 0)
                windowHeader.addLast(new HeaderElement(read.getAlignmentEnd() - elementsToAdd));

            // update stopLocation accordingly
            stopLocation = read.getAlignmentEnd();
        }

        // Special case for leading insertions before the beginning of the sliding read
        if (ReadUtils.readStartsWithInsertion(read).getFirst() && read.getAlignmentStart() == startLocation) {
            windowHeader.addFirst(new HeaderElement(startLocation - 1));          // create a new first element to the window header with no bases added
            startLocation--;                                                      // moving the start of the window back one position to the dummy head
            locationIndex = 1;                                                    // This allows the first element (I) to look at locationIndex - 1 in the subsequent switch and do the right thing.
        }

        Iterator<HeaderElement> headerElementIterator = windowHeader.listIterator(locationIndex);
        HeaderElement headerElement;
        for (CigarElement cigarElement : cigar.getCigarElements()) {
            switch (cigarElement.getOperator()) {
                case H:
                case S:                                                           // nothing to add to the window
                    break;
                case I:                                                           // insertions are added to the base to the left (previous element)
                    windowHeader.get(locationIndex - 1).addInsertionToTheRight();
                    readBaseIndex += cigarElement.getLength();
                    break;                                                        // just ignore the insertions at the beginning of the read
                case D:
                    int nDeletionsToAdd = cigarElement.getLength();
                    while (nDeletionsToAdd-- > 0) {                               // deletions are added to the baseCounts with the read mapping quality as it's quality score
                        headerElement = headerElementIterator.next();
                        headerElement.addBase((byte) 'D', (byte) read.getMappingQuality(), read.getMappingQuality());
                        locationIndex++;
                    }
                    break;
                case M:
                case P:
                case EQ:
                case X:
                    int nBasesToAdd = cigarElement.getLength();
                    while (nBasesToAdd-- > 0) {
                        headerElement = headerElementIterator.next();
                        headerElement.addBase(bases[readBaseIndex], quals[readBaseIndex], read.getMappingQuality());
                        readBaseIndex++;
                        locationIndex++;
                    }
                    break;
            }
        }
    }

    /**
     * Hard clips the beginning of the read if necessary. This function is used over and over
     * as the sliding window slides forward.
     *
     * @param read        the read to be clipped
     * @param refNewStart the new beginning of the sliding window
     * @return the original read if refNewStart < read.getAlignmentStart, null if refNewStart > read.getAlignmentEd, or the clipped read
     */
    protected GATKSAMRecord clipStart(GATKSAMRecord read, int refNewStart) {
        if (refNewStart > read.getAlignmentEnd())
            return null;
        if (refNewStart <= read.getAlignmentStart())
            return read;
        return ReadClipper.hardClipByReferenceCoordinatesLeftTail(read, refNewStart - 1);
    }

    /**
     * The element that describes the header of the sliding window.
     * <p/>
     * Each site has a header element containing the counts of each base, it's reference based location and whether or not the site
     * has insertions (to it's right). It also contains information about the bases that have been filtered out due to mapping or base
     * quality.
     */
    protected class HeaderElement {
        protected BaseCounts consensusBaseCounts;      // How many A,C,G,T (and D's) are in this site.
        protected BaseCounts filteredBaseCounts;       // How many A,C,G,T (and D's) were filtered out in this site.
        protected int insertionsToTheRight;            // How many reads in this site had insertions to the immediate right
        protected int location;                        // Genome location of this site (the sliding window knows which contig we're at
        protected LinkedList<Integer> mappingQuality;  // keeps the mapping quality of each read that contributed to this element (site)

        public int getLocation() {
            return location;
        }

        /**
         * Creates a new HeaderElement with the default location 0
         */
        public HeaderElement() {
            this(0);
        }

        /**
         * Creates a new HeaderElement with the following default values:
         * - empty consensusBaseCounts
         * - empty filteredBaseCounts
         * - 0 insertions to the right
         * - empty mappingQuality list
         *
         * @param location the reference location for the new element
         */
        public HeaderElement(int location) {
            this(new BaseCounts(), new BaseCounts(), 0, location, new LinkedList<Integer>());
        }

        /**
         * Creates a new HeaderElement with all given parameters
         *
         * @param consensusBaseCounts  the BaseCounts object for the running consensus synthetic read
         * @param filteredBaseCounts   the BaseCounts object for the filtered data synthetic read
         * @param insertionsToTheRight number of insertions to the right of this HeaderElement
         * @param location             the reference location of this reference element
         * @param mappingQuality       the list of mapping quality values of all reads that contributed to this HeaderElement
         */
        public HeaderElement(BaseCounts consensusBaseCounts, BaseCounts filteredBaseCounts, int insertionsToTheRight, int location, LinkedList<Integer> mappingQuality) {
            this.consensusBaseCounts = consensusBaseCounts;
            this.filteredBaseCounts = filteredBaseCounts;
            this.insertionsToTheRight = insertionsToTheRight;
            this.location = location;
            this.mappingQuality = mappingQuality;
        }

        /**
         * Whether or not the site represented by this HeaderElement is variant according to the
         * definitions of variant by insertion, deletion and mismatches.
         *
         * @return true if site is variant by any definition. False otherwise.
         */
        public boolean isVariant() {
            return hasConsensusData() && (isVariantFromInsertions() || isVariantFromMismatches() || isVariantFromDeletions());
        }

        /**
         * Adds a new base to the HeaderElement updating all counts accordingly
         *
         * @param base           the base to add
         * @param qual           the base quality
         * @param mappingQuality the mapping quality of the read this base belongs to
         */
        public void addBase(byte base, byte qual, int mappingQuality) {
            if (qual >= MIN_BASE_QUAL_TO_COUNT && mappingQuality >= MIN_MAPPING_QUALITY)
                consensusBaseCounts.incr(base, qual); // If the base passes filters, it is included in the consensus base counts
            else
                filteredBaseCounts.incr(base, qual);  // If the base fails filters, it is included with the filtered data base counts

            this.mappingQuality.add(mappingQuality);  // Filtered or not, the RMS mapping quality includes all bases in this site
        }

        /**
         * Adds an insertions to the right of the HeaderElement and updates all counts accordingly.
         * All insertions should be added to the right of the element.
         */
        public void addInsertionToTheRight() {
            insertionsToTheRight++;
        }

        /**
         * Does this HeaderElement contain consensus data?
         *
         * @return whether or not this HeaderElement contains consensus data
         */
        public boolean hasConsensusData() {
            return consensusBaseCounts.totalCount() > 0;
        }

        /**
         * Does this HeaderElement contain filtered data?
         *
         * @return whether or not this HeaderElement contains filtered data
         */
        public boolean hasFilteredData() {
            return filteredBaseCounts.totalCount() > 0;
        }

        /**
         * A HeaderElement is empty if it has no consensus or filtered data
         *
         * @return whether or not this HeaderElement has no data
         */
        public boolean isEmpty() {
            return (!hasFilteredData() && !hasConsensusData());
        }

        /**
         * The RMS of the mapping qualities of all reads that contributed to this HeaderElement
         *
         * @return the RMS of the mapping qualities of all reads that contributed to this HeaderElement
         */
        public double getRMS() {
            return MathUtils.rms(mappingQuality);
        }

        /**
         * Whether or not the HeaderElement is variant due to excess insertions
         *
         * @return whether or not the HeaderElement is variant due to excess insertions
         */
        private boolean isVariantFromInsertions() {
            int numberOfBases = consensusBaseCounts.totalCount();
            if (numberOfBases == 0 && insertionsToTheRight > 0)
                return true;  // we only have insertions
            else if (numberOfBases == 0)
                return false; // we don't have anything

            // if we have bases and insertions, check the ratio
            return ((double) insertionsToTheRight / numberOfBases) > MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT;
        }

        /**
         * Whether or not the HeaderElement is variant due to excess deletions
         *
         * @return whether or not the HeaderElement is variant due to excess insertions
         */
        private boolean isVariantFromDeletions() {
            return consensusBaseCounts.baseIndexWithMostCounts() == BaseIndex.D || consensusBaseCounts.baseCountProportion(BaseIndex.D) > MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT;
        }

        /**
         * Whether or not the HeaderElement is variant due to excess mismatches
         *
         * @return whether or not the HeaderElement is variant due to excess insertions
         */
        private boolean isVariantFromMismatches() {
            BaseIndex mostCommon = consensusBaseCounts.baseIndexWithMostCountsWithoutIndels();
            double mostCommonProportion = consensusBaseCounts.baseCountProportionWithoutIndels(mostCommon);
            return mostCommonProportion != 0.0 && mostCommonProportion < (1 - MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT);
        }

    }

    /**
     * The types of synthetic to use in the finalizeAndAdd method
     */
    private enum ConsensusType {
        CONSENSUS,
        FILTERED,
        BOTH
    }

    /**
     * workaround function to convert an object array to a GATKSAMRecord array since we
     * can't typecast collections
     *
     * @param array a GATKSAMRecord array that has been typecasted to Object array
     * @return a new array with the typecasts to GATKSAMRecord
     */
    private static GATKSAMRecord [] convertArray (Object [] array ) {
        GATKSAMRecord [] result = new GATKSAMRecord [array.length];
        for (int i=0; i<array.length; i++)
            result[i] = (GATKSAMRecord) array[i];
        return result;
    }

}
