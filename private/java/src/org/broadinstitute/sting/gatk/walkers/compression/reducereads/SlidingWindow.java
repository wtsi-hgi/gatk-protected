package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.MathUtils;
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
    final private LinkedList<GATKSAMRecord> readsInWindow;
    final private LinkedList<HeaderElement> windowHeader;
    protected int contextSize;                                                                                          // the largest context size (between mismatches and indels)
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
    protected double MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT;                                                        // proportion has to be greater than this value to trigger variant region due to mismatches
    protected double MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT;                                                      // proportion has to be greater than this value to trigger variant region due to deletions
    protected int MIN_BASE_QUAL_TO_COUNT;                                                                               // qual has to be greater than or equal to this value
    protected int MIN_MAPPING_QUALITY;

    protected ReduceReadsWalker.DownsampleStrategy downsampleStrategy;

    /**
     * The types of synthetic reads to use in the finalizeAndAdd method
     */
    private enum ConsensusType {
        CONSENSUS,
        FILTERED,
        BOTH
    }

    public int getStopLocation() {
        return stopLocation;
    }

    public String getContig() {
        return contig;
    }

    public int getContigIndex() {
        return contigIndex;
    }

    public int getStartLocation() {
        return windowHeader.isEmpty() ? -1 : windowHeader.peek().getLocation();
    }


    public SlidingWindow(String contig, int contigIndex, int contextSize, SAMFileHeader header, GATKSAMReadGroupRecord readGroupAttribute, int windowNumber, final double minAltProportionToTriggerVariant, final double minIndelProportionToTriggerVariant, int minBaseQual, int minMappingQuality, int downsampleCoverage, final ReduceReadsWalker.DownsampleStrategy downsampleStrategy) {
        this.stopLocation = -1;
        this.contextSize = contextSize;
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
     * 
     * Reads are assumed to be in order, therefore, when a read is added the sliding window can
     * assume that no more reads will affect read.getUnclippedStart() - contextSizeMismatches. The window
     * slides forward to that position and returns all reads that may have been finalized in the
     * sliding process.
     *
     * @param read the read
     * @return a list of reads that have been finished by sliding the window.
     */
    public List<GATKSAMRecord> addRead(GATKSAMRecord read) {
        updateHeaderCounts(read, false);                                                                                // update the window header counts
        readsInWindow.add(read);                                                                                        // add read to sliding reads
        return slideWindow(read.getUnclippedStart());
    }

    /**
     * returns the next complete or incomplete variant region between 'from' (inclusive) and 'to' (exclusive)
     *
     * @param from         beginning window header index of the search window (inclusive)
     * @param to           end window header index of the search window (exclusive)
     * @param variantSite  boolean array with true marking variant regions
     * @return null if nothing is variant, start/stop if there is a complete variant region, start/-1 if there is an incomplete variant region.
     */
    private Pair<Integer, Integer> getNextVariantRegion(int from, int to, boolean[] variantSite) {
        boolean foundStart = false;
        int variantRegionStartIndex = 0;
        for (int i=from; i<to; i++) {
            if (variantSite[i] && !foundStart) {
                variantRegionStartIndex = i;
                foundStart = true;
            }
            else if(!variantSite[i] && foundStart) {
                return(new Pair<Integer, Integer>(variantRegionStartIndex, i-1));
            }
        }
        return (foundStart) ? new Pair<Integer, Integer>(variantRegionStartIndex, -1) : null;
    }

    /**
     * Creates a list with all the complete and incomplete variant regions within 'from' (inclusive) and 'to' (exclusive)
     *
     * @param from         beginning window header index of the search window (inclusive)
     * @param to           end window header index of the search window (exclusive)
     * @param variantSite  boolean array with true marking variant regions
     * @return a list with start/stops of variant regions following getNextVariantRegion description
     */
    private List<Pair<Integer, Integer>> getAllVariantRegions(int from, int to, boolean[] variantSite) {
        List<Pair<Integer,Integer>> regions = new LinkedList<Pair<Integer, Integer>>();
        int index = from;
        while(index < to) {
            Pair<Integer,Integer> result = getNextVariantRegion(index, to, variantSite);
            if (result == null)
                break;

            regions.add(result);
            if (result.getSecond() < 0)
                break;
            index = result.getSecond() + 1;
        }
        return regions;
    }


    /**
     * Determines if the window can be slid given the new incoming read.
     *
     * We check from the start of the window to the (unclipped) start of the new incoming read if there
     * is any variant.
     * If there are variant sites, we check if it's time to close the variant region.
     *
     * @param incomingReadUnclippedStart the incoming read's start position. Must be the unclipped start!
     * @return all reads that have fallen to the left of the sliding window after the slide
     */
    protected List<GATKSAMRecord> slideWindow(int incomingReadUnclippedStart) {
        List<GATKSAMRecord> finalizedReads = new LinkedList<GATKSAMRecord>();

        if (incomingReadUnclippedStart - contextSize > getStartLocation()) {
            int readStartHeaderIndex = incomingReadUnclippedStart - getStartLocation();
            boolean[] variantSite = markSites(getStartLocation() + readStartHeaderIndex);
            int breakpoint = Math.max(readStartHeaderIndex - contextSize - 1, 0);                                       // this is the limit of what we can close/send to consensus (non-inclusive)

            List<Pair<Integer,Integer>> regions = getAllVariantRegions(0, breakpoint, variantSite);
            finalizedReads = closeVariantRegions(regions, false);

            List<GATKSAMRecord> readsToRemove = new LinkedList<GATKSAMRecord>();
            for (GATKSAMRecord read : readsInWindow) {                                                                  // todo -- unnecessarily going through all reads in the window !! Optimize this (But remember reads are not sorted by alignment end!)
                if (read.getAlignmentEnd() < getStartLocation()) {
                    readsToRemove.add(read);
                }
            }
            for (GATKSAMRecord read : readsToRemove) {
                readsInWindow.remove(read);
            }
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

        boolean[] markedSites = new boolean[stop - getStartLocation() + contextSize + 1];

        Iterator<HeaderElement> headerElementIterator = windowHeader.iterator();
        for (int i = getStartLocation(); i < stop; i++) {
            if (headerElementIterator.hasNext()) {
                HeaderElement headerElement = headerElementIterator.next();

                if (headerElement.isVariantFromDeletions(MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT) || headerElement.isVariantFromInsertions(MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT) || headerElement.isVariantFromMismatches(MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT))
                    markVariantRegion(markedSites, i - getStartLocation());

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
     */
    protected void markVariantRegion(boolean[] markedSites, int variantSiteLocation) {
        int from = (variantSiteLocation < contextSize) ? 0 : variantSiteLocation - contextSize;
        int to = (variantSiteLocation + contextSize + 1 > markedSites.length) ? markedSites.length : variantSiteLocation + contextSize + 1;
        for (int i = from; i < to; i++)
            markedSites[i] = true;
    }

    /**
     * Adds bases to the running consensus or filtered data accordingly
     * 
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
        if (read != null)
            list.add(read);

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
     * 
     * Different from the addToConsensus method, this method assumes a contiguous sequence of filteredData
     * bases.
     *
     * @param start the first header index to add to consensus
     * @param end   the first header index NOT TO add to consensus
     */
    private void addToFilteredData(int start, int end) {
        if (filteredDataConsensus == null)
            filteredDataConsensus = new SyntheticRead(header, readGroupAttribute, contig, contigIndex, filteredDataReadName + filteredDataConsensusCounter++, windowHeader.get(start).getLocation(), GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG);

        ListIterator<HeaderElement> headerElementIterator = windowHeader.listIterator(start);
        for (int index = start; index < end; index++) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("Requested to create a filtered data synthetic read from " + start + " to " + end + " but " + index + " does not exist");

            HeaderElement headerElement = headerElementIterator.next();
            if (headerElement.hasConsensusData())
                throw new ReviewedStingException("Found consensus data inside region to add to filtered data.");

            if (!headerElement.hasFilteredData())
                throw new ReviewedStingException("No filtered data in " + index);

            BaseIndex base = headerElement.getFilteredBaseCounts().baseIndexWithMostCounts();
            byte count = (byte) Math.min(headerElement.getFilteredBaseCounts().countOfMostCommonBase(), Byte.MAX_VALUE);
            byte qual = headerElement.getFilteredBaseCounts().averageQualsOfMostCommonBase();
            filteredDataConsensus.add(base, count, qual, headerElement.getRMS());
        }
    }

    /**
     * Adds bases to the filtered data synthetic read.
     * 
     * Different from the addToConsensus method, this method assumes a contiguous sequence of filteredData
     * bases.
     *
     * @param start the first header index to add to consensus
     * @param end   the first header index NOT TO add to consensus
     */
    private void addToRunningConsensus(int start, int end) {
        if (runningConsensus == null)
            runningConsensus = new SyntheticRead(header, readGroupAttribute, contig, contigIndex, consensusReadName + consensusCounter++, windowHeader.get(start).getLocation(), GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG);

        Iterator<HeaderElement> headerElementIterator = windowHeader.listIterator(start);
        for (int index = start; index < end; index++) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("Requested to create a running consensus synthetic read from " + start + " to " + end + " but " + index + " does not exist");

            HeaderElement headerElement = headerElementIterator.next();
            if (!headerElement.hasConsensusData())
                throw new ReviewedStingException("No CONSENSUS data in " + index);

            BaseIndex base = headerElement.getConsensusBaseCounts().baseIndexWithMostCounts();
            byte count = (byte) Math.min(headerElement.getConsensusBaseCounts().countOfMostCommonBase(), Byte.MAX_VALUE);
            byte qual = headerElement.getConsensusBaseCounts().averageQualsOfMostCommonBase();
            runningConsensus.add(base, count, qual, headerElement.getRMS());
        }
    }


    /**
     * Finalizes a variant region, any adjacent synthetic reads.
     *
     * @param start the first window header index in the variant region (inclusive)
     * @param stop  the last window header index of the variant region (inclusive)
     * @return all reads contained in the variant region plus any adjacent synthetic reads
     */
    @Requires("start <= stop")
    protected List<GATKSAMRecord> closeVariantRegion(int start, int stop) {
        List<GATKSAMRecord> allReads = new LinkedList<GATKSAMRecord>();

        int refStart = windowHeader.get(start).getLocation();                                                           // All operations are reference based, not read based
        int refStop = windowHeader.get(stop).getLocation();

        for (GATKSAMRecord read : readsInWindow) {                                                                      // Keep all reads that overlap the variant region
            if (read.getAlignmentStart() <= refStop && read.getAlignmentEnd() >= refStart) {
                allReads.add(read);
                updateHeaderCounts(read, true);                                                                         // Remove this read from the window header entirely
            }
        }

        List<GATKSAMRecord> result = (downsampleCoverage > 0) ? downsampleVariantRegion(allReads, refStart, refStop) : allReads;
        result.addAll(addToSyntheticReads(0, start));
        result.addAll(finalizeAndAdd(ConsensusType.BOTH));

        for (GATKSAMRecord read : result) {
            readsInWindow.remove(read);                                                                                 // todo -- not optimal, but needs to be done so the next region doesn't try to remove the same reads from the header counts.
        }

        return result;                                                                                                  // finalized reads will be downsampled if necessary
    }


    private List<GATKSAMRecord> closeVariantRegions(List<Pair<Integer, Integer>> regions, boolean forceClose) {
        List<GATKSAMRecord> allReads = new LinkedList<GATKSAMRecord>();
        if (!regions.isEmpty()) {
            int lastStop = -1;
            for (Pair<Integer, Integer> region : regions) {
                int start = region.getFirst();
                int stop = region.getSecond();
                if (stop < 0 && forceClose)
                    stop = windowHeader.size() - 1;
                if (stop >= 0) {
                    allReads.addAll(closeVariantRegion(start, stop));
                    lastStop = stop;
                }
            }
            for (int i = 0; i <= lastStop; i++)                                                                         // clean up the window header elements up until the end of the variant region.
                windowHeader.remove();                                                                                  // todo -- can't believe java doesn't allow me to just do windowHeader = windowHeader.get(stop). Should be more efficient here!
        }
        return allReads;
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
        int [] coverageDistribution = ReadUtils.getCoverageDistributionOfReads(allReads, refStart, refStop);            // the full coverage distribution array for the variant region with all the reads
        switch (downsampleStrategy) {
            case Normal:
                readList.addAll(downsampleVariantRegionNormally(mappings, refStart, coverageDistribution));             // todo -- maybe return the set to avoid going through the list every time?
                break;
            case Adaptive:
                readList.addAll(downsampleVariantRegionAdaptively(mappings, refStart, coverageDistribution));           // todo -- maybe return the set to avoid going through the list every time?
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
        int [] downsampledCoverageDistribution = new int [coverageDistribution.length];                                 // the downsampled distribution array with only the selected reads
        HashSet<GATKSAMRecord> downsampledReads = new HashSet<GATKSAMRecord>();

        int middle = downsampledCoverageDistribution.length / 2;                                                        // start covering with randomly selected reads from the middle of the window
        for (int i=middle; i < downsampledCoverageDistribution.length; i++) {                                           // cover the entire window forward and backward
            downsampleLocus(mappings, i, refStart, downsampleCoverage, downsampledCoverageDistribution, downsampledReads);
            if (middle - i >= 0)
                downsampleLocus(mappings, middle - i, refStart, downsampleCoverage, downsampledCoverageDistribution, downsampledReads);
        }
        return downsampledReads;
    }

    protected Set<GATKSAMRecord> downsampleVariantRegionAdaptively(final Pair<HashMap<Integer, HashSet<GATKSAMRecord>>, HashMap<GATKSAMRecord, Boolean[]>> mappings, int refStart, int [] coverageDistribution ) {
        int [] downsampledCoverageDistribution = new int [coverageDistribution.length];
        HashSet<GATKSAMRecord> downsampledReads = new HashSet<GATKSAMRecord>();
        
        int minCoverage = MathUtils.arrayMin(coverageDistribution);                                                     // find the base with the least coverage in the region
        int transform = Math.max(minCoverage-downsampleCoverage, 0);                                                    // define our transformation to be subtraction that takes the minimum coverage point and brings it to the downsampleCoverage level

        for (int i=0; i < coverageDistribution.length; i++) {
            int goalCoverage = coverageDistribution[i] - transform;                                                     // find our goal coverage for this locus and downsample it
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
        HashMap<Integer, HashSet<GATKSAMRecord>> locusToReadMap = mappings.getFirst();                                  // a map from every locus in the window to the reads that cover it
        HashMap<GATKSAMRecord, Boolean[]> readToLocusMap = mappings.getSecond();                                        // a map from every read in the window to the loci it covers
        HashSet readsOnThisLocus = locusToReadMap.get(refStart + locusIndex);                                           // Get the reads that are represented in this loci
        readsOnThisLocus.removeAll(downsampledReads);                                                                   // Remove all reads that have already been chosen by previous loci
        int numberOfReadsToAdd = coverageGoal - downsampledCoverageDistribution[locusIndex];                            // we need to add this many reads to achieve the minimum downsampleCoverage
        if (numberOfReadsToAdd > 0) {                                                                                   // no need to add reads if we're already covered or over covered.
            GATKSAMRecord [] readArray = convertArray(readsOnThisLocus.toArray());                                      // convert to array so we can get a random subset
            GATKSAMRecord [] selectedReads = convertArray(MathUtils.randomSubset(readArray, numberOfReadsToAdd));       // get a random subset of the reads with the exact (or less if not available) number of reads we need to hit the downsample coverage
            downsampledReads.addAll(Arrays.asList(selectedReads));                                                      // add the selected reads to the final set of reads

            for (GATKSAMRecord selectedRead : selectedReads) {                                                          // update the coverage distribution with the newly added reads
                Boolean [] lociAffected = readToLocusMap.get(selectedRead);                                             // get the boolean array to know which loci the read affected

                for (int j=0; j<downsampledCoverageDistribution.length; j++)
                    downsampledCoverageDistribution[j] += lociAffected[j] ? 1 : 0;                                      // if it affects, increase coverage by one (no reduced reads in this pileup)
            }
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
            List<Pair<Integer,Integer>> regions = getAllVariantRegions(0, windowHeader.size(), variantSite);
            finalizedReads = closeVariantRegions(regions, true);

            if (!windowHeader.isEmpty()) {
                finalizedReads.addAll(addToSyntheticReads(0, windowHeader.size() - 1));
                finalizedReads.addAll(finalizeAndAdd(ConsensusType.BOTH));                                              // if it ended in running consensus, finish it up
            }

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
    protected void updateHeaderCounts(GATKSAMRecord read, boolean removeRead) {
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();
        Cigar cigar = read.getCigar();

        int readBaseIndex = 0;
        int startLocation = getStartLocation();
        int locationIndex = startLocation < 0 ? 0 : read.getAlignmentStart() - startLocation;

        if (removeRead && locationIndex < 0)
            throw new ReviewedStingException("read is behind the Sliding Window. read: " + read + " cigar: " + read.getCigarString() + " window: " + startLocation + "," + stopLocation);

        if (!removeRead) {                                                                                              // we only need to create new header elements if we are adding the read, not when we're removing it
            if (locationIndex < 0) {                                                                                    // Do we need to add extra elements before the start of the header? -- this may happen if the previous read was clipped and this alignment starts before the beginning of the window
                for (int i = 1; i <= -locationIndex; i++)
                    windowHeader.addFirst(new HeaderElement(startLocation - i));

                startLocation = read.getAlignmentStart();                                                               // update start location accordingly
                locationIndex = 0;
            }

            if (stopLocation < read.getAlignmentEnd()) {                                                                // Do we need to add extra elements to the header?
                int elementsToAdd = (stopLocation < 0) ? read.getAlignmentEnd() - read.getAlignmentStart() + 1 : read.getAlignmentEnd() - stopLocation;
                while (elementsToAdd-- > 0)
                    windowHeader.addLast(new HeaderElement(read.getAlignmentEnd() - elementsToAdd));

                stopLocation = read.getAlignmentEnd();                                                                  // update stopLocation accordingly
            }

            // Special case for leading insertions before the beginning of the sliding read
            if (ReadUtils.readStartsWithInsertion(read).getFirst() && (read.getAlignmentStart() == startLocation || startLocation < 0)) {
                windowHeader.addFirst(new HeaderElement(read.getAlignmentStart() - 1));                                 // create a new first element to the window header with no bases added
                locationIndex = 1;                                                                                      // This allows the first element (I) to look at locationIndex - 1 in the subsequent switch and do the right thing.
            }
        }

        Iterator<HeaderElement> headerElementIterator = windowHeader.listIterator(locationIndex);
        HeaderElement headerElement;
        for (CigarElement cigarElement : cigar.getCigarElements()) {
            switch (cigarElement.getOperator()) {
                case H:
                case S:                                                                                                 // nothing to add to the window
                    break;
                case I:

                    if (removeRead && locationIndex == 0) {                                                             // special case, if we are removing a read that starts in insertion and we don't have the previous header element anymore, don't worry about it.
                        break;
                    }

                    headerElement = windowHeader.get(locationIndex - 1);                                                // insertions are added to the base to the left (previous element)
                    if (removeRead) {
                        headerElement.removeInsertionToTheRight();
                    }
                    else
                        headerElement.addInsertionToTheRight();
                    readBaseIndex += cigarElement.getLength();
                    break;                                                                                              // just ignore the insertions at the beginning of the read
                case D:
                    int nDeletions = cigarElement.getLength();
                    while (nDeletions-- > 0) {                                                                          // deletions are added to the baseCounts with the read mapping quality as it's quality score
                        headerElement = headerElementIterator.next();
                        if (removeRead)
                            headerElement.removeBase((byte) 'D', (byte) read.getMappingQuality(), read.getMappingQuality(), MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY);
                        else
                            headerElement.addBase((byte) 'D', (byte) read.getMappingQuality(), read.getMappingQuality(), MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY);

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
                        if (removeRead)
                            headerElement.removeBase(bases[readBaseIndex], quals[readBaseIndex], read.getMappingQuality(), MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY);
                        else
                            headerElement.addBase(bases[readBaseIndex], quals[readBaseIndex], read.getMappingQuality(), MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY);

                        readBaseIndex++;
                        locationIndex++;
                    }
                    break;
            }
        }
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

