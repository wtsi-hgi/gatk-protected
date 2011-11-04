package org.broadinstitute.sting.gatk.walkers.reducereads;

import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 8/3/11
 * Time: 2:24 PM
 */
public class SlidingWindow {

    protected final static String RG_POSTFIX = ".ReducedReads";

    // Sliding Window data
    protected LinkedList<SlidingRead> SlidingReads;
    protected LinkedList<HeaderElement> windowHeader;
    protected int contextSize;
    protected int contextSizeIndels;
    protected int startLocation;
    protected int stopLocation;

    // Running consensus data
    protected RunningConsensus runningConsensus;
    protected int consensusCounter;
    protected String contig;
    protected int contigIndex;
    protected SAMFileHeader header;
    protected Object readGroupAttribute;
    protected String readName;

    // Additional parameters
    protected double MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT;    // proportion has to be greater than this value to trigger variant region due to mismatches
    protected double MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT;  // proportion has to be greater than this value to trigger variant region due to deletions
    protected int MIN_BASE_QUAL_TO_COUNT;                           // qual has to be greater than or equal to this value
    protected int MAX_QUAL_COUNT;                                   // to avoid blowing up the qual field of a consensus site
    protected int MIN_MAPPING_QUALITY;


    public SlidingWindow(String contig, int contigIndex, int contextSize, int contextSizeIndels, SAMFileHeader header, Object readGroupAttribute, int windowNumber, final double minAltProportionToTriggerVariant, final double minIndelProportionToTriggerVariant, int minBaseQual, int maxQualCount, int minMappingQuality) {
        this.startLocation = -1;
        this.stopLocation = -1;
        this.contextSize = contextSize;
        this.contextSizeIndels = contextSizeIndels;

        this.MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT = minAltProportionToTriggerVariant;
        this.MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT = minIndelProportionToTriggerVariant;
        this.MAX_QUAL_COUNT = maxQualCount;
        this.MIN_BASE_QUAL_TO_COUNT = minBaseQual;
        this.MIN_MAPPING_QUALITY = minMappingQuality;

        this.windowHeader = new LinkedList<HeaderElement>();
        this.SlidingReads = new LinkedList<SlidingRead>();

        this.consensusCounter = 0;

        this.contig = contig;
        this.contigIndex = contigIndex;
        this.header = header;
        this.readGroupAttribute = readGroupAttribute;
        this.readName = "Consensus-" + windowNumber + "-";

        this.runningConsensus = null;
    }

    public int getStartLocation() {
        return startLocation;
    }

    public int getStopLocation() {
        return stopLocation;
    }

    public String getContig() {
        return contig;
    }

    public List<GATKSAMRecord> addRead( GATKSAMRecord read ) {
        // If this is the first read in the window, update startLocation
        if (startLocation < 0)
            startLocation = read.getAlignmentStart();

        List<GATKSAMRecord> finalizedReads = slideIfPossible(read.getUnclippedStart());

        // update the window header counts
        updateHeaderCounts(read);

        // add read to sliding reads
        SlidingReads.add(new SlidingRead(read));

        return finalizedReads;
    }

    public int getContigIndex() {
        return contigIndex;
    }

    /**
     * Determines if the window can be slid given the new incoming read.
     *
     * We check from the start of the window to the (unclipped) start of the new incoming read if there
     * is any variant.
     *
     * If there are no variant sites, we slide the left side of the window to the unclipped start of the
     * new incoming read minus the context size (we can only guarantee that the reads before this point
     * will not be in any variant region ever)
     *
     * If there are variant sites, we check if it's time to close the variant region.
     *
     * @param incomingReadUnclippedStart the incoming read's start position. Must be the unclipped start!
     * @return all reads that have fallen to the left of the sliding window after the slide
     */
    protected List<GATKSAMRecord> slideIfPossible (int incomingReadUnclippedStart) {

        List<GATKSAMRecord> finalizedReads = new LinkedList<GATKSAMRecord>();

        // No point doing any calculation if the read's unclipped start is behind
        // the start of the window
        if (incomingReadUnclippedStart - contextSize > startLocation) {
            int positionShift = incomingReadUnclippedStart - startLocation;
            boolean [] variantSite = markSites(startLocation + positionShift);

            // this is the limit of what we can close/send to consensus (non-inclusive)
            int breakpoint = Math.max(positionShift - contextSize, 0);
            int start = 0;
            int i = 0;

            while (i < breakpoint) {
                while (i<breakpoint && !variantSite[i]) i++;
                finalizedReads.addAll(addToConsensus(start, i));

                start = i;
                while(i<breakpoint && variantSite[i]) i++;
                if ( i < breakpoint || (i-1 > 0 && variantSite[i-1] && !variantSite[i]) ) {
                    finalizedReads.addAll(closeVariantRegion(start, i-1));
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
     * @param stop check the window from start to 'upto'
     * @return a boolean array with 'true' marking variant regions and false marking consensus sites
     */
    protected boolean [] markSites(int stop) {

        boolean [] markedSites = new boolean [stop - startLocation + contextSize + 1];

        Iterator<HeaderElement> headerElementIterator = windowHeader.iterator();
        for (int i = startLocation; i < stop; i++) {
            if (headerElementIterator.hasNext()) {
                HeaderElement headerElement = headerElementIterator.next();
                if (headerElement.isVariant())
                    markVariantRegion(markedSites, i - startLocation);
            }
            else
                break;
        }
        return markedSites;
    }

    /**
     * performs the variant region marks (true) around a variant site
     *
     * @param markedSites the boolean array to bear the marks
     * @param variantSiteLocation the location where a variant site was found
     */
    protected void markVariantRegion (boolean [] markedSites, int variantSiteLocation) {
        int from = (variantSiteLocation < contextSize) ? 0 : variantSiteLocation - contextSize;
        int to = (variantSiteLocation + contextSize + 1 > markedSites.length) ? markedSites.length : variantSiteLocation + contextSize + 1;
        for (int i=from; i<to; i++)
            markedSites[i] = true;
    }

    /**
     * If adding a sequence with gaps, it will finalize multiple consensus reads and keep the last running consensus
     *
     * @param start the first header index to add to consensus
     * @param end the first header index NOT TO add to consensus
     * @return a list of consensus reads generated by this call. Empty list if no consensus was generated.
     */
    protected List<GATKSAMRecord> addToConsensus (int start, int end) {
        LinkedList<GATKSAMRecord> consensusList = new LinkedList<GATKSAMRecord>();
        if (start < end) {
            if (runningConsensus == null)
                runningConsensus = new RunningConsensus(header, readGroupAttribute, contig, contigIndex, readName + consensusCounter++, windowHeader.get(start).location, MIN_BASE_QUAL_TO_COUNT);

            int i = 0;
            for (HeaderElement wh : windowHeader) {
                if (i == end)
                    break;
                if (i >= start) {
                    // Element may be an empty element representing a gap between the reads
                    if (wh.isEmpty()) {
                        GATKSAMRecord consensus = finalizeConsensus();        // we finalize the running consensus and start a new one with the remaining
                        if(consensus != null)
                            consensusList.add(consensus);
                        consensusList.addAll(addToConsensus(i + 1, end)); // and start a new one starting at the next position
                        break;                                            // recursive call takes care of the rest of this loop, we are done.
                    }
                    else {
                        byte base  = wh.baseCounts.baseWithMostCounts();
                        byte count = (byte) Math.min(wh.baseCounts.countOfMostCommonBase(), MAX_QUAL_COUNT);
                        byte qual  = wh.baseCounts.averageQualsOfMostCommonBase();
                        runningConsensus.add(base, count, qual, wh.getRMS());
                    }
                }
                i++;
            }
        }
        return consensusList;
    }

    @Requires("start <= end")
    protected List<GATKSAMRecord> closeVariantRegion(int start, int end) {
        List<GATKSAMRecord> finalizedReads = new LinkedList<GATKSAMRecord>();

        // Finalize consensus if there is one to finalize before the Variant Region
        if (runningConsensus != null) {
            GATKSAMRecord consensus = finalizeConsensus();
            if (consensus != null)
                finalizedReads.add(consensus);
        }
        // Clipping operations are reference based, not read based
        int refStart = windowHeader.get(start).location;
        int refEnd = windowHeader.get(end).location;         // update to refStart + end?

        for ( SlidingRead read: SlidingReads ) {
            GATKSAMRecord SAM = read.trimToVariableRegion(refStart, refEnd);
            //System.out.println("HardClippedEnds:  (" + refStart +","+ refStop +") " + SAM.getCigarString() + "\t" + SAM.getAlignmentStart() + "\t" + SAM.getAlignmentEnd());
            if ( SAM.getReadLength() > 0 ) {
                finalizedReads.add(SAM);
            }
        }
        return finalizedReads;
    }

    /**
     * Slides the window to the new position
     *
     * @param shift the new starting location of the window
     */
    protected void slide (int shift) {
        if (shift > 0) {
            LinkedList<SlidingRead> newSlidingReads = new LinkedList<SlidingRead>();

            // drop base baseCounts out of the window
            for (int i=0; i<shift; i++)
                windowHeader.remove();
            startLocation += shift;

            // clip reads to new window
            int refLeftPosition = windowHeader.peekFirst().location;  // should I just get startLocation ?  (works if there are no gaps)
            for ( SlidingRead slidingRead: SlidingReads ) {
                SlidingRead sr = slidingRead.clipStart(refLeftPosition);
                if (sr != null)
                    newSlidingReads.add(sr);
            }
            SlidingReads = newSlidingReads;
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
            boolean [] variantSite = markSites(stopLocation+1);

            // close everything (+1 to include the last site) -- consensus or variant region
            int sitesToClose = stopLocation - startLocation + 1;
            int start = 0;
            int i = 0;

            while (i < sitesToClose) {
                while (i<sitesToClose && !variantSite[i]) i++;
                finalizedReads.addAll(addToConsensus(start, i));
                start = i;

                // close all variant regions regardless of having enough for context size
                // on the last one
                while(i<sitesToClose && variantSite[i]) i++;
                if (start <= i-1)
                    finalizedReads.addAll(closeVariantRegion(start, i-1));
                start = i;
            }
            // if it ended in consensus, finish it up
            if (runningConsensus != null) {
                GATKSAMRecord consensus = finalizeConsensus();
                if (consensus != null)
                    finalizedReads.add(consensus);
            }
        }
        return finalizedReads;
    }

    /**
     * generates the SAM record for the consensus read and resets the runningConsensus
     * object (to null).
     *
     * @return the read contained in the running consensus
     */
    @Requires("runningConsensus != null")
    protected GATKSAMRecord finalizeConsensus() {
        GATKSAMRecord finalizedRead = null;
        if (runningConsensus.size() > 0)
            finalizedRead = runningConsensus.close();
        else
            consensusCounter--;

        runningConsensus = null;
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
            for(int i = 1; i <= -locationIndex; i++)
                windowHeader.addFirst(newHeaderElement(startLocation - i));

            // update start location accordingly
            startLocation = read.getAlignmentStart();
            locationIndex = 0;
        }

        // Do we need to add extra elements to the header?
        if (stopLocation < read.getAlignmentEnd()) {
            int elementsToAdd = (stopLocation < 0) ? read.getAlignmentEnd() - read.getAlignmentStart() + 1 : read.getAlignmentEnd() - stopLocation;
            while (elementsToAdd-- > 0)
                windowHeader.addLast(newHeaderElement(read.getAlignmentEnd() - elementsToAdd));

            // update stopLocation accordingly
            stopLocation = read.getAlignmentEnd();
        }

        // Reads that don't pass the minimum mapping quality filter are not added to the
        // consensus, or count towards a variant region so no point in keeping track of
        // their base counts.
        if (read.getMappingQuality() < MIN_MAPPING_QUALITY) {
            Iterator<HeaderElement> headerElementIterator = windowHeader.listIterator(locationIndex);
            for (int i = read.getAlignmentStart(); i <= read.getAlignmentEnd(); i++) {
                if (!headerElementIterator.hasNext())
                    throw new ReviewedStingException("No header element for position: + " + read.getReferenceName() + ":" + i);

                HeaderElement headerElement = headerElementIterator.next();
                headerElement.addMappingQuality(read.getMappingQuality());
            }
        }
        else {
            // todo -- perhaps rewrite this iteration using list iterator to save time searching for each index.
            // todo -- they should be consecutive (as far as I can tell)
            for (CigarElement cigarElement : cigar.getCigarElements()) {
                switch (cigarElement.getOperator()) {
                    case H:
                    case S:
                        // nothing to add to the window
                        break;
                    case I:
                        // insertions are added to the base to the left (previous element) with the quality score of the first inserted base
                        if (locationIndex > 0) {
                            windowHeader.get(locationIndex - 1).addInsertionToTheRight();     // check if it's the first element in the read!
                            readBaseIndex += cigarElement.getLength();
                        }

                        // just ignore the insertions at the beginning of the read
                        break;
                    case D:
                        // deletions are added to the baseCounts with the read mapping quality as it's quality score
                        int nDeletionsToAdd = cigarElement.getLength();
                        while(nDeletionsToAdd-- > 0) {
                            windowHeader.get(locationIndex).addBase((byte) 'D', (byte) read.getMappingQuality(), read.getMappingQuality());
                            locationIndex++;
                        }
                        break;
                    case M:
                    case P:
                    case EQ:
                    case X:
                        int nBasesToAdd = cigarElement.getLength();
                        while(nBasesToAdd-- > 0) {
                            windowHeader.get(locationIndex).addBase(bases[readBaseIndex], quals[readBaseIndex], read.getMappingQuality());
                            readBaseIndex++;
                            locationIndex++;
                        }
                        break;
                }
            }
        }
    }

    protected HeaderElement newHeaderElement(int index) {
        return new HeaderElement(index);
    }

    /**
     * The element the composes the header of the sliding window.
     *
     * Each site has a header element containing the counts of each
     * base, it's reference based location and whether or not the site
     * has insertions (to it's right)
     */
    protected class HeaderElement {
        protected BaseCounts baseCounts;               // How many A,C,G,T (and D's) are in this site.
        protected int insertionsToTheRight;            // How many reads in this site had insertions to the immediate right
        protected int location;                        // Genome location of this site (the sliding window knows which contig we're at
        protected LinkedList<Integer> mappingQuality;   // keeps the mapping quality of each read that contributed to this element (site)


        public HeaderElement() {
            this.baseCounts = new BaseCounts();
            this.insertionsToTheRight = 0;
            this.location = 0;
            this.mappingQuality = new LinkedList<Integer>();
        }

        public int getLocation() { return location; }

        public HeaderElement(int location) {
            this();
            this.location = location;
        }

        public boolean isVariant() {
            return isVariantFromInsertions() || isVariantFromMismatches() || isVariantFromDeletions();
        }

        public void addBase(byte base, byte qual, int mappingQuality) {
            if ( qual >= MIN_BASE_QUAL_TO_COUNT )  {
                baseCounts.incr(base, qual);
            }
            this.mappingQuality.add(mappingQuality);
        }

        protected boolean isVariantFromInsertions() {
            return ((double) insertionsToTheRight / baseCounts.totalCount()) > MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT;
        }

        protected boolean isVariantFromDeletions() {
            return baseCounts.baseWithMostCounts() == BaseIndex.D.getByte() || baseCounts.baseCountProportion(BaseIndex.D) > MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT;
        }

        protected boolean isVariantFromMismatches() {
            return baseCounts.baseCountProportion(baseCounts.baseWithMostCounts()) < (1 - MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT);
        }

        public void addInsertionToTheRight() {
            insertionsToTheRight++;
        }

        public boolean isEmpty() {
            return baseCounts.totalCount() == 0;
        }

        public double getRMS() {
            return MathUtils.rms(mappingQuality);
        }

        public void addMappingQuality(int mappingQuality) {
            this.mappingQuality.add(mappingQuality);
        }
    }


}

