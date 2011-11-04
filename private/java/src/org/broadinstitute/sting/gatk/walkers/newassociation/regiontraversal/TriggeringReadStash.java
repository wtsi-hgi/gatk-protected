package org.broadinstitute.sting.gatk.walkers.newassociation.regiontraversal;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.walkers.reducereads.*;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.AlignmentStartWithNoTiesComparator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 10/11/11
 * Time: 7:48 PM
 * To change this template use File | Settings | File Templates.
 */
public class TriggeringReadStash extends ReduceReadsStash {
    int nSamples;

    protected TriggeringReadStash(TriggeringMultiSampleCompressor compressor) {
        super(compressor);
    }

    public TriggeringReadStash(SAMFileHeader header) {
        // todo -- take inputs
        this(new TriggeringMultiSampleCompressor(header,80,80,100,20,0.15,0.05,15,100));
        nSamples = SampleUtils.getSAMFileSamples(header).size();
    }

    @Override
    public Iterable<GATKSAMRecord> compress(GATKSAMRecord read) {
        return super.compress(read);
    }
}

class TriggeringMultiSampleCompressor extends MultiSampleCompressor {

    public TriggeringMultiSampleCompressor(SAMFileHeader header,
                                           final int contextSize,
                                           final int contextSizeIndels,
                                           final int downsampleCoverage,
                                           final int minMappingQuality,
                                           final double minAltProportionToTriggerVariant,
                                           final double minIndelProportionToTriggerVariant,
                                           final int minBaseQual,
                                           final int maxQualCount) {
        super(header,contextSize,contextSizeIndels,downsampleCoverage,minMappingQuality,minAltProportionToTriggerVariant,minIndelProportionToTriggerVariant,minBaseQual,maxQualCount);
        compressorsPerSample.clear(); // out with the bad, in with the good
        for ( String name : SampleUtils.getSAMFileSamples(header) ) {
            compressorsPerSample.put(name,
                    new TriggeringSingleSampleCompressor(name, header.getReadGroup(name), contextSize, contextSizeIndels, downsampleCoverage,
                                    minMappingQuality, minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, maxQualCount,this));
        }
    }

    protected boolean[] markRegion(int startLoc, int stopLoc, int contextSize) {
        boolean[] multiSampleRegion = new boolean[stopLoc-startLoc+contextSize+1];
        for ( Map.Entry<String,SingleSampleCompressor> compressor : compressorsPerSample.entrySet() ) {
            ((TriggeringSingleSampleCompressor) compressor.getValue()).markWherePossible(multiSampleRegion,startLoc,stopLoc,contextSize);
        }

        return multiSampleRegion;
    }

    @Override
    public Iterable<GATKSAMRecord> addAlignment(GATKSAMRecord read) {
        String sample = read.getReadGroup().getSample();
        SingleSampleCompressor compressor = compressorsPerSample.get(sample);
        if ( compressor == null )
            throw new ReviewedStingException("No compressor for sample " + sample);
        TreeSet<GATKSAMRecord> toReturn = new TreeSet<GATKSAMRecord>(new AlignmentStartWithNoTiesComparator());
        for ( Map.Entry<String,SingleSampleCompressor> compressorEntry : compressorsPerSample.entrySet() ) {
            if ( compressorEntry.getKey().equals(sample) ) {
                toReturn.addAll( (Collection<GATKSAMRecord>) compressorEntry.getValue().addAlignment(read));
            } else {
                toReturn.addAll( (Collection<GATKSAMRecord>) ( (TriggeringSingleSampleCompressor) compressorEntry.getValue()).registerAlignment(read));
            }
        }

        return toReturn;
    }
}

class TriggeringSingleSampleCompressor extends SingleSampleCompressor {
    protected TriggeringMultiSampleCompressor parentMultiSampleCompressor;
    public TriggeringSingleSampleCompressor(final String sampleName,
                                            final SAMReadGroupRecord readGroupRecord,
                                            final int contextSize,
                                            final int contextSizeIndels,
                                            final int downsampleCoverage,
                                            final int minMappingQuality,
                                            final double minAltProportionToTriggerVariant,
                                            final double minIndelProportionToTriggerVariant,
                                            final int minBaseQual,
                                            final int maxQualCount,
                                            final TriggeringMultiSampleCompressor parent) {
        super(sampleName,readGroupRecord,contextSize,contextSizeIndels,downsampleCoverage,minMappingQuality,minAltProportionToTriggerVariant,minIndelProportionToTriggerVariant,minBaseQual,maxQualCount);
        parentMultiSampleCompressor = parent;
    }

    @Override
    protected void instantiateSlidingWindow(GATKSAMRecord read) {
        slidingWindow = new TriggeringSlidingWindow(read.getReferenceName(), read.getReferenceIndex(), contextSize, contextSizeIndels, read.getHeader(), read.getAttribute("RG"), slidingWindowCounter, minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, maxQualCount, minMappingQuality,this);
    }

    protected boolean[] markRegion(int startLoc, int stopLoc, int contextSize) {
        return parentMultiSampleCompressor.markRegion(startLoc, stopLoc, contextSize);
    }

    protected void markWherePossible(boolean[] variantArray, int startLoc, int stopLoc, int contextSize) {
        if ( slidingWindow != null )
            ( (TriggeringSlidingWindow) slidingWindow).markWherePossible(variantArray, startLoc, stopLoc, contextSize);
    }

    public Iterable<GATKSAMRecord> registerAlignment(GATKSAMRecord read) {
       return slidingWindow == null ? new ArrayList<GATKSAMRecord>(0) : ( (TriggeringSlidingWindow) slidingWindow).registerAlignment(read);
    }
}

class TriggeringSlidingWindow extends SlidingWindow {
    protected final static double MIN_LARGE_INSERT_READS_TO_TRIGGER_VARIANT = 0.10;
    protected final static double MIN_CLIPPED_READS_TO_TRIGGER_VARIANT = 0.10;
    protected TriggeringSingleSampleCompressor parentCompressor;

    public TriggeringSlidingWindow(String contig, int contigIndex, int contextSize, int contextSizeIndels,
                                   SAMFileHeader header, Object readGroupAttribute, int windowNumber,
                                   final double minAltProportionToTriggerVariant,
                                   final double minIndelProportionToTriggerVariant, int minBaseQual,
                                   int maxQualCount, int minMappingQuality,
                                   TriggeringSingleSampleCompressor parent) {
        super(contig, contigIndex, contextSize, contextSizeIndels, header, readGroupAttribute, windowNumber, minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, maxQualCount, minMappingQuality);
        parentCompressor = parent;
    }

    @Override
    protected boolean [] markSites(int stop) {
        return parentCompressor.markRegion(startLocation,stop,contextSize);
        /* old code
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
        */
    }

    protected void markWherePossible(boolean[] arrayToMark, int locArrayStart, int locArrayStop, int conSize) {
        int myStartOffset = Math.max(locArrayStart-startLocation,0);
        int arrayStartOffset = Math.max(startLocation-locArrayStart,0);
        Iterator<HeaderElement> headerElements = windowHeader.iterator();
        // may need to discard header elements
        while ( myStartOffset > 0 && headerElements.hasNext()) {
            headerElements.next();
            myStartOffset--;
        }
        // for remaining header elements (to array stop at most), mark if there are variant sites, plus contexts
        int index = 0;
        while ( headerElements.hasNext() && index < arrayToMark.length && index + arrayStartOffset + startLocation < locArrayStop ) {
            TriggeringHeaderElement elem = ( TriggeringHeaderElement ) headerElements.next();
            if ( elem.isVariant() ) {
                int from = Math.max(arrayStartOffset+index-conSize,0);
                int to = Math.min(arrayStartOffset+index+conSize,arrayToMark.length);
                for ( int i = from; i < to; i ++ ) {
                    arrayToMark[i] = true;
                }
            }
        }
    }

    @Override
    protected List<GATKSAMRecord> addToConsensus(int start, int end) {
        // reads that would die from the slide to start get added in
        List<GATKSAMRecord> consensus = new LinkedList<GATKSAMRecord>();
        return consensus;
    }

    @Override
    protected GATKSAMRecord finalizeConsensus() {
        runningConsensus = null;
        return null;
    }

    @Override
    protected HeaderElement newHeaderElement(int location) {
        return new TriggeringHeaderElement(location);
    }

    @Override
        protected void updateHeaderCounts(GATKSAMRecord read) {
        // Reads that don't pass the minimum mapping quality filter are not added to the
        // consensus, or count towards a variant region so no point in keeping track of
        // their base counts.
        if (read.getMappingQuality() < MIN_MAPPING_QUALITY)
            return;

        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();
        Cigar cigar = read.getCigar();
        boolean largeInsert = read.getAttribute("LI") != null;

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

        // todo -- perhaps rewrite this iteration using list iterator to save time searching for each index.
        // todo -- they should be consecutive (as far as I can tell)
        Boolean firstElement = null;
        for (CigarElement cigarElement : cigar.getCigarElements()) {
            firstElement = firstElement == null;
            switch (cigarElement.getOperator()) {
                case H:
                case S:
                    // nothing to add to the window
                    int basesToAdd = cigarElement.getLength();
                    int innerLocIndex = firstElement ? Math.max(0,locationIndex-cigarElement.getLength()) : locationIndex;
                    while ( basesToAdd -- > 0 ) {
                        if ( innerLocIndex + read.getAlignmentStart() < stopLocation && quals[readBaseIndex] >= 20 ) {
                            ((TriggeringHeaderElement) windowHeader.get(innerLocIndex)).addClippedBase();
                        }
                        if ( innerLocIndex + read.getAlignmentStart() < stopLocation && largeInsert ) {
                            ((TriggeringHeaderElement) windowHeader.get(innerLocIndex)).addLargeInsert();
                        }
                        readBaseIndex++;
                        innerLocIndex++;
                    }
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
                        if ( largeInsert )
                            ((TriggeringHeaderElement) windowHeader.get(locationIndex)).addLargeInsert();
                        readBaseIndex++;
                        locationIndex++;
                    }
                    break;
            }
        }
    }

    public List<GATKSAMRecord> registerAlignment( GATKSAMRecord read ) {
        // If this is the first read in the window, update startLocation
        if (startLocation < 0)
            startLocation = read.getAlignmentStart();

        List<GATKSAMRecord> finalizedReads = slideIfPossible(read.getUnclippedStart());

        return finalizedReads;
    }

    protected class TriggeringHeaderElement extends HeaderElement {

        protected int readsLargeInsert = 0;
        protected int highQClippedBases = 0;

        public TriggeringHeaderElement() {
            super();
        }

        public TriggeringHeaderElement(int location) {
            super(location);
        }

        public void addLargeInsert() {
            readsLargeInsert++;
        }

        public void addClippedBase() {
            highQClippedBases++;
        }

        @Override
        public boolean isVariant() {
            return super.isVariant() || isVariantFromClippedReads() || isVariantFromInsertSize();
        }

        public boolean isVariantFromInsertSize() {
            return ( readsLargeInsert/totalReads() ) > MIN_LARGE_INSERT_READS_TO_TRIGGER_VARIANT;
        }

        public boolean isVariantFromClippedReads() {
            return (highQClippedBases/totalReads()) > MIN_CLIPPED_READS_TO_TRIGGER_VARIANT;
        }

        protected double totalReads() {
            return baseCounts.totalCount() + highQClippedBases;
        }
    }
}
