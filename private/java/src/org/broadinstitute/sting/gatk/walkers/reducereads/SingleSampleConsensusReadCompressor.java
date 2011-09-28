package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;

import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author depristo
 * @version 0.1
 */
public class SingleSampleConsensusReadCompressor implements ConsensusReadCompressor {
    protected static final Logger logger = Logger.getLogger(SingleSampleConsensusReadCompressor.class);

    private final int readContextSize;
    private final int downsampleCoverage;
    private int minMappingQuality;
    private int slidingWindowCounter;

    private final String sampleName;
    private final SAMReadGroupRecord reducedReadGroup;

    private SlidingWindow slidingWindow;
    private double minAltProportionToTriggerVariant;
    private double minIndelProportionToTriggerVariant;
    private int minBaseQual;
    private int maxQualCount;


    public SingleSampleConsensusReadCompressor(final String sampleName,
                                               final SAMReadGroupRecord readGroupRecord,
                                               final int readContextSize,
                                               final int downsampleCoverage,
                                               final int minMappingQuality,
                                               final double minAltProportionToTriggerVariant,
                                               final double minIndelProportionToTriggerVariant,
                                               final int minBaseQual,
                                               final int maxQualCount) {
        this.sampleName = sampleName;
        this.reducedReadGroup = readGroupRecord;
        this.readContextSize = readContextSize;
        this.downsampleCoverage = downsampleCoverage;
        this.minMappingQuality = minMappingQuality;
        this.slidingWindowCounter = 0;
        this.minAltProportionToTriggerVariant = minAltProportionToTriggerVariant;
        this.minIndelProportionToTriggerVariant = minIndelProportionToTriggerVariant;
        this.minBaseQual = minBaseQual;
        this.maxQualCount = maxQualCount;
    }

    public SAMReadGroupRecord getReducedReadGroup() {
        return reducedReadGroup;
    }

    /**
     * @{inheritDoc}
     */
    @Override
    public Iterable<SAMRecord> addAlignment( SAMRecord read ) {
        List<SAMRecord> result = new LinkedList<SAMRecord>();
        int position = read.getUnclippedStart();

        // create a new window if:
        if ((slidingWindow != null) &&
            ( ( read.getReferenceIndex() != slidingWindow.getContigIndex() ) ||     // this is a brand new contig
              (position - readContextSize > slidingWindow.getStopLocation()) ) ) {  // this read is too far away from the end of the current sliding window

            // close the current sliding window
            result.addAll(slidingWindow.close());
            slidingWindow = null;                                                   // so we create a new one on the next if
        }

        if ( slidingWindow == null) {       // this is the first read
            slidingWindow = new SlidingWindow(read.getReferenceName(), read.getReferenceIndex(), readContextSize,
                                              read.getHeader(), read.getAttribute("RG"), slidingWindowCounter,
                                              minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, maxQualCount, minMappingQuality);
            slidingWindowCounter++;
        }

        result.addAll(slidingWindow.addRead(read));
        return result;
    }

    @Override
    public List<SAMRecord> close() {
        return (slidingWindow != null) ? slidingWindow.close() : new LinkedList<SAMRecord>();
    }

}

