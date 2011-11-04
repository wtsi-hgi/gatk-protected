package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.SAMReadGroupRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.sam.AlignmentStartWithNoTiesComparator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.TreeSet;

/**
 *
 * @author depristo
 * @version 0.1
 */
public class SingleSampleConsensusReadCompressor implements ConsensusReadCompressor {
    protected static final Logger logger = Logger.getLogger(SingleSampleConsensusReadCompressor.class);

    protected final int contextSize;
    protected final int contextSizeIndels;
    protected final int downsampleCoverage;
    protected int minMappingQuality;
    protected int slidingWindowCounter;

    protected final String sampleName;
    protected final SAMReadGroupRecord reducedReadGroup;

    protected SlidingWindow slidingWindow;
    protected double minAltProportionToTriggerVariant;
    protected double minIndelProportionToTriggerVariant;
    protected int minBaseQual;
    protected int maxQualCount;


    public SingleSampleConsensusReadCompressor(final String sampleName,
                                               final SAMReadGroupRecord readGroupRecord,
                                               final int contextSize,
                                               final int contextSizeIndels,
                                               final int downsampleCoverage,
                                               final int minMappingQuality,
                                               final double minAltProportionToTriggerVariant,
                                               final double minIndelProportionToTriggerVariant,
                                               final int minBaseQual,
                                               final int maxQualCount) {
        this.sampleName = sampleName;
        this.reducedReadGroup = readGroupRecord;
        this.contextSize = contextSize;
        this.contextSizeIndels = contextSizeIndels;
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
    public Iterable<GATKSAMRecord> addAlignment( GATKSAMRecord read ) {
        TreeSet<GATKSAMRecord> result = new TreeSet<GATKSAMRecord>(new AlignmentStartWithNoTiesComparator());
        int position = read.getUnclippedStart();

        // create a new window if:
        if ((slidingWindow != null) &&
            ( ( read.getReferenceIndex() != slidingWindow.getContigIndex() ) ||     // this is a brand new contig
              (position - contextSize > slidingWindow.getStopLocation()) ) ) {  // this read is too far away from the end of the current sliding window

            // close the current sliding window
            result.addAll(slidingWindow.close());
            slidingWindow = null;                                                   // so we create a new one on the next if
        }

        if ( slidingWindow == null) {       // this is the first read
            instantiateSlidingWindow(read);
            slidingWindowCounter++;
        }

        result.addAll(slidingWindow.addRead(read));
        return result;
    }

    protected void instantiateSlidingWindow(GATKSAMRecord read) {
        slidingWindow = new SlidingWindow(read.getReferenceName(), read.getReferenceIndex(), contextSize, contextSizeIndels, read.getHeader(), read.getAttribute("RG"), slidingWindowCounter, minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, maxQualCount, minMappingQuality);
    }

    @Override
    public Iterable<GATKSAMRecord> close() {
        return (slidingWindow != null) ? slidingWindow.close() : new TreeSet<GATKSAMRecord>();
    }

}

