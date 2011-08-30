package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;

import java.util.LinkedList;
import java.util.List;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 *
 * @author depristo
 * @version 0.1
 */
public class SingleSampleConsensusReadCompressor implements ConsensusReadCompressor {
    protected static final Logger logger = Logger.getLogger(SingleSampleConsensusReadCompressor.class);
    /*
    private static final boolean DEBUG = false;
    private static final boolean INVERT = false;
    private static final boolean PRINT_CONSENSUS_READS = false;
    private static final int CYCLES_BEFORE_RETRY = 1000;
    private static final double MAX_FRACTION_DISAGREEING_BASES = 0.1;
    private static final ClippingRepresentation VARIABLE_READ_REPRESENTATION = ClippingRepresentation.SOFTCLIP_BASES;
    private static final double MIN_FRACT_BASES_FOR_VARIABLE_READ = 0.33;  // todo -- should be variable
    private static final int MIN_BASES_IN_VARIABLE_SPAN_TO_INCLUDE_READ = 10;
    */

    private final static String RG_POSTFIX = ".ReducedReads";

    // todo  -- should merge close together spans
    // TODO WE WANT TO PUT ALL functions is SlidingWindow
    // TODO comment out unused code
    private final int readContextSize;
    private final int AverageDepthAtVariableSites;
    private int minMappingQuality;
    private int slidingWindowCounter;

    private final SAMReadGroupRecord reducedReadGroup;

    private SlidingWindow slidingWindow;
    private double minAltProportionToTriggerVariant;
    private int minBaseQual;
    private int maxQualCount;


    public SingleSampleConsensusReadCompressor(final String sampleName,
                                               final int readContextSize,
                                               final int AverageDepthAtVariableSites,
                                               final int minMappingQuality,
                                               final double minAltProportionToTriggerVariant,
                                               final int minBaseQual,
                                               final int maxQualCount) {
        this.readContextSize = readContextSize;
        this.AverageDepthAtVariableSites = AverageDepthAtVariableSites;
        this.reducedReadGroup = createReducedReadGroup(sampleName);
        this.minMappingQuality = minMappingQuality;
        this.slidingWindowCounter = 0;
        this.minAltProportionToTriggerVariant = minAltProportionToTriggerVariant;
        this.minBaseQual = minBaseQual;
        this.maxQualCount = maxQualCount;
    }

    /**
     * Helper function to create a read group for these reduced reads
     * @param sampleName
     * @return
     */
    private static final SAMReadGroupRecord createReducedReadGroup(final String sampleName) {
        SAMReadGroupRecord rg = new SAMReadGroupRecord(sampleName + RG_POSTFIX);
        rg.setSample(sampleName);
        return rg;
    }

    public SAMReadGroupRecord getReducedReadGroup() {
        return reducedReadGroup;
    }

    // ------------------------------------------------------------------------------------------
    //
    // public interface functions
    //
    // ------------------------------------------------------------------------------------------

    /**
     * @{inheritDoc}
     */
    @Override
    public Iterable<SAMRecord> addAlignment( SAMRecord read ) {

        /**
         * Steps to adding alignment:
         *
         * 1 - Slide the current window if we can.
         * 2 - If necessary, close and create a new window.
         * 3 - Add read to the window
         *
         *
         * When adding read to the window, handle 3 types of reads :
         *  1 - Read is fully contained by current sliding window
         *    * add read to the sliding window -- simple case
         *  2 - Read is partially contained by current sliding window (left side is inside, right side is outside)
         *    * increase the sliding window to accommodate the new read
         *  3 - Read is completely outside of the sliding window
         *    * close current sliding window and create a new one
         *
         * Notes:
         *  - The starting position of the incoming read must be the unclipped end because reads may have been trimmed
         *  and the next read start can be smaller than the current incoming read's alignmentStart, but never smaller
         *  than the unclippedStart.
         *
         *  - The size of the sliding window has to include the context size (to the right).
         *
         */

        List<SAMRecord> result = new LinkedList<SAMRecord>();
        int position = read.getUnclippedStart();

        // create a new window if:
        if ((slidingWindow != null) &&
            ( (!read.getReferenceName().equals(slidingWindow.getContig())) ||     // this is a brand new contig
              (position - readContextSize > slidingWindow.getStopLocation()))) {  // this read is too far away from the end of the current sliding window

            // close the current sliding window
            result.addAll(slidingWindow.close());
            slidingWindow = null;  // so we create a new one on the next if
        }

        if ( slidingWindow == null) {       // this is the first read
            slidingWindow = new SlidingWindow(read.getReferenceName(), read.getReferenceIndex(), readContextSize,
                                              read.getHeader(), read.getAttribute("RG"), slidingWindowCounter,
                                              minAltProportionToTriggerVariant, minBaseQual, maxQualCount, minMappingQuality);
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

