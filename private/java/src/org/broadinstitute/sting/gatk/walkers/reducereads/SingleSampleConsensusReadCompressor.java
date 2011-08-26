package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

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
    private SAMFileHeader header;
    private final int readContextSize;
    private final int AverageDepthAtVariableSites;
    private final int QualityEquivalent;
    private final int minMapQuality;
    private int consensusCounter = 0;

    private final SAMReadGroupRecord reducedReadGroup;
    private String contig = null;

    private SlidingWindow slidingWindow;


    public SingleSampleConsensusReadCompressor(final String sampleName,
                                               final int readContextSize,
                                               final int AverageDepthAtVariableSites,
                                               final int qualityEquivalent,
                                               final int minMapQuality) {
        this.readContextSize = readContextSize;
        this.AverageDepthAtVariableSites = AverageDepthAtVariableSites;
        this.reducedReadGroup = createReducedReadGroup(sampleName);
        this.QualityEquivalent = qualityEquivalent;
        this.minMapQuality = minMapQuality;
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
        if ( contig == null ||                                                                                                  // this is the first read
             !read.getReferenceName().equals(contig) ||                                                                         // this is a brand new contig
             ( position - readContextSize > slidingWindow.getStopLocation() && slidingWindow.getStopLocation() != -1 ))   {     // this read is too far away from the end of the current sliding window

            // if the current sliding window exists, close it
            if (slidingWindow != null)
                result.addAll(slidingWindow.close());
            slidingWindow = new SlidingWindow(read.getReadGroup().getSample(), read.getReferenceName(), read.getHeader(), readContextSize);
        }
        result.addAll(slidingWindow.addRead(read));

        return result;
    }

    @Override
    public List<SAMRecord> close() {
        return (slidingWindow != null) ? slidingWindow.close() : new LinkedList<SAMRecord>();
    }

}

