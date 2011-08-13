package org.broadinstitute.sting.gatk.walkers.reducereads;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.clipreads.ClippingOp;
import org.broadinstitute.sting.utils.clipreads.ClippingRepresentation;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import sun.management.counter.Variability;

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
    // TODO WE WANT TO PUT ALL function is SLlidingWindow
    // TODO comment out unused code
    private SAMFileHeader header;
    private final int readContextSize;
    private final int AverageDepthAtVariableSites;
    private final int QualityEquivalent;
    private final int minMapQuality;
    private int consensusCounter = 0;

    private final SAMReadGroupRecord reducedReadGroup;
    private String contig = null;
    private final GenomeLocParser glParser;

    private SlidingWindow slidingWindow;


    public SingleSampleConsensusReadCompressor(final String sampleName,
                                               final int readContextSize,
                                               final GenomeLocParser glParser,
                                               final int AverageDepthAtVariableSites,
                                               final int qualityEquivalent,
                                               final int minMapQuality) {
        this.readContextSize = readContextSize;
        this.glParser = glParser;
        this.slidingWindow = new SlidingWindow("SampleName",contig, header);
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
        if ( contig == null )
            contig = read.getReferenceName();
        if ( ! read.getReferenceName().equals(contig) )
            throw new ReviewedStingException("ConsensusRead system doesn't support multiple contig processing right now");


        if ( header == null )
            header = read.getHeader();

        /*
        if ( ! waitingReads.isEmpty() && read.getAlignmentStart() < waitingReads.peek().getAlignmentStart() )
            throw new ReviewedStingException(
                    String.format("Adding read %s starting at %d before current queue head start position %d",
                            read.getReadName(), read.getAlignmentStart(), waitingReads.peek().getAlignmentStart()));
        */
        List<SAMRecord> result = new LinkedList<SAMRecord>();
        /*
        if ( retryTimer == 0 ) {

            if ( chunkReadyForConsensus(read) ) {
                result = consensusReads(false);
            }
        } else {
            //logger.info("Retry: " + retryTimer);
            retryTimer--;
        }

        waitingReads.add(read);
        */

        int position = read.getAlignmentStart();
        logger.info(String.format("Setting position to %d", position));slidingWindow.addRead(read);

        // did adding the read create variance?
        List<VariableRegion> variableRegions = slidingWindow.getVariableRegions(readContextSize);
        for ( VariableRegion variableRegion : variableRegions ) {
            logger.info(String.format("Found a variable region : %d - %d", variableRegion.start, variableRegion.end) );
            if ( (position - readContextSize) >= variableRegion.start ) {
                result.addAll(slidingWindow.finalizeConsensusRead(variableRegion));
            }
            if ( position > variableRegion.end ) {
                result.addAll(slidingWindow.finalizeVariableRegion(variableRegion));
            }
        }
        if ( variableRegions.isEmpty() )
            slidingWindow.compressWindow(position-readContextSize);
        else
            slidingWindow.compressWindow(Math.min( variableRegions.get(0).start, position - readContextSize ));
        return result;
    }

    @Override
    public Iterable<SAMRecord> close() {
        // nothing needs ot happen
        LinkedList<SAMRecord> result = new LinkedList<SAMRecord>();
        for ( VariableRegion variableRegion : slidingWindow.getVariableRegions(readContextSize) ) {
            logger.info(String.format("Variable region at close() : %d - %d", variableRegion.start, variableRegion.end) );
            result.addAll(slidingWindow.finalizeVariableRegion(variableRegion));

        }
        logger.info(String.format("Finalizing LAST Consensus Read at %d", slidingWindow.getEnd()) );
        result.addAll(slidingWindow.finalizeConsensusRead(new VariableRegion(-1,slidingWindow.getEnd() + 1)));
        return result;
    }

    // ------------------------------------------------------------------------------------------
    //
    // NO private implementation functions
    //
    // ------------------------------------------------------------------------------------------


}

