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

package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckReadFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentReadFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.clipreads.ClippingOp;
import org.broadinstitute.sting.utils.clipreads.ClippingRepresentation;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: April 7, 2011
 */

@ReadFilters({UnmappedReadFilter.class,NotPrimaryAlignmentReadFilter.class,DuplicateReadFilter.class,FailsVendorQualityCheckReadFilter.class})
public class ReduceReadsWalker extends ReadWalker<SAMRecord, ConsensusReadCompressor> {

    @Output
    protected StingSAMFileWriter out;

    @Argument(fullName = "contextSize", shortName = "CS", doc = "", required = false)
    protected int contextSize = 10;

    @Argument(fullName = "AverageDepthAtVariableSites", shortName = "ADAV", doc = "", required = false)
    protected int AverageDepthAtVariableSites = 500;

    @Argument(fullName = "ReadQualityEquivalent", shortName = "QE", doc = "", required = false)
    protected int QUALITY_EQUIVALENT = 20;

    @Argument(fullName = "MinimumMappingQuality", shortName = "MM", doc = "", required = false)
    protected int MIN_MAPPING_QUALITY = 20;

    protected int totalReads = 0;
    int nCompressedReads = 0;

    Iterator<GenomeLoc> intervalIterator = null;
    GenomeLoc currentInterval = null;

    MultiSampleConsensusReadCompressor compressor;


    /**
     * Hard clips the read around the edges of the interval it overlaps with.
     * Note: If read overlaps more than one interval, it will be hard clipped at the end of the first interval it overlaps.
     *   (maybe split in two reads and treat it specially in the future?)
     *
     * @param read the read to be hard clipped to the interval.
     * @return a shallow copy of the read hard clipped to the interval
     */
    private SAMRecord hardClipReadToInterval(SAMRecord read) {
        ReadClipper clipper = new ReadClipper(read);

        boolean clipRead = false;
        boolean foundInterval = false;

        while (!foundInterval) {
            foundInterval = true;        // we will only need to look again if there is no overlap.
            ReadUtils.ReadAndIntervalOverlap overlapType = ReadUtils.getReadAndIntervalOverlapType(read, currentInterval);
            switch (overlapType) {
                case NO_OVERLAP:         // check the next interval
                    foundInterval = false;
                    break;

                case LEFT_OVERLAP:       // clip the left end of the read
                    clipper.addOp( new ClippingOp( 0 , currentInterval.getStart() - read.getUnclippedStart() - 1 ));
                    clipRead = true;
                    break;

                case RIGHT_OVERLAP:      // clip the right end of the read
                    clipper.addOp( new ClippingOp( currentInterval.getStop() - read.getUnclippedStart() + 1 , read.getReadLength() - 1 ));
                    clipRead = true;
                    break;

                case FULL_OVERLAP:       // clip both left and right ends of the read
                    clipper.addOp( new ClippingOp( 0 , currentInterval.getStart() - read.getUnclippedStart() - 1 ));
                    clipper.addOp( new ClippingOp( currentInterval.getStop() - read.getUnclippedStart() + 1 , read.getReadLength() - 1 ));
                    clipRead = true;
                    break;

                case CONTAINED:          // don't do anything to the read
                    break;
            }

            // If there is no overlap, we need to get the next interval
            // Because the reads are sorted we should only traverse the interval list once for the entire genome.
            if (!foundInterval) {
                if (intervalIterator.hasNext())
                    currentInterval = intervalIterator.next();
                else
                    throw new ReviewedStingException("Read is over the last requested interval. Either the reads are not sorted or the GATK Engine is not filtering reads outside the requested interval");
            }

        }
        return (clipRead) ? clipper.clipRead(ClippingRepresentation.HARDCLIP_BASES) : read;
    }





    @Override
    public void initialize() {
        super.initialize();

        compressor = new MultiSampleConsensusReadCompressor(getToolkit().getSAMFileHeader(),
                contextSize, getToolkit().getGenomeLocParser(),
                AverageDepthAtVariableSites, QUALITY_EQUIVALENT, MIN_MAPPING_QUALITY);

        //todo -- should be TRUE
        out.setPresorted(false);

        for ( SAMReadGroupRecord rg : compressor.getReducedReadGroups())
            out.getFileHeader().addReadGroup(rg);

        // Keep track of the interval list so we can filter out reads that are not within the
        // requested intervals
        intervalIterator = getToolkit().getIntervals().iterator();
        currentInterval = intervalIterator.next();

    }

    @Override
    public SAMRecord map( ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker ) {
        totalReads++;

        // If the user provided a list of intervals, hard clip the reads to the intervals
        return (!getToolkit().getIntervals().isEmpty()) ? hardClipReadToInterval(read) : read;
    }


    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     * @return SAMFileWriter, set to the BAM output file if the command line option was set, null otherwise
     */
    @Override
    public ConsensusReadCompressor reduceInit() {
        return compressor;
    }

    /**
     * given a read and a output location, reduce by emitting the read
     * @param read the read itself
     * @return the SAMFileWriter, so that the next reduce can emit to the same source
     */
    public ConsensusReadCompressor reduce( SAMRecord read, ConsensusReadCompressor comp ) {
        // write out compressed reads as they become available
        for ( SAMRecord consensusRead : comp.addAlignment(read)) {
            out.addAlignment(consensusRead);
            System.out.println(String.format("Output Read: %d-%d, Cigar: %s, NAME: %s", consensusRead.getAlignmentStart(), consensusRead.getAlignmentEnd(), consensusRead.getCigarString(), consensusRead.getReadName()));                nCompressedReads++;
        }
        return comp;
    }

    @Override
    public void onTraversalDone( ConsensusReadCompressor compressor ) {
        //compressor.writeConsensusBed(bedOut);
        // write out any remaining reads

        for ( SAMRecord consensusRead : compressor.close() ) {
            out.addAlignment(consensusRead);
            System.out.println(String.format("Output Read: %d-%d, CIGAR: %s, NAME: %s", consensusRead.getAlignmentStart(), consensusRead.getAlignmentEnd(), consensusRead.getCigarString(), consensusRead.getReadName()));
            nCompressedReads++;
        }

        double percent = (100.0 * nCompressedReads) / totalReads;
        logger.info("Compressed reads : " + nCompressedReads + String.format(" (%.2f%%)", percent));
        logger.info("Total reads      : " + totalReads);
    }
    
}
