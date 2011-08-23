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
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.sam.SimplifyingSAMFileWriter;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: April 7, 2011
 */

@ReadFilters({UnmappedReadFilter.class,NotPrimaryAlignmentFilter.class,DuplicateReadFilter.class,FailsVendorQualityCheckFilter.class})
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

    @Argument(fullName = "minimum_tail_qualities", shortName = "mtq", doc = "", required = false)
    protected byte minTailQuality = 2;

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
        SAMRecord clippedRead = read;
        ReadClipper clipper = new ReadClipper(read);
        boolean foundInterval = false;

        while (!foundInterval) {
            foundInterval = true;        // we will only need to look again if there is no overlap.
            ReadUtils.ReadAndIntervalOverlap overlapType = ReadUtils.getReadAndIntervalOverlapType(read, currentInterval);
            switch (overlapType) {
                case NO_OVERLAP_CONTIG:         // check the next interval
                case NO_OVERLAP_RIGHT:          // check the next interval
                    foundInterval = false;
                    break;

                case NO_OVERLAP_LEFT:           // read used to overlap but got hard clipped and doesn't overlap anymore
                    return new SAMRecord(read.getHeader());

                case OVERLAP_LEFT:              // clip the left end of the read
                    clippedRead = clipper.hardClipByReferenceCoordinates(-1 , currentInterval.getStart() - 1);
                    break;

                case OVERLAP_RIGHT:             // clip the right end of the read
                    clippedRead = clipper.hardClipByReferenceCoordinates(currentInterval.getStop() + 1, -1);
                    break;

                case OVERLAP_LEFT_AND_RIGHT:    // clip both left and right ends of the read
                    clippedRead = clipper.hardClipBothEndsByReferenceCoordinates(currentInterval.getStart()-1, currentInterval.getStop()+1);
                    break;

                case OVERLAP_CONTAINED:         // don't do anything to the read
                    break;
            }

            // If there is no overlap due to contig or because the read was to the right of the current interval, we need to get the next interval
            // Because the reads are sorted we should only traverse the interval list once for the entire genome.
            if (!foundInterval) {
                if (intervalIterator.hasNext())
                    currentInterval = intervalIterator.next();
                else
                    return new SAMRecord(read.getHeader());
            }

        }
        return clippedRead;
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
        read = SimplifyingSAMFileWriter.simplifyRead(read);

        ReadClipper clipper = new ReadClipper(read);
        //byte [] bases = read.getReadBases();
        //System.out.printf("\nOriginal: %s %s %d %d\n", read, read.getCigar(), read.getAlignmentStart(), read.getAlignmentEnd());
        /*

        for (byte base : bases) System.out.print((char)base);
        System.out.println();
        */
        SAMRecord filteredRead = clipper.hardClipLowQualEnds(minTailQuality);

        if (clipper.wasClipped()) {
            //System.out.printf("******* WAS CLIPPED ***********\n");
            //bases = filteredRead.getReadBases();
            //for (byte base : bases) System.out.print((char)base);
            //System.out.println();
        }

        SAMRecord clippedRead = read;
        if (filteredRead.getReadLength() > 0 && !getToolkit().getIntervals().isEmpty())
            clippedRead = hardClipReadToInterval(filteredRead);

        //System.out.printf("Result: %s %d %d  => %s %d %d => %s %d %d\n", read.getCigar(), read.getAlignmentStart(), read.getAlignmentEnd(), filteredRead.getCigar(), filteredRead.getAlignmentStart(), filteredRead.getAlignmentEnd(), clippedRead.getCigar(), clippedRead.getAlignmentStart(), clippedRead.getAlignmentEnd());
        /*
        bases = clippedRead.getReadBases();
        for (byte base : bases) System.out.print((char)base);
        System.out.println();
        */
        return clippedRead;

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
        if (read.getReadLength() != 0) {
            // write out compressed reads as they become available
            for ( SAMRecord consensusRead : comp.addAlignment(read)) {
                out.addAlignment(consensusRead);
                //System.out.println(String.format("Output Read: %d-%d, Cigar: %s, NAME: %s", consensusRead.getAlignmentStart(), consensusRead.getAlignmentEnd(), consensusRead.getCigarString(), consensusRead.getReadName()));                nCompressedReads++;
                nCompressedReads++;
            }
        }
        return comp;
    }

    @Override
    public void onTraversalDone( ConsensusReadCompressor compressor ) {
        //compressor.writeConsensusBed(bedOut);
        // write out any remaining reads

        for ( SAMRecord consensusRead : compressor.close() ) {
            out.addAlignment(consensusRead);
            //System.out.println(String.format("Output Read: %d-%d, CIGAR: %s, NAME: %s", consensusRead.getAlignmentStart(), consensusRead.getAlignmentEnd(), consensusRead.getCigarString(), consensusRead.getReadName()));
            nCompressedReads++;
        }

        double percent = (100.0 * nCompressedReads) / totalReads;
        logger.info("Compressed reads : " + nCompressedReads + String.format(" (%.2f%%)", percent));
        logger.info("Total reads      : " + totalReads);
    }
    
}
