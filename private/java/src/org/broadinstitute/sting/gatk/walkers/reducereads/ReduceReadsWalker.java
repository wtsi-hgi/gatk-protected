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
import net.sf.samtools.SAMUtils;
import net.sf.samtools.util.SequenceUtil;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.sam.SimplifyingSAMFileWriter;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: April 7, 2011
 */

@PartitionBy(PartitionType.INTERVAL)
@ReadFilters({UnmappedReadFilter.class,NotPrimaryAlignmentFilter.class,DuplicateReadFilter.class,FailsVendorQualityCheckFilter.class})
public class ReduceReadsWalker extends ReadWalker<List<SAMRecord>, ConsensusReadCompressor> {

    @Output
    protected StingSAMFileWriter out;

    @Argument(fullName = "context_size", shortName = "cs", doc = "", required = false)
    protected int contextSize = 20;

    @Argument(fullName = "minimum_mapping_quality", shortName = "minmap", doc = "", required = false)
    protected int minMappingQuality = 20;

    @Argument(fullName = "minimum_tail_qualities", shortName = "mintail", doc = "", required = false)
    protected byte minTailQuality = 2;

    @Argument(fullName = "minimum_alt_proportion_to_trigger_variant", shortName = "minvar", doc = "", required = false)
    protected double minAltProportionToTriggerVariant = 0.05;

    @Argument(fullName = "minimum_del_proportion_to_trigger_variant", shortName = "mindel", doc = "", required = false)
    protected double minIndelProportionToTriggerVariant = 0.01;

    @Argument(fullName = "minimum_base_quality_to_consider", shortName = "minqual", doc = "", required = false)
    protected int minBaseQual = 20;

    @Argument(fullName = "maximum_consensus_base_qual", shortName = "maxqual", doc = "", required = false)
    protected byte maxQualCount = SAMUtils.MAX_PHRED_SCORE;

    @Hidden
    @Argument(fullName = "", shortName = "dl", doc = "", required = false)
    protected int debugLevel = 0;

    @Hidden //todo -- not yet implemented
    @Argument(fullName = "downsample_coverage", shortName = "ds", doc = "", required = false)
    protected int downsampleCoverage = 500;


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
    private List<SAMRecord> hardClipReadToInterval(SAMRecord read) {
        List<SAMRecord> clippedReads = new LinkedList<SAMRecord>();
        ReadClipper clipper = new ReadClipper(read);
        boolean foundInterval = false;

        while (!foundInterval) {
            foundInterval = true;                 // we will only need to look again if there is no overlap.
            ReadUtils.ReadAndIntervalOverlap overlapType = ReadUtils.getReadAndIntervalOverlapType(read, currentInterval);
            switch (overlapType) {
                case NO_OVERLAP_CONTIG:           // check the next interval
                case NO_OVERLAP_RIGHT:            // check the next interval
                    foundInterval = false;
                    break;

                case NO_OVERLAP_LEFT:             // read shouldn't be here in the first place
                    throw new ReviewedStingException("Read does not overlap interval, read out of order? " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd() + " + " + read.getCigarString());

                case NO_OVERLAP_HARDCLIPPED_LEFT: // read used to overlap but got hard clipped and doesn't overlap anymore
                case NO_OVERLAP_HARDCLIPPED_RIGHT:// read used to overlap but got hard clipped and doesn't overlap anymore
                    clippedReads.add(new SAMRecord(read.getHeader()));
                    break;

                case OVERLAP_LEFT:                // clip the left end of the read
                    clippedReads.add(clipper.hardClipByReferenceCoordinatesLeftTail(currentInterval.getStart() - 1));
                    break;

                case OVERLAP_RIGHT:               // clip the right end of the read
                    clippedReads.add(clipper.hardClipByReferenceCoordinatesRightTail(currentInterval.getStop() + 1));
                    break;

                case OVERLAP_LEFT_AND_RIGHT:      // clip both left and right ends of the read
                    clippedReads.add(clipper.hardClipBothEndsByReferenceCoordinates(currentInterval.getStart()-1, currentInterval.getStop()+1));
                    break;

                case OVERLAP_CONTAINED:           // don't do anything to the read
                    clippedReads.add(read);
                    break;
            }

            // If there is no overlap due to contig or because the read was to the right of the current interval, we need to get the next interval
            // Because the reads are sorted we should only traverse the interval list once for the entire genome.
            if (!foundInterval) {
                if (intervalIterator.hasNext())
                    currentInterval = intervalIterator.next();
                else {
                    clippedReads.add(new SAMRecord(read.getHeader()));
                    break;
                }
            }
        }
        return clippedReads;
    }

    @Override
    public void initialize() {
        super.initialize();

        compressor = new MultiSampleConsensusReadCompressor(getToolkit().getSAMFileHeader(), contextSize, downsampleCoverage, minMappingQuality, minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, maxQualCount);

        //todo -- should be TRUE
        out.setPresorted(false);

        for ( SAMReadGroupRecord rg : compressor.getReducedReadGroups())
            out.getFileHeader().addReadGroup(rg);

        // Keep track of the interval list so we can filter out reads that are not within the
        // requested intervals
        if (getToolkit().getIntervals() != null) {
            intervalIterator = getToolkit().getIntervals().iterator();
            currentInterval = intervalIterator.next();
        }

        if ( maxQualCount > SAMUtils.MAX_PHRED_SCORE )
            throw new UserException.BadArgumentValue("maximum_consensus_base_qual", "Maximum allowed quality score in a SAM file is " + SAMUtils.MAX_PHRED_SCORE);

    }

    @Override
    public List<SAMRecord> map( ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker ) {
        List<SAMRecord> mappedReads;
        totalReads++;
        read = SimplifyingSAMFileWriter.simplifyRead(read);

        if (debugLevel == 1) System.out.printf("Original: %s %s %d %d\n", read, read.getCigar(), read.getAlignmentStart(), read.getAlignmentEnd());

        ReadClipper clipper = new ReadClipper(read);
        SAMRecord clippedRead = clipper.hardClipLowQualEnds(minTailQuality);

        clipper = new ReadClipper(clippedRead);
        clippedRead = clipper.hardClipSoftClippedBases();

        clipper = new ReadClipper(clippedRead);
        clippedRead = clipper.hardClipLeadingInsertions();

        if (intervalIterator != null && clippedRead.getReadLength() > 0)
            mappedReads = hardClipReadToInterval(clippedRead);
        else
            mappedReads = new LinkedList<SAMRecord>();

        if (debugLevel == 1) {
            System.out.println("Output Reads:");
            for (SAMRecord mappedRead : mappedReads)
                System.out.printf("%s %d %d\n\n", mappedRead.getCigar(), mappedRead.getAlignmentStart(), mappedRead.getAlignmentEnd());
        }

        return mappedReads;

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
     * @param mappedReads all reads mapped from the original read
     * @return the SAMFileWriter, so that the next reduce can emit to the same source
     */
    public ConsensusReadCompressor reduce( List<SAMRecord> mappedReads, ConsensusReadCompressor compressor ) {
        // write out compressed reads as they become available
        for (SAMRecord read : mappedReads)
            if (read.getReadLength() != 0)
                for ( SAMRecord compressedRead : compressor.addAlignment(read))
                    outputRead(compressedRead);
        return compressor;
    }

    @Override
    public void onTraversalDone( ConsensusReadCompressor compressor ) {
        // write out any remaining reads
        for ( SAMRecord compressedRead : compressor.close() )
            outputRead(compressedRead);
    }

    private void checkForHighMismatch(SAMRecord read) {
        final int start = read.getAlignmentStart();
        final int stop = read.getAlignmentEnd();
        final byte[] ref = getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(read.getReferenceName(), start, stop).getBases();
        final int nm = SequenceUtil.countMismatches(read, ref, start - 1);
        final int readLen = read.getReadLength();
        final double nmFraction = nm / (1.0*readLen);
        if ( nmFraction > 0.4 && readLen > 20 && read.getAttribute(ReadUtils.REDUCED_READ_QUALITY_TAG) != null)
            throw new ReviewedStingException("BUG: High mismatch fraction found in read " + read.getReadName() + " position: " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd());
    }

    private boolean isConsensus(SAMRecord read) {
        return read.getAttribute(ReadUtils.REDUCED_READ_QUALITY_TAG) != null;
    }

    private void outputRead(SAMRecord read) {
        if (debugLevel == 2)
            checkForHighMismatch(read);

        if (isConsensus(read))
            nCompressedReads++;
        else
            totalReads++;

        out.addAlignment(read);

    }


}
