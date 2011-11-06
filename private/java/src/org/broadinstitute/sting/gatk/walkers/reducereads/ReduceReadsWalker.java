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
import org.broadinstitute.sting.utils.GenomeLocComparator;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

@PartitionBy(PartitionType.INTERVAL)
@ReadFilters({UnmappedReadFilter.class,NotPrimaryAlignmentFilter.class,DuplicateReadFilter.class,FailsVendorQualityCheckFilter.class})
public class ReduceReadsWalker extends ReadWalker<List<GATKSAMRecord>, ReduceReadsStash> {

    @Output
    protected StingSAMFileWriter out;

    @Argument(fullName = "context_size", shortName = "cs", doc = "", required = false)
    protected int contextSize = 10;

    @Argument(fullName = "context_size_indels", shortName = "csindel", doc = "", required = false)
    protected int contextSizeIndels = 50;

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

    @Hidden
    @Argument(fullName = "", shortName = "dr", doc = "", required = false)
    protected String debugRead = "";

    @Hidden //todo -- not yet implemented
    @Argument(fullName = "downsample_coverage", shortName = "ds", doc = "", required = false)
    protected int downsampleCoverage = 500;


    protected int totalReads = 0;
    int nCompressedReads = 0;

    HashMap<String, Long> readNameHash;  // This hash will keep the name of the original read the new compressed name (a number).
    Long nextReadNumber = 1L;             // The next number to use for the compressed read name.

    SortedSet<GenomeLoc> intervalList;


    /**
     * Hard clips the read around the edges of the interval it overlaps with.
     * Note: If read overlaps more than one interval, it will be hard clipped at the end of the first interval it overlaps.
     *   (maybe split in two reads and treat it specially in the future?)
     *
     * @param read the read to be hard clipped to the interval.
     * @return a shallow copy of the read hard clipped to the interval
     */
    private List<GATKSAMRecord> hardClipReadToInterval(GATKSAMRecord read) {
        List<GATKSAMRecord> clippedReads = new LinkedList<GATKSAMRecord>();
        ReadClipper clipper = new ReadClipper(read);

        GenomeLoc intervalOverlapped = null;       // marks the interval to which the original read overlapped (so we can cut all previous intervals from the list)

        boolean originalRead = true;               // false if this is the right tail of the original read
        boolean overlap;                           // keeps track of the interval that overlapped the original read
        boolean doneClipping;                      // triggers an early exit if we are done clipping this read

        if (intervalList.isEmpty())
            clippedReads.add(read);                // if we don't have intervals (wgs) the read goes in unchanged

        for (GenomeLoc interval : intervalList) {

            if ( read.getReadLength() == 0 )       // nothing to do with an empty read (could have been fully clipped before)
                break;

            GATKSAMRecord clippedRead = null;                 // this will hold the read clipped to the interval to be added in the end of the switch

            switch (ReadUtils.getReadAndIntervalOverlapType(read, interval)) {
                case NO_OVERLAP_RIGHT:             // no reads on this interval, check the next interval if this is the original read
                    if (!originalRead)             // something went wrong if this is the tail of the read
                        throw new ReviewedStingException("tail of the read should never NO_OVERLAP_RIGHT the following interval. " + read.getReadName() + " -- " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd() + " x " + interval.getLocation().toString());
                    overlap = false;
                    doneClipping = false;
                    break;


                case NO_OVERLAP_HARDCLIPPED_RIGHT: // read used to overlap but got hard clipped and doesn't overlap anymore
                    if (originalRead) {
                        overlap = true;            // effectively, we have found the read's location and now we are going to try and match it's tail (which happens to be the entire read).
                        clippedRead = new GATKSAMRecord(read.getHeader());
                    }
                    else
                        overlap = false;

                    doneClipping = false;
                    break;

                case NO_OVERLAP_CONTIG:            // read is in a different contig
                    if (originalRead) {            // the original read can be in a bigger contig, but not on a smaller one.
                        if (read.getReferenceIndex() < interval.getContigIndex())
                            throw new ReviewedStingException("read is behind interval list. (contig) " + read.getReadName() + " -- " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd() + " x " + interval.getLocation().toString());
                        else {
                            overlap = false;
                            doneClipping = false;
                        }
                    }                              // tail read CANNOT be in a different contig.
                    else {
                        if (read.getReferenceIndex() < interval.getContigIndex()) {
                            overlap = false;
                            doneClipping = true;
                        }
                        else
                            throw new ReviewedStingException("Tail read is in bigger contig than interval traversal. " + read.getReadName() + " -- " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd() + " x " + interval.getLocation().toString());

                    }
                    break;

                case NO_OVERLAP_LEFT:
                    if (originalRead)              // if this is the first read this should never happen.
                        throw new ReviewedStingException("original read cannot be behind the first interval. (position) " + read.getReadName() + " -- " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd() + " x " + interval.getLocation().toString());

                    overlap = false;
                    doneClipping = true;
                    break;

                case NO_OVERLAP_HARDCLIPPED_LEFT:  // read used to overlap but got hard clipped and doesn't overlap anymore
                    overlap = originalRead;        // if this is the original read, we should not advance the interval list, the original overlap was here.
                    doneClipping = true;
                    break;

                case OVERLAP_LEFT:                 // clip the left tail of the read
                    clippedRead = clipper.hardClipByReferenceCoordinatesLeftTail(interval.getStart() - 1);

                    overlap = true;
                    doneClipping = true;
                    break;

                case OVERLAP_RIGHT:                // clip the right tail of the read and try to match it to the next interval
                    clippedRead = clipper.hardClipByReferenceCoordinatesRightTail(interval.getStop() + 1);
                    read = clipper.hardClipByReferenceCoordinatesLeftTail(interval.getStop());
                    clipper = new ReadClipper(read);

                    overlap = true;
                    doneClipping = false;
                    break;

                case OVERLAP_LEFT_AND_RIGHT:       // clip both left and right ends of the read
                    clippedRead = clipper.hardClipBothEndsByReferenceCoordinates(interval.getStart()-1, interval.getStop()+1);
                    read = clipper.hardClipByReferenceCoordinatesLeftTail(interval.getStop());
                    clipper = new ReadClipper(read);

                    overlap = true;
                    doneClipping = false;
                    break;

                case OVERLAP_CONTAINED:            // don't do anything to the read
                    clippedRead = read;

                    overlap = true;
                    doneClipping = true;
                    break;

                default:
                    throw new ReviewedStingException("interval overlap returned an unknown / unhandled state. If new state was added to intervalOverlap, it should be handled by hardClipReadToInterval.");
            }

            if (overlap && originalRead)
                intervalOverlapped = interval;

            if (clippedRead != null) {
                originalRead = false;

                if (clippedRead.getReadLength() > 0)
                    clippedReads.add(clippedRead);     // if the read overlaps the interval entirely within a deletion, it will be entirely clipped off
            }

            if (doneClipping)
                break;
        }

        if (intervalOverlapped != null)
            intervalList = intervalList.tailSet(intervalOverlapped);

        return clippedReads;
    }


    @Override
    public void initialize() {
        super.initialize();

        readNameHash = new HashMap<String, Long>();

        intervalList = new TreeSet<GenomeLoc>(new GenomeLocComparator ());
        intervalList.addAll(getToolkit().getIntervals());

        out.setPresorted(false);

//        for ( SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups())
//            out.getFileHeader().addReadGroup(rg);

        if ( maxQualCount > SAMUtils.MAX_PHRED_SCORE )
            throw new UserException.BadArgumentValue("maximum_consensus_base_qual", "Maximum allowed quality score in a SAM file is " + SAMUtils.MAX_PHRED_SCORE);
    }

    @Override
    public List<GATKSAMRecord> map( ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker ) {
        totalReads++;
        if (!debugRead.isEmpty() && read.getReadName().contains(debugRead))
            System.out.println("Found debug read!");

        if (debugLevel == 1) System.out.printf("\nOriginal: %s %s %d %d\n", read, read.getCigar(), read.getAlignmentStart(), read.getAlignmentEnd());

        ReadClipper clipper = new ReadClipper(read);
        GATKSAMRecord clippedRead = clipper.hardClipLowQualEnds(minTailQuality);

        clipper = new ReadClipper(clippedRead);
        clippedRead = clipper.hardClipSoftClippedBases();

        clipper = new ReadClipper(clippedRead);
        clippedRead = clipper.hardClipLeadingInsertions();

        List<GATKSAMRecord> mappedReads = hardClipReadToInterval(clippedRead);

        if (debugLevel == 1) {
            for (GATKSAMRecord mappedRead : mappedReads)
                System.out.printf("MAPPED: %s %d %d\n", mappedRead.getCigar(), mappedRead.getAlignmentStart(), mappedRead.getAlignmentEnd());
        }

        return mappedReads;

    }


    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     * @return SAMFileWriter, set to the BAM output file if the command line option was set, null otherwise
     */
    @Override
    public ReduceReadsStash reduceInit() {
        return new ReduceReadsStash(new MultiSampleCompressor(getToolkit().getSAMFileHeader(), contextSize, contextSizeIndels, downsampleCoverage, minMappingQuality, minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, maxQualCount));
    }

    /**
     * given a read and an output location, reduce by emitting the read
     */
    public ReduceReadsStash reduce( List<GATKSAMRecord> mappedReads, ReduceReadsStash stash ) {
        boolean firstRead = true;
        for (GATKSAMRecord read : mappedReads) {
            boolean originalRead = firstRead && ReadUtils.getReadAndIntervalOverlapType(read, intervalList.first()) == ReadUtils.ReadAndIntervalOverlap.OVERLAP_CONTAINED;

            if (read.getReadLength() == 0)
                throw new ReviewedStingException("Empty read sent to reduce, this should never happen! " + read.getReadName() + " -- " + read.getCigar() + " -- " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd() );

            if (originalRead) {
                List<GATKSAMRecord> readsReady = new LinkedList<GATKSAMRecord>();
                readsReady.addAll(stash.getAllReadsBefore(read));
                readsReady.add(read);

                for (GATKSAMRecord readReady : readsReady) {
                    if (debugLevel == 1) System.out.println("REDUCE: " + readReady.getCigar() + " " + readReady.getAlignmentStart() + " " + readReady.getAlignmentEnd());

                    for ( GATKSAMRecord compressedRead : stash.compress(readReady))
                        outputRead(compressedRead);

                }
            }
            else
                stash.add(read);

            firstRead = false;
        }

        return stash;
    }

    @Override
    public void onTraversalDone(ReduceReadsStash stash) {

        // output any remaining reads in the compressor
        for ( GATKSAMRecord read : stash.close() )
            outputRead(read);
    }

    private void checkForHighMismatch(GATKSAMRecord read) {
        final int start = read.getAlignmentStart();
        final int stop = read.getAlignmentEnd();
        final byte[] ref = getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(read.getReferenceName(), start, stop).getBases();
        final int nm = SequenceUtil.countMismatches(read, ref, start - 1);
        final int readLen = read.getReadLength();
        final double nmFraction = nm / (1.0*readLen);
        if ( nmFraction > 0.4 && readLen > 20 && read.getAttribute(GATKSAMRecord.REDUCED_READ_QUALITY_TAG) != null)
            throw new ReviewedStingException("BUG: High mismatch fraction found in read " + read.getReadName() + " position: " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd());
    }

    private boolean isConsensus(GATKSAMRecord read) {
        return read.getAttribute(GATKSAMRecord.REDUCED_READ_QUALITY_TAG) != null;
    }

    private void outputRead(GATKSAMRecord read) {
        if (debugLevel == 2)
            checkForHighMismatch(read);

        if (isConsensus(read))
            nCompressedReads++;
        else
            totalReads++;

        if (debugLevel == 1) System.out.println("BAM: " + read.getCigar() + " " + read.getAlignmentStart() + " " + read.getAlignmentEnd());

        out.addAlignment(compressReadName(read));
    }

    private GATKSAMRecord compressReadName(GATKSAMRecord read) {
        GATKSAMRecord compressedRead;

        try {
            compressedRead = (GATKSAMRecord) read.clone();
        } catch (CloneNotSupportedException e) {
            throw new ReviewedStingException("Where did the clone go?");
        }

        String name = read.getReadName();
        String compressedName = isConsensus(read) ? "C" : "";
        if (readNameHash.containsKey(name))
            compressedName += readNameHash.get(name).toString();
        else {
            readNameHash.put(name, nextReadNumber);
            compressedName += nextReadNumber.toString();
            nextReadNumber++;
        }

        compressedRead.setReadName(compressedName);
        return compressedRead;
    }
}
