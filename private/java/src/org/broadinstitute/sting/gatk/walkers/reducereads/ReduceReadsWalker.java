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
import org.broadinstitute.sting.commandline.Hidden;
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
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.clipreads.ClippingOp;
import org.broadinstitute.sting.utils.clipreads.ClippingRepresentation;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;

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

    @Output(fullName="bedOut", shortName = "bedOut", doc="BED output", required = false)
    protected PrintStream bedOut = null;

    @Argument(fullName = "contextSize", shortName = "CS", doc = "", required = false)
    protected int contextSize = 10;

    @Argument(fullName = "AverageDepthAtVariableSites", shortName = "ADAV", doc = "", required = false)
    protected int AverageDepthAtVariableSites = 500;

    @Argument(fullName = "ReadQualityEquivalent", shortName = "QE", doc = "", required = false)
    protected int QUALITY_EQUIVALENT = 20;

    @Argument(fullName = "MinimumMappingQuality", shortName = "MM", doc = "", required = false)
    protected int MIN_MAPPING_QUALITY = 20;

    @Hidden
    @Argument(fullName = "INCLUDE_RAW_READS", shortName = "IRR", doc = "", required = false)
    protected boolean INCLUDE_RAW_READS = false;

    @Hidden
    @Argument(fullName = "useRead", shortName = "UR", doc = "", required = false)
    protected Set<String> readNamesToUse;



    protected int totalReads = 0;
    int nCompressedReads = 0;

    GenomeLocSortedSet intervals = null;
    Iterator<GenomeLoc> i = null;
    GenomeLoc interval = null;

    MultiSampleConsensusReadCompressor compressor;

    @Override
    public void initialize() {
        super.initialize();

        compressor = new MultiSampleConsensusReadCompressor(getToolkit().getSAMFileHeader(),
                contextSize, getToolkit().getGenomeLocParser(),
                AverageDepthAtVariableSites, QUALITY_EQUIVALENT, MIN_MAPPING_QUALITY);

        out.setPresorted(false);

        for ( SAMReadGroupRecord rg : compressor.getReducedReadGroups())
            out.getFileHeader().addReadGroup(rg);
        intervals = getToolkit().getIntervals();
        i = intervals.iterator();
        interval = i.next();

    }

    @Override
    public SAMRecord map( ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker ) {
        totalReads++;
        return read; // all the work is done in the reduce step for this walker
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

    private void findIntersectingInterval (SAMRecord read) {
        while ( i.hasNext() ){
            if ( (interval.getContig() == read.getReferenceName()) && (interval.getStart() < read.getUnclippedEnd()) && (interval.getStop() > read.getUnclippedStart()) )
                     break;
            else {
                    interval = i.next();
            }
            // TODO what happens if interval is not found?
        }
    }

    /**
     * given a read and a output location, reduce by emitting the read
     * @param read the read itself
     * @return the SAMFileWriter, so that the next reduce can emit to the same source
     */
    public ConsensusReadCompressor reduce( SAMRecord read, ConsensusReadCompressor comp ) {
        if ( readNamesToUse == null || readNamesToUse.contains(read.getReadName()) ) {
            if ( INCLUDE_RAW_READS )
                out.addAlignment(read);
            findIntersectingInterval(read);

            ReadClipper clipper = new ReadClipper(read);
            if(read.getUnclippedEnd() > interval.getStop())
                clipper.addOp( new ClippingOp( interval.getStop() - read.getUnclippedStart() + 1 , read.getReadLength() - 1 ));
            if(read.getUnclippedStart() < interval.getStart())
                clipper.addOp( new ClippingOp( 0 , interval.getStart() - read.getUnclippedStart() - 1 ));

            // must be +2 here because findSpan needs to include the last loci of the interval

            read = clipper.clipRead(ClippingRepresentation.HARDCLIP_BASES);

            // write out compressed reads as they become available
            for ( SAMRecord consensusRead : comp.addAlignment(read)) {
                out.addAlignment(consensusRead);
                System.out.println(String.format("Output Read: %d-%d, Cigar: %s, NAME: %s", consensusRead.getAlignmentStart(), consensusRead.getAlignmentEnd(), consensusRead.getCigarString(), consensusRead.getReadName()));                nCompressedReads++;
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
            System.out.println(String.format("Output Read: %d-%d, CIGAR: %s, NAME: %s", consensusRead.getAlignmentStart(), consensusRead.getAlignmentEnd(), consensusRead.getCigarString(), consensusRead.getReadName()));
            nCompressedReads++;
        }

        double percent = (100.0 * nCompressedReads) / totalReads;
        logger.info("Compressed reads : " + nCompressedReads + String.format(" (%.2f%%)", percent));
        logger.info("Total reads      : " + totalReads);
    }
    
}
