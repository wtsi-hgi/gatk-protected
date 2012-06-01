package org.broadinstitute.sting.gatk.walkers.diagnostics;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 5/31/12
 * Time: 3:24 PM
 * To change this template use File | Settings | File Templates.
 */
@ReadFilters({UnmappedReadFilter.class,NotPrimaryAlignmentFilter.class,DuplicateReadFilter.class,FailsVendorQualityCheckFilter.class})
public class CoveredRegionWalker extends ReadWalker<GATKSAMRecord, CoveredRegionWalker.IntervalStats> {

    @Output(required = false)
    PrintStream out = System.out;

    @Argument(fullName="outputAsBED", shortName="bed", doc="Output as BED file", required=false)
    boolean outputAsBED = false;

    
    long coveredBases = 0;
    long totalBases = 0;
    long numReads = 0;
    
    @Override
    public void initialize() {
    }

    protected class ReadInterval {
        public String contig;
        public int start;
        public int stop;
        
        public ReadInterval(String c, int a, int b) {
            contig = c;
            start = a;
            stop = b;
        }
    }
    protected class IntervalStats {
        private ArrayList<ReadInterval> readIntervals;
        
        public IntervalStats() {
            readIntervals = new ArrayList<ReadInterval>();
        }
        
        public void addReadToInterval(GATKSAMRecord read) {
            final ReadInterval  newInterval = new ReadInterval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());

            if (readIntervals.isEmpty()) {
                readIntervals.add(newInterval);
                return;
            }

            final ReadInterval lastInterval = readIntervals.get(readIntervals.size()-1);
            if (!read.getReferenceName().equals(lastInterval.contig) || read.getAlignmentStart() > lastInterval.stop+1) {
                readIntervals.add(newInterval);
            }
            else {
                // read in same contig as last element and interval spanned by read can be merged with previous interval
                readIntervals.set(readIntervals.size()-1, new ReadInterval(lastInterval.contig, lastInterval.start, newInterval.stop));
            }
            
        }
        
        public String toString() {
            String result = ""; // = "Contig\tstart\tstop\n";
            
            for (ReadInterval interval: readIntervals) {
                if (outputAsBED)
                    result += String.format("%s\t%d\t%d\n", interval.contig, interval.start-1, interval.stop );
                else
                    result += String.format("%s:%d-%d\n", interval.contig, interval.start, interval.stop );
            }
            return result;
        }
    }
    @Override
    public GATKSAMRecord map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {
        return read;
    }

    @Override
    public IntervalStats reduceInit() {
        return new IntervalStats();
    }

    public IntervalStats reduce(GATKSAMRecord read, IntervalStats stats) {
        int numCoveredBases = read.getAlignmentEnd() - read.getAlignmentStart()+1;
        int unclippedSpan = read.getUnclippedEnd() - read.getUnclippedStart()+1;
   //     System.out.format("Read:%s Contig:%s, Aligned span:%d Unclipped span:%d\n", read.getReadName(), read.getReferenceName(),
     //           numCoveredBases, unclippedSpan);
        totalBases += unclippedSpan;
        coveredBases += numCoveredBases;
        numReads++;
        
        stats.addReadToInterval(read);
        return stats;
    }

    @Override
    public void onTraversalDone(IntervalStats result) {
        out.println(result);
        
        System.out.println("Total read span including only aligned Bases:"+coveredBases);
        System.out.println("Total read span including soft-clipped Bases:" + totalBases);
        System.out.println("Number of reads:"+numReads);
    }
}