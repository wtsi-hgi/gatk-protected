package org.broadinstitute.sting.gatk.walkers.diagnostics;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckReadFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentReadFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.ContigComparator;


import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 7/23/11
 * Time: 1:18 PM
 * To change this template use File | Settings | File Templates.
 */

@ReadFilters({UnmappedReadFilter.class,NotPrimaryAlignmentReadFilter.class,DuplicateReadFilter.class,FailsVendorQualityCheckReadFilter.class})
public class ContigStatsWalker extends ReadWalker<SAMRecord, ContigStatsWalker.ContigStats> {

    @Output(required = false)
    PrintStream out = System.out;

    protected class ContigStats {
        Map<String, Long> contigCount;
        long totalReads;

        public ContigStats() {
            this.contigCount = new HashMap<String, Long>();
            this.totalReads = 0;
        }

        public Map<String, Long> getContigStat() {
            return contigCount;
        }

        public Long getTotalReads() {
            return totalReads;
        }

        public Long getReadsAtContig(String contig) {
            return contigCount.get(contig);
        }

        public double getPercentReadsAtContig(String contig) {
            return (double) contigCount.get(contig) / totalReads;
        }

        public void addReadToContig(String contig) {
            Long count = 1L;
            if (contigCount.containsKey(contig))
                count += contigCount.get(contig);
            contigCount.put(contig, count);
            totalReads++;
        }

        public Long getExpectedReadsAtContig (String contig) {
            int contigLength = getToolkit().getSAMFileHeader().getSequenceDictionary().getSequence(contig).getSequenceLength();
            long totalLength = 0;
            for (SAMSequenceRecord c : getToolkit().getSAMFileHeader().getSequenceDictionary().getSequences())
                totalLength += c.getSequenceLength();
            double contigProportion = (double) contigLength/totalLength;
            return Math.round(contigProportion * totalReads);
        }

        public double getExpectedPercetReadsAtContig(String contig) {
            return (double) getExpectedReadsAtContig(contig)/totalReads;
        }

        public String getContigStats (String contig) {
            return contig + "\t" +
                    getToolkit().getSAMFileHeader().getSequenceDictionary().getSequence(contig).getSequenceLength() + "\t" +
                    getReadsAtContig(contig) + "\t" +
                    getExpectedReadsAtContig(contig) + "\t" +
                    String.format("%.3f", getPercentReadsAtContig(contig)) + "\t" +
                    String.format("%.3f", getExpectedPercetReadsAtContig(contig)) + "\n";
        }

        public String toString() {
            String result = "contig" + "\t" + "size" + "\t" + "reads" + "\t" + "exp_reads" + "\t" + "enrichment" +"\t" + "exp_enrichment\n";

            // create a sorted set for the contigs to always be output in the same order (using Contig Comparator for ordering)
            SortedSet<String> contigs = new TreeSet<String>(new ContigComparator());
            contigs.addAll(contigCount.keySet());

            for (String contig : contigs)
                result += getContigStats(contig);
            return result;
        }
    }

    @Override
    public void initialize() {
    }

    @Override
    public SAMRecord map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        return read;
    }

    @Override
    public ContigStats reduceInit() {
        return new ContigStats();
    }

    public ContigStats reduce(SAMRecord read, ContigStats stats) {
        stats.addReadToContig(read.getReferenceName());
        return stats;
    }

    @Override
    public void onTraversalDone(ContigStats result) {
        out.println(result);
    }
}
