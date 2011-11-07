package org.broadinstitute.sting.gatk.walkers.reducereads;

import org.broadinstitute.sting.utils.sam.AlignmentStartWithNoTiesComparator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

public class ReduceReadsStash {
        protected MultiSampleCompressor compressor;
        SortedSet<GATKSAMRecord> outOfOrderReads;

        public ReduceReadsStash(MultiSampleCompressor compressor) {
            this.compressor = compressor;
            this.outOfOrderReads = new TreeSet<GATKSAMRecord>(new AlignmentStartWithNoTiesComparator());
        }

        public List<GATKSAMRecord> getAllReadsBefore(GATKSAMRecord read) {
            List<GATKSAMRecord> result = new LinkedList<GATKSAMRecord>();
            GATKSAMRecord newHead = null;

            for (GATKSAMRecord stashedRead : outOfOrderReads) {
                if (ReadUtils.compareSAMRecords(stashedRead, read) <= 0) {
                    result.add(stashedRead);
                }
                else {
                    newHead = stashedRead;
                    break;
                }
            }

            if (result.size()  > 0) {
                if (result.size() == outOfOrderReads.size())
                    outOfOrderReads.clear();
                else
                    outOfOrderReads = outOfOrderReads.tailSet(newHead);
            }

            return result;
        }

        public Iterable<GATKSAMRecord> compress(GATKSAMRecord read) {
            return compressor.addAlignment(read);
        }

        public void add(GATKSAMRecord read) {
            outOfOrderReads.add(read);
        }

        public SortedSet<GATKSAMRecord> getAllReads() {
            return outOfOrderReads;
        }

        public Iterable<GATKSAMRecord> close() {
            LinkedList<GATKSAMRecord> result = new LinkedList<GATKSAMRecord>();

            // compress all the stashed reads (in order)
            for (GATKSAMRecord read : outOfOrderReads)
                for (GATKSAMRecord compressedRead : compressor.addAlignment(read))
                    result.add(compressedRead);

            // output any remaining reads from the compressor
            for (GATKSAMRecord read : compressor.close())
                result.add(read);

            return result;
        }


}