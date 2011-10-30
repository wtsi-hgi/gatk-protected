package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.sam.AlignmentStartWithNoTiesComparator;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

public class ReduceReadsStash {
        protected MultiSampleConsensusReadCompressor compressor;
        SortedSet<SlidingRead> outOfOrderReads;

        public ReduceReadsStash(MultiSampleConsensusReadCompressor compressor) {
            this.compressor = compressor;
            this.outOfOrderReads = new TreeSet<SlidingRead>(new SlidingReadComparator());
        }

        public List<SlidingRead> getAllReadsBefore(SlidingRead slidingRead) {
            List<SlidingRead> result = new LinkedList<SlidingRead>();
            SlidingRead newHead = null;

            for (SlidingRead stashedRead : outOfOrderReads) {
                if (ReadUtils.compareSAMRecords(stashedRead.getRead(), slidingRead.getRead()) <= 0) {
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

        public Iterable<SAMRecord> compress(SlidingRead slidingRead) {
            return compressor.addAlignment(slidingRead);
        }

        public void add(SlidingRead read) {
            outOfOrderReads.add(read);
        }

        public SortedSet<SlidingRead> getAllReads() {
            return outOfOrderReads;
        }

        public Iterable<SAMRecord> close() {
            LinkedList<SAMRecord> result = new LinkedList<SAMRecord>();

            // compress all the stashed reads (in order)
            for (SlidingRead read : outOfOrderReads)
                for (SAMRecord compressedRead : compressor.addAlignment(read))
                    result.add(compressedRead);

            // output any remaining reads from the compressor
            for (SAMRecord read : compressor.close())
                result.add(read);

            return result;
        }


}