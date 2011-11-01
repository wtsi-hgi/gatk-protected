package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.sam.AlignmentStartWithNoTiesComparator;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

public class ReduceReadsStash {
        protected MultiSampleConsensusReadCompressor compressor;
        SortedSet<SAMRecord> outOfOrderReads;

        public ReduceReadsStash(MultiSampleConsensusReadCompressor compressor) {
            this.compressor = compressor;
            this.outOfOrderReads = new TreeSet<SAMRecord>(new AlignmentStartWithNoTiesComparator());
        }

        public List<SAMRecord> getAllReadsBefore(SAMRecord read) {
            List<SAMRecord> result = new LinkedList<SAMRecord>();
            SAMRecord newHead = null;

            for (SAMRecord stashedRead : outOfOrderReads) {
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

        public Iterable<SAMRecord> compress(SAMRecord read) {
            return compressor.addAlignment(read);
        }

        public void add(SAMRecord read) {
            outOfOrderReads.add(read);
        }

        public SortedSet<SAMRecord> getAllReads() {
            return outOfOrderReads;
        }

        public Iterable<SAMRecord> close() {
            LinkedList<SAMRecord> result = new LinkedList<SAMRecord>();

            // compress all the stashed reads (in order)
            for (SAMRecord read : outOfOrderReads)
                for (SAMRecord compressedRead : compressor.addAlignment(read))
                    result.add(compressedRead);

            // output any remaining reads from the compressor
            for (SAMRecord read : compressor.close())
                result.add(read);

            return result;
        }


}