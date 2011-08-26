package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

import java.util.LinkedList;
import java.util.List;

/**
 * Running Consensus is a read that is composed as a sliding window travels over the reads
 * and keeps track of all the bases that are outside of variant regions.
 *
 * Consensus reads have qual fields that correspond to the number of reads that had the base
 * and passed the minimum quality threshold.
 *
 * The mapping quality of a consensus read is the RMS of the mapping qualities of all reads
 * that compose the consensus
 *
 * @author Mauricio Carneiro
 * @since 8/26/11
 */
public class RunningConsensus {
    private double rms;
    private List<Byte> counts;
    private List<Byte> bases;

    // Information to produce a SAMRecord
    private SAMFileHeader header;
    private Object readGroupRecord;



    public RunningConsensus (SAMFileHeader header, Object readGroupRecord) {
        counts = new LinkedList<Byte>();
        bases = new LinkedList<Byte>();
        rms = 0.0;
        this.header = header;
        this.readGroupRecord = readGroupRecord;
    }

    public void add(byte base, byte count) {
        counts.add(count);
        bases.add(base);
    }

    public SAMRecord close () {
        SAMRecord samRecord = new SAMRecord(header);
        samRecord.setAttribute("RG", readGroupRecord);
        return samRecord;
    }
}
