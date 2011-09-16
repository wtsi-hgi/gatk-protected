package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.*;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.Iterator;
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
    private List<Byte> counts;
    private List<Byte> bases;
    private double rms;   // todo -- implement this

    // Information to produce a SAMRecord
    private SAMFileHeader header;
    private Object readGroupAttribute;
    private String contig;
    private int contigIndex;
    private String readName;
    private Integer refStart;
    final private int consensusBaseQuality;

    /**
     * Initialize your running consensus if you don't know yet what the first base and it's name
     * are going to be.
     *
     * @param header
     * @param readGroupAttribute
     * @param contig
     * @param contigIndex
     */
    public RunningConsensus (SAMFileHeader header, Object readGroupAttribute, String contig, int contigIndex, int consensusBaseQuality) {
        this(header, readGroupAttribute, contig, contigIndex, null, null, consensusBaseQuality);
    }

    /**
     * Full initialization of the running consensus if you have all the information and are ready to
     * start adding to the running consensus.
     *
     * @param header
     * @param readGroupAttribute
     * @param contig
     * @param contigIndex
     * @param readName
     * @param refStart
     */
    public RunningConsensus (SAMFileHeader header, Object readGroupAttribute, String contig, int contigIndex, String readName, Integer refStart, int consensusBaseQuality) {
        counts = new LinkedList<Byte>();
        bases = new LinkedList<Byte>();
        rms = 0.0;

        this.header = header;
        this.readGroupAttribute = readGroupAttribute;
        this.contig = contig;
        this.contigIndex = contigIndex;
        this.readName = readName;
        this.refStart = refStart;
        this.consensusBaseQuality = consensusBaseQuality;
    }

    /**
     * Easy access to keep adding to a running consensus that has already been
     * initialized with the correct read name and refStart
     *
     * @param base
     * @param count
     */
    public void add(byte base, byte count, double rms) {
        counts.add(count);
        bases.add(base);
        this.rms += rms;
    }

    public SAMRecord close () {
        SAMRecord samRecord = new SAMRecord(header);
        samRecord.setAttribute("RG", readGroupAttribute);
        samRecord.setAttribute(ReadUtils.REDUCED_READ_QUALITY_TAG, consensusBaseQuality);
        samRecord.setReferenceName(contig);
        samRecord.setReferenceIndex(contigIndex);
        samRecord.setReadPairedFlag(false);
        samRecord.setReadUnmappedFlag(false);
        samRecord.setAlignmentStart(refStart);
        samRecord.setCigar(buildCigar());
        samRecord.setReadName(readName);
        samRecord.setBaseQualities(convertBaseQualities());
        samRecord.setReadBases(convertReadBases());
        samRecord.setMappingQuality((int) Math.ceil(rms/bases.size()));
        return samRecord;
    }

    public int size () {
        return bases.size();
    }

    private byte [] convertBaseQualities() {
        return listToByteArray(counts);
    }

    private byte [] convertReadBases() {
        return listToByteArray(bases);
    }

    private byte [] listToByteArray(List<Byte> list) {
        byte [] array = new byte[list.size()];
        int i = 0;
        Iterator<Byte> basesIterator = bases.listIterator();
        for (Byte element : list) {
            Byte b = basesIterator.next();
            switch (BaseIndex.byteToBase(b)) {
                case D:  // do not add deletions to the consensus SAM record
                    break;
                default:
                    array[i++] = element;
                    break;
            }
        }
        if (i < list.size()) {
            byte[] arrayWithoutDeletions = new byte[i];
            System.arraycopy(array, 0, arrayWithoutDeletions, 0, i);
            return arrayWithoutDeletions;
        }
        return array;
    }

    private Cigar buildCigar() {
        LinkedList<CigarElement> cigarElements = new LinkedList<CigarElement>();
        CigarOperator cigarOperator = null;
        int length = 0;
        for (Byte b : bases) {
            CigarOperator op;
            switch (BaseIndex.byteToBase(b)) {
                case D:
                    op = CigarOperator.DELETION;
                    break;
                case I:
                    op = CigarOperator.INSERTION;
                    break;
                default:
                    op = CigarOperator.MATCH_OR_MISMATCH;
                    break;
            }
            if (cigarOperator == null)
                cigarOperator = op;

            else if (cigarOperator != op) {      // need to treat 1st case
                cigarElements.add(new CigarElement(length, cigarOperator));
                cigarOperator = op;
                length = 0;
            }
            length++;
        }
        if (length > 0)
            cigarElements.add(new CigarElement(length, cigarOperator));
        return new Cigar(cigarElements);
    }
}
