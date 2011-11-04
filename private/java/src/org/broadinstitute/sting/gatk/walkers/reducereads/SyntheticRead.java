package org.broadinstitute.sting.gatk.walkers.reducereads;

import com.google.java.contract.Requires;
import net.sf.samtools.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Running Consensus is a read that is compressed as a sliding window travels over the reads
 * and keeps track of all the bases that are outside of variant regions.
 *
 * Consensus reads have qual fields that correspond to the number of reads that had the base
 * and passed the minimum quality threshold.
 *
 * The mapping quality of a consensus read is the average RMS of the mapping qualities of all reads
 * that compose the consensus
 *
 * @author Mauricio Carneiro
 * @since 8/26/11
 */
public class SyntheticRead {
    private List<Byte> bases;
    private List<Byte> counts;
    private List<Byte> quals;
    private double mappingQuality;          // the average of the rms of the mapping qualities of all the reads that contributed to this consensus

    // Information to produce a GATKSAMRecord
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
    public SyntheticRead(SAMFileHeader header, Object readGroupAttribute, String contig, int contigIndex, int consensusBaseQuality) {
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
    public SyntheticRead(SAMFileHeader header, Object readGroupAttribute, String contig, int contigIndex, String readName, Integer refStart, int consensusBaseQuality) {
        bases = new LinkedList<Byte>();
        counts = new LinkedList<Byte>();
        quals = new LinkedList<Byte>();
        mappingQuality = 0.0;

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
    @Requires("count < Byte.MAX_VALUE")
    public void add(byte base, byte count, byte qual, double mappingQuality) {
        counts.add(count);
        bases.add(base);
        quals.add(qual);
        this.mappingQuality += mappingQuality;
    }

    public GATKSAMRecord close () {
        GATKSAMRecord samRecord = new GATKSAMRecord(header);
        samRecord.setReferenceName(contig);
        samRecord.setReferenceIndex(contigIndex);
        samRecord.setReadPairedFlag(false);
        samRecord.setReadUnmappedFlag(false);
        samRecord.setAlignmentStart(refStart);
        samRecord.setCigar(buildCigar());
        samRecord.setReadName(readName);
        samRecord.setBaseQualities(convertBaseQualities());
        samRecord.setReadBases(convertReadBases());
        samRecord.setMappingQuality((int) Math.ceil(mappingQuality /bases.size()));
        samRecord.setAttribute("RG", readGroupAttribute);
        samRecord.setAttribute(GATKSAMRecord.REDUCED_READ_QUALITY_TAG, convertBaseCounts());
        return samRecord;
    }

    public int size () {
        return bases.size();
    }

    private byte [] convertBaseQualities() {
        return listToByteArray(quals);
    }

    private byte [] convertBaseCounts() {
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
                    System.out.println("BUG CATCHER: " + readName + " " + contig + ":" + refStart + "-" + refStart + bases.size());
                    throw new ReviewedStingException("Trying to create a deletion in the consensus");
                case I:
                    throw new ReviewedStingException("Trying to create an insertion in the consensus");
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
