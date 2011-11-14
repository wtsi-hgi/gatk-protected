package org.broadinstitute.sting.gatk.walkers.reducereads;

import com.google.java.contract.Requires;
import net.sf.samtools.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
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
    private List<BaseIndex> bases;
    private List<Byte> counts;
    private List<Byte> quals;
    private double mappingQuality;          // the average of the rms of the mapping qualities of all the reads that contributed to this consensus
    private String readTag;

    // Information to produce a GATKSAMRecord
    private SAMFileHeader header;
    private GATKSAMReadGroupRecord readGroupRecord;
    private String contig;
    private int contigIndex;
    private String readName;
    private Integer refStart;

    /**
     * Full initialization of the running consensus if you have all the information and are ready to
     * start adding to the running consensus.
     *
     * @param header GATKSAMRecord file header
     * @param readGroupRecord Read Group for the GATKSAMRecord
     * @param contig the read's contig name
     * @param contigIndex the read's contig index
     * @param readName the read's name
     * @param refStart the alignment start (reference based)
     * @param readTag
     */
    public SyntheticRead(SAMFileHeader header, GATKSAMReadGroupRecord readGroupRecord, String contig, int contigIndex, String readName, Integer refStart, String readTag) {
        bases = new LinkedList<BaseIndex>();
        counts = new LinkedList<Byte>();
        quals = new LinkedList<Byte>();
        mappingQuality = 0.0;

        this.readTag = readTag;
        this.header = header;
        this.readGroupRecord = readGroupRecord;
        this.contig = contig;
        this.contigIndex = contigIndex;
        this.readName = readName;
        this.refStart = refStart;
    }

    /**
     * Easy access to keep adding to a running consensus that has already been
     * initialized with the correct read name and refStart
     *
     * @param base
     * @param count
     */
    @Requires("count < Byte.MAX_VALUE")
    public void add(BaseIndex base, byte count, byte qual, double mappingQuality) {
        counts.add(count);
        bases.add(base);
        quals.add(qual);
        this.mappingQuality += mappingQuality;
    }

    public GATKSAMRecord close () {
        GATKSAMRecord read = new GATKSAMRecord(header);
        read.setReferenceName(contig);
        read.setReferenceIndex(contigIndex);
        read.setReadPairedFlag(false);
        read.setReadUnmappedFlag(false);
        read.setAlignmentStart(refStart);
        read.setCigar(buildCigar());
        read.setReadName(readName);
        read.setBaseQualities(convertBaseQualities());
        read.setReadBases(convertReadBases());
        read.setMappingQuality((int) Math.ceil(mappingQuality / bases.size()));
        read.setReadGroup(readGroupRecord);
        read.setAttribute(readTag, convertBaseCounts());
        return read;
    }

    public int size () {
        return bases.size();
    }


    private byte [] convertBaseQualities() {
        return convertVariableGivenBases(bases, quals);
    }

    private byte [] convertBaseCounts() {
        return convertVariableGivenBases(bases, counts);
    }

    private byte [] convertReadBases() {
        byte [] readArray = new byte[getReadLengthWithNoDeletions(bases)];
        int i = 0;
        for (BaseIndex baseIndex : bases)
            if (baseIndex != BaseIndex.D)
                readArray[i++] = baseIndex.getByte();

        return readArray;
    }

    private Cigar buildCigar() {
        LinkedList<CigarElement> cigarElements = new LinkedList<CigarElement>();
        CigarOperator cigarOperator = null;
        int length = 0;
        for (BaseIndex b : bases) {
            CigarOperator op;
            switch (b) {
                case D:
                    op = CigarOperator.DELETION;
                    break;
                case I:
                    throw new ReviewedStingException("Trying to create an insertion in a synthetic read. This operation is currently unsupported.");
                default:
                    op = CigarOperator.MATCH_OR_MISMATCH;
                    break;
            }
            if (cigarOperator == null)
                cigarOperator = op;

            else if (cigarOperator != op) {                                 // if this is a new operator, we need to close the previous one
                cigarElements.add(new CigarElement(length, cigarOperator)); // close previous operator
                cigarOperator = op;
                length = 0;
            }

            length++;                                                       // add this element to the length, either if we just created a new cigarElement, or if it has the same operator as before
        }
        if (length > 0)
            cigarElements.add(new CigarElement(length, cigarOperator));
        return new Cigar(cigarElements);
    }

    /**
     * Shared functionality for all conversion utilities
     *
     * @param bases
     * @param variable
     * @return a converted variable given the bases and skipping deletions
     */

    private static byte [] convertVariableGivenBases (List<BaseIndex> bases, List<Byte> variable) {
        byte [] variableArray = new byte[getReadLengthWithNoDeletions(bases)];
        int i = 0;
        Iterator<Byte> variableIterator = variable.iterator();
        for (BaseIndex baseIndex : bases) {
            byte count = variableIterator.next();
            if (baseIndex != BaseIndex.D)
                variableArray[i++] = count;
        }
        return variableArray;

    }

    /**
     * Shared functionality for all conversion utilities
     *
     * @param bases
     * @return the length of the read with no deletions
     */
    private static int getReadLengthWithNoDeletions(List<BaseIndex> bases) {
        int readLength = bases.size();
        for (BaseIndex baseIndex : bases)
            if (baseIndex == BaseIndex.D)
                readLength--;
        return readLength;
    }


}
