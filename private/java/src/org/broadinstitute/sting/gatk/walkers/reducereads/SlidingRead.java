package org.broadinstitute.sting.gatk.walkers.reducereads;

import com.google.java.contract.Requires;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;

import java.util.BitSet;

public class SlidingRead {
    protected SAMRecord read;
    protected BitSet baseIsMismatch;

    public SlidingRead(SAMRecord read, ReferenceContext referenceContext) {
        this.read = read;
        baseIsMismatch = AlignmentUtils.mismatchesInRefWindow(read, referenceContext);
        debugRead();
    }

    private SlidingRead(SAMRecord read, BitSet baseIsMismatch) {
        this.read = read;
        this.baseIsMismatch = baseIsMismatch;
    }

    public static SlidingRead createDummySlidingRead(SAMRecord read) {
        BitSet mismatches = new BitSet(read.getReadLength());
        mismatches.set(0, read.getReadLength());
        return new SlidingRead(read, mismatches);
    }

    public SAMRecord getRead() {
        return read;
    }

    public boolean getBaseIsMismatch(int i) {
        return baseIsMismatch.get(i);
    }

    public int getReadLength() {
        return read.getReadLength();
    }

    /**
     * Outputs a read trimmed to the [refStart, refStop] interval (inclusive) with '=' for all bases that match
     * the reference.
     *
     * @param refStart start of the read in reference coordinates (inclusive)
     * @param refStop end of the read in reference coordinates (inclusive)
     * @return a trimmed read with matches ('=') and mismatches.
     */
    public SAMRecord trimToVariableRegion(int refStart, int refStop) {
        int start = read.getAlignmentStart();
        int stop = read.getAlignmentEnd();

        ReadClipper clipper = new ReadClipper(read);

        // check to see if read is contained in region
        if ( start <= refStop && stop >= refStart) {
            if ( start < refStart && stop > refStop )
                return clipper.hardClipBothEndsByReferenceCoordinates(refStart-1, refStop+1);
            else if ( start < refStart )
                return(clipper.hardClipByReferenceCoordinatesLeftTail(refStart-1));
            else if ( stop > refStop )
                return (clipper.hardClipByReferenceCoordinatesRightTail(refStop+1));
            return read;
        }
        else
            return new SAMRecord(read.getHeader());
    }

    public SlidingRead clipStart(int refNewStart) {
        if (refNewStart > read.getAlignmentEnd())
            return null;
        if (refNewStart <= read.getAlignmentStart())
            return this;

        ReadClipper readClipper = new ReadClipper(read);
        SAMRecord clippedRead = readClipper.hardClipByReferenceCoordinatesLeftTail(refNewStart-1);
        int newStart = read.getReadLength() - clippedRead.getReadLength();
        BitSet mismatches = baseIsMismatch.get(newStart, read.getReadLength());

        return new SlidingRead(clippedRead, mismatches);
    }

    /**
     * @return Returns the read with '=' where it matches the reference and the base where it doesn't
     */
    public void compress() {
        byte [] compressedReadBases = new byte[read.getReadLength()];
        byte [] originalReadBases = read.getReadBases();

        for (int i = 0; i < read.getReadLength(); i++) {
            if (baseIsMismatch.get(i))
                compressedReadBases[i] = originalReadBases[i];
            else
                compressedReadBases[i] = BaseIndex.EQ.getByte();
        }
        read.setReadBases(compressedReadBases);
    }

    private static void debugRead (SAMRecord read, BitSet bitSet) {
        byte [] bases = read.getReadBases();
        System.out.print("Read: ");
        for (int i=0; i < read.getReadLength(); i++)
            System.out.print(String.format("%c", bases[i]));
        System.out.println();
        System.out.print("Byte: ");
        for (int i=0; i < read.getReadLength(); i++) {
            if (bitSet.get(i))
                System.out.print("1");
            else
                System.out.print("0");
        }
        System.out.println(" - " + bitSet.size());
    }

    private void debugRead() {
        debugRead(read, baseIsMismatch);
    }

}

