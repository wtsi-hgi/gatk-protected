package org.broadinstitute.sting.gatk.walkers.techdev;

import net.sf.samtools.SAMFileWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.LinkedList;

/**
 * Splits and chops reads to simulate a sequencing of smaller read lengths.
 *
 * User: carneiro
 * Date: 1/23/13
 * Time: 6:10 PM
 */

public class SplitReads extends ReadWalker<Integer, Integer> {
    @Argument (shortName = "split", required = false, doc = "how many times to split the read (and all subsequent reads) -- Warning: this option overrides chop")
    private int splits = 0;

    @Argument (shortName = "chop", required = false, doc = "chops the read at the given index. Negative values means no chopping")
    private int chopIndex = -1;

    @Output(doc = "output bam file")
    SAMFileWriter out;

    @Override
    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
        LinkedList<GATKSAMRecord> reads = splitRead(read, splits, chopIndex);
        for (GATKSAMRecord splitRead : reads)
                out.addAlignment(splitRead);
        return reads.size();
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    protected static LinkedList<GATKSAMRecord> splitRead(GATKSAMRecord read, int splits, int chopIndex) {
        LinkedList<GATKSAMRecord> reads = new LinkedList<GATKSAMRecord>();
        if (splits > 0) {
            reads.add(read);
            reads = splitReads(reads, splits);

            // fix mate alignment start information
            final int originalAlignmentStart = read.getAlignmentStart();
            for (GATKSAMRecord splitRead : reads) {
                splitRead.setMateAlignmentStart(splitRead.getMateAlignmentStart() + splitRead.getAlignmentStart() - originalAlignmentStart);
            }
        } else if (chopIndex >= 0 ) {
            if (chopIndex < read.getReadLength()) {
                reads.add(ReadClipper.hardClipByReadCoordinates(read, chopIndex, read.getReadLength() - 1));
            }
        } else {
            reads.add(read);
        }
        return reads;
    }

    private static LinkedList<GATKSAMRecord> splitReads(LinkedList<GATKSAMRecord> reads, int splits) {
        LinkedList<GATKSAMRecord> result = reads;
        if (splits > 0) {
            result = splitReads(reads, splits-1);
            LinkedList<GATKSAMRecord> tempResult = new LinkedList<GATKSAMRecord>();
            for (GATKSAMRecord read : result) {
                tempResult.addAll(0, splitRead(read));
            }
            result = tempResult;
        }

        return result;
    }

    private static LinkedList<GATKSAMRecord> splitRead(GATKSAMRecord read) {
        LinkedList<GATKSAMRecord> result = new LinkedList<GATKSAMRecord>();
        int readLength = read.getReadLength();
        if (readLength == 1) {
            result.addFirst(read);
        } else {
            int splitPoint = readLength / 2;
            GATKSAMRecord splitLeft = ReadClipper.hardClipByReadCoordinates(read, 0, splitPoint-1);
            GATKSAMRecord splitRight = ReadClipper.hardClipByReadCoordinates(read, splitPoint, readLength - 1);
            splitLeft.setReadName(splitLeft.getReadName() + "L");
            splitRight.setReadName(splitRight.getReadName() + "R");
            result.addFirst(splitLeft);
            result.addFirst(splitRight);
        }
        return result;
    }
}
