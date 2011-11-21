package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

public class SyntheticReadUnitTest extends BaseTest {
    final SAMFileHeader artificialSAMHeader = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1);
    final GATKSAMReadGroupRecord artificialGATKRG = new GATKSAMReadGroupRecord("synthetic");
    final String artificialContig = "1";
    final int artificialContigIndex = 0;
    final String artificialReadName = "synth";
    final int artificialRefStart = 1;
    final double artificialMappingQuality = 60;
    final Random random = new Random(8854875);


@Test
public void testBaseCounts() {
    for (int i=0; i<10; i++) {
        SyntheticRead syntheticRead = new SyntheticRead(artificialSAMHeader, artificialGATKRG, artificialContig, artificialContigIndex, artificialReadName, artificialRefStart, GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG);
        ReadInfo readInfo = new ReadInfo();
        for (BaseInfo baseInfo : readInfo)
            syntheticRead.add(baseInfo.getBase(), baseInfo.getCounts(), baseInfo.getQual(), artificialMappingQuality);

        byte [] expectedBaseCounts = readInfo.getCounts();
        byte [] compressedBaseCounts = syntheticRead.convertBaseCounts();

        // This reverts back to the original quals, causing overflows on the byte.
        // this works because we're also overflowing in the very same way in the compression.
        byte firstBase = expectedBaseCounts[0];
        for (int j=1; j<expectedBaseCounts.length; j++)
            compressedBaseCounts[j] += firstBase;

        Assert.assertEquals(compressedBaseCounts, expectedBaseCounts);
    }
}

private class ReadInfo implements Iterable<BaseInfo> {
    private BaseIndex [] bases;
    private byte [] quals;
    private byte [] counts;

    private ReadInfo() {
        final BaseIndex[] values = new BaseIndex []{BaseIndex.A, BaseIndex.C, BaseIndex.G, BaseIndex.T};
        final int size = 1 + random.nextInt(99);

        bases = new BaseIndex[size];
        quals = new byte[size];
        counts = new byte[size];

        for (int i = 0; i < bases.length; i++) {
            bases[i] = values[random.nextInt(values.length)];
            quals[i] = (byte) random.nextInt(64);
            counts[i] = (byte) random.nextInt(Byte.MAX_VALUE);
        }
    }

    public byte[] getCounts() {
        return counts;
    }

    @Override
    public Iterator<BaseInfo> iterator() {
        List<BaseInfo> list = new ArrayList<BaseInfo>();
        for (int i=0; i<bases.length; i++)
            list.add(new BaseInfo(bases[i], quals[i], counts[i]));
        return list.iterator();
    }
}

private class BaseInfo {
    private BaseIndex base;
    private byte qual;
    private byte counts;

    private BaseInfo(BaseIndex base, byte qual, byte counts) {
        this.base = base;
        this.qual = qual;
        this.counts = counts;
    }

    public BaseIndex getBase() {
        return base;
    }


    public byte getQual() {
        return qual;
    }


    public byte getCounts() {
        return counts;
    }

}

}
