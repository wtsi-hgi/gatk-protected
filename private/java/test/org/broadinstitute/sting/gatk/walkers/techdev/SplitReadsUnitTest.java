package org.broadinstitute.sting.gatk.walkers.techdev;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

/**
 * User: carneiro
 * Date: 1/23/13
 * Time: 7:09 PM
 */
public class SplitReadsUnitTest {

    @Test (enabled = true)
    public void splitMultipleTimes() {
        GATKSAMRecord read = createRandomRead(4096);
        for (int splits = 10; splits > 0; splits--) {
            assertReadSplit(read, splits);
        }
    }

    @Test (enabled = true)
    public void splitWithOddLengths() {
        GATKSAMRecord read = createRandomRead(5439);
        for (int splits = 10; splits > 0; splits--) {
            assertReadSplit(read, splits);
        }
    }

    @Test (enabled = true)
    public void chopRead() {
        final int readLength = 500;
        GATKSAMRecord read = createRandomRead(readLength);
        for (int i = 0; i < readLength; i += 10) {
            LinkedList<GATKSAMRecord> result = SplitReads.splitRead(read, 0, i);
            Assert.assertEquals(result.size(), 1);
            GATKSAMRecord choppedRead = result.getFirst();
            Assert.assertEquals(choppedRead.getReadLength(), i);
            assertBases(choppedRead.getReadBases(), read.getReadBases(), 0);
        }
        LinkedList<GATKSAMRecord> result = SplitReads.splitRead(read, 0, 500);
        Assert.assertTrue(result.isEmpty());

    }

    private void assertReadSplit(GATKSAMRecord read, int nSplits) {
        final LinkedList<GATKSAMRecord> splittedReads = SplitReads.splitRead(read, nSplits, -1);

        final int originalLength = read.getReadLength();

        final int expectedReads = (int) Math.pow(2, nSplits);
        Assert.assertEquals(splittedReads.size(), expectedReads);

        final int expectedReadLength = originalLength/expectedReads;
        final boolean hasOddLengths = originalLength % expectedReads != 0;

        int readNumber = 0;
        HashSet<String> readNames = new HashSet<String>(expectedReads*2);
        for (GATKSAMRecord splitRead : splittedReads) {
            if (hasOddLengths) {
                Assert.assertTrue(splitRead.getReadLength() == expectedReadLength || splitRead.getReadLength() == expectedReadLength + 1);
            } else {
                Assert.assertEquals(splitRead.getReadLength(), expectedReadLength);
            }
            assertBases(splitRead.getReadBases(), read.getReadBases(), readNumber * expectedReadLength);
            readNames.add(splitRead.getReadName());
        }
        Assert.assertEquals(readNames.size(), expectedReads);
    }

    private void assertBases(byte[] actualBase, byte[] expectedBase, int startIndex) {
        for (int i = 0; i < actualBase.length; i++) {
            Assert.assertEquals(actualBase[i], expectedBase[startIndex + i]);
        }
    }

    private static GATKSAMRecord createRandomRead(int length) {
        List<CigarElement> cigarElements = new LinkedList<CigarElement>();
        cigarElements.add(new CigarElement(length, CigarOperator.M));
        Cigar cigar = new Cigar(cigarElements);
        return ArtificialSAMUtils.createArtificialRead(cigar);
    }
}
