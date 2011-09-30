package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.*;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 9/28/11
 * Time: 7:55 PM
 */
public class SlidingReadUnitTest extends BaseTest {

    // TODO: Add error messages on failed tests

    SAMRecord read, expected;
    SlidingRead slidingRead;
    final static String BASES = "ACTG";
    final static String QUALS = "!+5?";

    @BeforeClass
    public void init() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        read = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 1, BASES.length());
        read.setReadUnmappedFlag(true);
        read.setReadBases(new String(BASES).getBytes());
        read.setBaseQualityString(new String(QUALS));

        slidingRead = new SlidingRead(read);
    }

    @Test ( enabled = false )
    public void testClipStart() {
        logger.warn("Executing testClipStart");

        // Testing situations, before read, at start, after end
        Assert.assertEquals(slidingRead.clipStart(0), slidingRead);
        Assert.assertEquals(slidingRead.clipStart(1), slidingRead);
        Assert.assertEquals(slidingRead.clipStart(5), null);

        //clip 1 base
        expected = slidingRead.clipStart(2).getRead();
        Assert.assertEquals(expected.getReadBases(), BASES.substring(1,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(1,4));
        Assert.assertEquals(expected.getCigarString(), "1H3M");

        //clip 2 bases
        expected = slidingRead.clipStart(3).getRead();
        Assert.assertEquals(expected.getReadBases(), BASES.substring(2,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(2,4));
        Assert.assertEquals(expected.getCigarString(), "2H2M");
    }

    @Test ( enabled = false )
    public void testTrimToVariableRegion() {
        logger.warn("Executing testTrimToVariableRegion");

        // Read is fully contained
        Assert.assertEquals(slidingRead.trimToVariableRegion(0,5), slidingRead.getRead());

        // Read spans the entire Region
        Assert.assertEquals(slidingRead.trimToVariableRegion(1,4), slidingRead.getRead());

        // Region is entirely before read
        Assert.assertEquals(slidingRead.trimToVariableRegion(0,0), new SAMRecord( new SAMFileHeader() ));

        // Read is entirely after read
        Assert.assertEquals(slidingRead.trimToVariableRegion(5,6), new SAMRecord( new SAMFileHeader() ));

        // Read overlaps region on the left
        expected = slidingRead.trimToVariableRegion(0,2);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(0,2).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(0,2));
        Assert.assertEquals(expected.getCigarString(), "2M2H");

        // Read overlaps region on the right
        expected = slidingRead.trimToVariableRegion(3,5);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(2,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(2,4));
        Assert.assertEquals(expected.getCigarString(), "2H2M");

        // Read overlaps region on both left and right
        expected = slidingRead.trimToVariableRegion(2,3);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(1,3).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(1,3));
        Assert.assertEquals(expected.getCigarString(), "1H2M1H");

    }
}
