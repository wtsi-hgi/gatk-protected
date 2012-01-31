package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ReduceReadsIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String BAM = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String DELETION_BAM = validationDataLocation + "filtered_deletion_for_reduce_reads.bam";
    final String L = " -L 20:10,100,000-10,120,000 ";

    private void RRTest(String testName, String args, String md5) {
        String base = String.format("-T ReduceReads -R %s -I %s ", REF, BAM) + " -o %s ";
        WalkerTestSpec spec = new WalkerTestSpec(base + args, Arrays.asList(md5));
        executeTest(testName, spec);
    }

    @Test(enabled = true)
    public void testDefaultCompression() {
        RRTest("testDefaultCompression ", L, "faa8c118de4488e4397884b35f016660");
    }

    @Test(enabled = true)
    public void testMultipleIntervals() {
        String intervals = "-L 20:10,100,000-10,100,500 -L 20:10,200,000-10,200,500 -L 20:10,300,000-10,300,500 -L 20:10,400,000-10,500,000 -L 20:10,500,050-10,500,060 -L 20:10,600,000-10,600,015 -L 20:10,700,000-10,700,110";
        RRTest("testMultipleIntervals ", intervals, "30618b6d85895d16c455d8622637b4a0");
    }

    @Test(enabled = true)
    public void testHighCompression() {
        RRTest("testHighCompression ", " -csmm 10 -minvar 0.3 -mindel 0.3 " + L, "cc094b995e5729165250da43482fdad9");
    }

    @Test(enabled = true)
    public void testLowCompression() {
        RRTest("testLowCompression ", " -csmm 30 -minvar 0.01 -mindel 0.01 -minmap 5 -minqual 5 " + L, "abe92f675e5a8bd68a1bc2062489387e");
    }

    @Test(enabled = true)
    public void testIndelCompression() {
        RRTest("testIndelCompression ", " -csindel 50 -L 20:10,100,500-10,100,600 ", "f23e79cb5c5eb908da9f70fea5545cc4");
    }

    @Test(enabled = true)
    public void testNormalDownsampling() {
        RRTest("testNormalDownsampling ", " -ds 25 -dm Normal " + L, "ae3206c9cd66a22616f6c64d5823e1e1");
    }

    @Test(enabled = true)
    public void testAdaptiveDownsampling() {
        RRTest("testAdaptiveDownsampling ", " -ds 25 -dm Adaptive " + L, "d32a124f329b01168c1c46b8603fa461");
    }
    
    @Test(enabled = true)
    public void testFilteredDeletionCompression() {
        String base = String.format("-T ReduceReads -R %s -I %s ", REF, DELETION_BAM) + " -o %s ";
        executeTest("testFilteredDeletionCompression", new WalkerTestSpec(base, Arrays.asList("122e4e60c4412a31d0aeb3cce879e841")));
    }


}

