package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ReduceReadsIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String BAM = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String L = " -L 20:10,100,000-10,120,000 ";

    private void RRTest(String testName, String args, String md5) {
        String base = String.format("-T ReduceReads -R %s -I %s ", REF, BAM) + " -o %s ";
        WalkerTestSpec spec = new WalkerTestSpec(base + args, Arrays.asList(md5));
        executeTest(testName, spec);
    }

    @Test(enabled = true)
    public void testDefaultCompression() {
        RRTest("testDefaultCompression ", L, "8dfd85511c2431f2cf4fdbd8191f0aa1");
    }

    @Test(enabled = true)
    public void testMultipleIntervals() {
        String intervals = "-L 20:10,100,000-10,100,500 -L 20:10,200,000-10,200,500 -L 20:10,300,000-10,300,500 -L 20:10,400,000-10,500,000 -L 20:10,500,050-10,500,060 -L 20:10,600,000-10,600,015 -L 20:10,700,000-10,700,110";
        RRTest("testMultipleIntervals ", intervals, "691501a71c7d55c045a1ace60755883b");
    }

    @Test(enabled = true)
    public void testHighCompression() {
        RRTest("testHighCompression ", " -csmm 10 -minvar 0.3 -mindel 0.3 " + L, "878a00cdf4a681756979783d678460b9");
    }

    @Test(enabled = true)
    public void testLowCompression() {
        RRTest("testLowCompression ", " -csmm 30 -minvar 0.01 -mindel 0.01 -minmap 5 -minqual 5 " + L, "d3ef6ccc26cb23bda494978e422d640e");
    }

    @Test(enabled = true)
    public void testIndelCompression() {
        RRTest("testIndelCompression ", " -csindel 50 -L 20:10,100,500-10,100,600 ", "bc8873df0f7bfeb539e73321aa9a23b9");
    }

    @Test(enabled = true)
    public void testNormalDownsampling() {
        RRTest("testNormalDownsampling ", " -ds 25 " + L, "d5a45e855a1ba0cf51347ceda541b31b");
    }

    @Test(enabled = true)
    public void testAdaptieDownsampling() {
        RRTest("testAdaptiveDownsampling ", " -ds 25 -dm Adaptive " + L, "42d8f766da6b91a6e0d0623397771f57");
    }


}

