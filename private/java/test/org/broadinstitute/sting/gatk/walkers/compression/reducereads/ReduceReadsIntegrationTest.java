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
        RRTest("testDefaultCompression ", L, "3257cb50ac02a2f453fc69a581482944");
    }

    @Test(enabled = true)
    public void testMultipleIntervals() {
        String intervals = "-L 20:10,100,000-10,100,500 -L 20:10,200,000-10,200,500 -L 20:10,300,000-10,300,500 -L 20:10,400,000-10,500,000 -L 20:10,500,050-10,500,060 -L 20:10,600,000-10,600,015 -L 20:10,700,000-10,700,110";
        RRTest("testMultipleIntervals ", intervals, "20af948a7f2d4537bc03893c046051b5");
    }

    @Test(enabled = true)
    public void testHighCompression() {
        RRTest("testHighCompression ", " -csmm 10 -minvar 0.3 -mindel 0.3 " + L, "47c0643dcfdef279c54c22d48e0206c6");
    }

    @Test(enabled = true)
    public void testLowCompression() {
        RRTest("testLowCompression ", " -csmm 30 -minvar 0.01 -mindel 0.01 -minmap 5 -minqual 5 " + L, "26f5093226fba270bcbb9eb7a06c843e");
    }

    @Test(enabled = true)
    public void testIndelCompression() {
        RRTest("testIndelCompression ", " -csindel 50 -L 20:10,100,500-10,100,600 ", "096c73163dfdcd748d5549e853f78040");
    }

    @Test(enabled = true)
    public void testNormalDownsampling() {
        RRTest("testNormalDownsampling ", " -ds 25 -dm Normal " + L, "0770bc55ac9107e27751ffe51ff4f37f");
    }

    @Test(enabled = true)
    public void testAdaptieDownsampling() {
        RRTest("testAdaptiveDownsampling ", " -ds 25 -dm Adaptive " + L, "22e8e85e42d727599f8f3c32aa95ac10");
    }


}

