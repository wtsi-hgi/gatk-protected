package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ReduceReadsIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String BAM = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String L = "20:10,100,000-10,120,000";

    private void RRTest(String args, String md5) {
        String base = String.format("-T ReduceReads -R %s -I %s -L %s", REF, BAM, L) + " -o %s ";
        WalkerTestSpec spec = new WalkerTestSpec(base + args, Arrays.asList(md5));
        executeTest("testReduceReads1: args=" + args, spec);
    }

    @Test(enabled = true)
    public void testDefaultCompression() {
        RRTest("", "da5fe1fb4132ce9884ea118f1abda2aa");
    }

    @Test(enabled = true)
    public void testHighCompression() {
        RRTest(" -cs 10 -minvar 0.3 -mindel 0.3", "828e0e4fa5d973252b08cb73faf148fb");
    }

    @Test(enabled = true)
    public void testLowCompression() {
        RRTest(" -cs 30 -minvar 0.01 -mindel 0.01 -minmap 5 -minqual 5", "1cad9627d7d13ce9ad6281caae7170fe");
    }

}

