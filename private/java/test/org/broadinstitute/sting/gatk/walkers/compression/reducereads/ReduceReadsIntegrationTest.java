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
        RRTest("", "0f9072be5e03a47e325f8f4efe0f851f");
    }

    @Test(enabled = true)
    public void testHighCompression() {
        RRTest(" -csmm 10 -minvar 0.3 -mindel 0.3", "a68ccfa29fa037ed126caeacd3388c30");
    }

    @Test(enabled = true)
    public void testLowCompression() {
        RRTest(" -csmm 30 -minvar 0.01 -mindel 0.01 -minmap 5 -minqual 5", "8efe7819d764807718305759c2f9c141");
    }

    @Test(enabled = true)
    public void testIndelCompression() {
        String base = String.format("-T ReduceReads -R %s -I %s -L 20:10,100,500-10,100,600 ", REF, BAM) + " -o %s ";
        WalkerTestSpec spec = new WalkerTestSpec(base, Arrays.asList("bc8873df0f7bfeb539e73321aa9a23b9"));
        executeTest("testIndelCompression ", spec);
    }

}

