package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class HaplotypeCallerIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String BAM = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String INTERVALS_FILE = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.test.intervals";
    final String RECAL_FILE = validationDataLocation + "NA12878.kmer.7.recal_data.csv";

    private void RRTest(String bam, String args, String md5) {
        String base = String.format("-T HaplotypeCaller -R %s -I %s -L %s -recalFile %s", REF, bam, INTERVALS_FILE, RECAL_FILE) + " -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(base + args, Arrays.asList(md5));
        executeTest("testHaplotypeCaller: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerSingleSample() {
        RRTest(BAM, "", "7b2a527f182bb238798eeca74330dd3c");
    }
}

