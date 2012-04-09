package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * @author carneiro
 * @since 3/27/12
 */
public class BaseQualityScoreRecalibratorIntegrationTest extends WalkerTest {
    @Test(enabled = true)
    public void recalibrateTest() {
        String REF = "public/testdata/exampleFASTA.fasta";
        String BAM = "public/testdata/exampleBAM.bam";
        String DBSNP = "public/testdata/exampleDBSNP.vcf";
        String base = String.format("-T BaseQualityScoreRecalibrator -R %s -I %s -knownSites %s", REF, BAM, DBSNP) + " -o %s ";
        WalkerTestSpec spec = new WalkerTestSpec(base, Arrays.asList("7bca63fc8cb810b8682928092ac4344a"));
        executeTest("recalibrateTest", spec);
    }

    @Test(enabled = true)
    public void onTheFlyRecalibrationTest() {
        String REF = "public/testdata/exampleFASTA.fasta";
        String BAM = "public/testdata/exampleBAM.bam";
        String GRP = "public/testdata/exampleGRP.grp";
        String base = String.format("-T PrintReads -R %s -I %s -BQSR %s", REF, BAM, GRP) + " -o %s ";
        WalkerTestSpec spec = new WalkerTestSpec(base, Arrays.asList("058961e5d166005b70f300f20afc60b5"));
        executeTest("onTheFlyRecalibrationTest", spec);
    }

}
