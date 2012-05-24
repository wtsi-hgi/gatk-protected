package org.broadinstitute.sting.gatk.walkers.validation;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

/**
 * tests for the ROD system in general; from rod system validation to empty VCF files
 */
public class RodSystemValidationIntegrationTest extends WalkerTest {

    public static String baseTestString1KG() {
            return "-T RodSystemValidation -o %s -R " + b36KGReference;
        }

    @Test
    public void testSimpleVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " --eval:VCF3 " + testDir + "MultiSample.vcf", 1,
                Arrays.asList("4541e694cc735fb822b018bb9de25369"));
        executeTest("testSimpleVCFPileup", spec);
    }

    @Test
    public void testEmptyVCF() {
        File vcf = new File(testDir + "justHeader.vcf.idx");
        if (vcf.exists()) vcf.delete();

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " --eval:VCF3 " + testDir + "justHeader.vcf", 1,
                Arrays.asList("1660f76ae84e6e39ec1dfea96622cf5a"));
        executeTest("testEmptyVCF", spec);
    }


    @Test
    public void testComplexVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " --eval:VCF3 " + testDir + "MultiSample.vcf" +
                " --eval:VCF " + validationDataLocation + "NA12878.chr1_10mb_11mb.slx.indels.vcf4"
                , 1,
                Arrays.asList("777f4da740b1acece6966475f3b0c63a"));
        executeTest("testComplexVCFPileup", spec);
    }

    @Test
    public void testLargeComplexVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " --eval:VCF3 " + testDir + "MultiSample.vcf" +
                " --eval:VCF3 " + validationDataLocation + "CEU_hapmap_nogt_23.vcf" +
                " --eval:VCF3 " + validationDataLocation + "CEU_hapmap_nogt_23.vcf" +
                " -L 1 -L 2 -L 20"
                , 1,
                Arrays.asList("427f378b312dad056ac88369887b7d98"));
        executeTest("testLargeComplexVCFPileup", spec);
    }

    //@Test
    public void testBlockZippedVrsUnzippedVCF1() {
        final String vcfName = validationDataLocation + "bgzipped_vcfs/vcfexample.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " --eval:VCF " + vcfName +
                " --eval:VCF3 " + vcfName + ".gz" +
                " --PerLocusEqual"
                , 1,
                Arrays.asList("ab3da32eae65e8c15a9f4a787a190a37"));
        executeTest("testLargeComplexVCFPileup", spec);
    }
}
