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
                baseTestString1KG() + " --eval:VCF3 " + privateTestDir + "MultiSample.vcf", 1,
                Arrays.asList("96e2fc4987cfc2475a3284d6204ebe9e"));
        executeTest("testSimpleVCFPileup", spec);
    }

    @Test
    public void testEmptyVCF() {
        File vcf = new File(privateTestDir + "justHeader.vcf.idx");
        if (vcf.exists()) vcf.delete();

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " --eval:VCF3 " + privateTestDir + "justHeader.vcf", 1,
                Arrays.asList("ff7960ddcd22b4ef935162e3dfadd3df"));
        executeTest("testEmptyVCF", spec);
    }


    @Test
    public void testComplexVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " --eval:VCF3 " + privateTestDir + "MultiSample.vcf" +
                " --eval:VCF " + validationDataLocation + "NA12878.chr1_10mb_11mb.slx.indels.vcf4"
                , 1,
                Arrays.asList("7012d9ee9c4f87b6cd3d60eca87640f7"));
        executeTest("testComplexVCFPileup", spec);
    }

    @Test
    public void testLargeComplexVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " --eval:VCF3 " + privateTestDir + "MultiSample.vcf" +
                " --eval:VCF3 " + validationDataLocation + "CEU_hapmap_nogt_23.vcf" +
                " --eval:VCF3 " + validationDataLocation + "CEU_hapmap_nogt_23.vcf" +
                " -L 1 -L 2 -L 20"
                , 1,
                Arrays.asList("96dadb0b54f2797adb13bf9bcae548cc"));
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
