package org.broadinstitute.sting.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * check that we're getting the expected results from the RODs for reads for a variety of input types
 */
public class ValidateRODForReadsIntegrationTest extends WalkerTest {

    private final String vcfFile = privateTestDir + "rodForReadsVCFCheck.withRG.vcf";

    public static String baseTestString() {
        return "-T ValidateRODForReads -o %s -R " + publicTestDir + "exampleFASTA.fasta" + " -I " + publicTestDir + "exampleBAM.bam";
    }

    @Test
         public void testSimpleVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -V " + vcfFile, 1,
                Arrays.asList("f7919e9dc156fb5d3ad0541666864ea5"));
        executeTestParallel("testSimpleVCFPileup", spec, ParallelTestType.NANO_SCHEDULED);
    }

    @Test
    public void testComplexRODs() {
        WalkerTestSpec spec = new WalkerTestSpec(
            "-T ValidateRODForReads -o %s -R "
                    + b37KGReference + " -I "
                    + validationDataLocation + "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam"
                    + " -L 20:10,000,000-11,000,000 "
                    + " -V " + validationDataLocation + "NA12878.HiSeq.WGS.b37_decoy.indel.recalibrated.vcf"
                    + " -V " + validationDataLocation + "NA12878.omni.vcf",
                1,
                Arrays.asList("8060cde53b9de9032ca54f738d2e7d19"));
        executeTestParallel("testComplexRODs", spec, ParallelTestType.NANO_SCHEDULED);
    }
}
