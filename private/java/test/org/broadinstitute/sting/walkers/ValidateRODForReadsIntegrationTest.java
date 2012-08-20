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
        executeTest("testSimpleVCFPileup", spec);
    }
}
