package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class HaplotypeCallerIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String NA12878_BAM = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String CEUTRIO_BAM = validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam";
    final String INTERVALS_FILE = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.test.intervals";
    //final String RECAL_FILE = validationDataLocation + "NA12878.kmer.8.subset.recal_data.bqsr";

    private void HCTest(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s -L %s", REF, bam, INTERVALS_FILE) + " --no_cmdline_in_header -o %s --allowMissingVCFHeaders";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        spec.disableShadowBCF();
        //
        // TODO TODO TODO TODO TODO TODO TODO TODO
        // TODO TODO TODO TODO TODO TODO TODO TODO
        //
        // TODO WHEN THE HC EMITS VALID VCF HEADERS ENABLE BCF AND REMOVE ALLOWMISSINGVCFHEADERS ARGUMENTS
        //
        // TODO TODO TODO TODO TODO TODO TODO TODO
        // TODO TODO TODO TODO TODO TODO TODO TODO
        // TODO TODO TODO TODO TODO TODO TODO TODO
        //
        executeTest("testHaplotypeCaller: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerMultiSample() {
        HCTest(CEUTRIO_BAM, "", "e3c10fe7c12485a566377fdec21890c5");
    }

    @Test
    public void testHaplotypeCallerSingleSample() {
        HCTest(NA12878_BAM, "", "5efb574b67aac6480d1f34908d88bd32");
    }

    @Test
    public void testHaplotypeCallerMultiSampleGGA() {
        HCTest(CEUTRIO_BAM, "-gt_mode GENOTYPE_GIVEN_ALLELES -alleles " + validationDataLocation + "combined.phase1.chr20.raw.indels.sites.vcf", "9c1db2ce10973cef5675a0f7f8636f28");
    }
}

