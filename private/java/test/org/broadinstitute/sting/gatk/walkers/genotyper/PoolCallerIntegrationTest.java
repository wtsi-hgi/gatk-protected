package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;

import java.util.Arrays;
import org.testng.annotations.Test;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 4/5/12
 * Time: 11:28 AM
 * To change this template use File | Settings | File Templates.
 */
public class PoolCallerIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String CEUTRIO_BAM = "/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37.list";
    final String REFSAMPLE_CALLS = comparisonDataLocation + "Unvalidated/mtDNA/NA12878.snp.vcf";
    final String REFSAMPLE_NAME = "NA12878";
    final String INTERVALS = "MT";
    final String NA12891_CALLS = comparisonDataLocation + "Unvalidated/mtDNA/NA12891.snp.vcf";

    private void PC_MT_Test(String bam, String args, String name, String md5) {
        final String base = String.format("-T PoolCaller -R %s -I %s -L %s -reference %s -refsample %s -glm POOLSNP -ignoreLane -pnrm POOL",
                REF, bam, INTERVALS, REFSAMPLE_CALLS, REFSAMPLE_NAME) + " --no_cmdline_in_header -o %s";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testPoolCaller:"+name+" args=" + args, spec);
    }

    @Test
    public void testMT_SNP_DISCOVERY_sp4() {
        // todo- force maxAlleles = 1 for now since multiallelic exact model not ready
        PC_MT_Test(CEUTRIO_BAM, " -maxAlleles 1 -sp 4", "MT_SNP_DISCOVERY_sp4","0c69e10ef24606737e561640e7092acb");
    }

    @Test
    public void testMT_SNP_GGA_sp10() {

        PC_MT_Test(CEUTRIO_BAM, String.format(" -maxAlleles 1 -sp 10 -gt_mode GENOTYPE_GIVEN_ALLELES -alleles %s",NA12891_CALLS), "MT_SNP_GGA_sp10", "0a6f1c4e12e9fb5d43b1af34717199a0");
    }

}
