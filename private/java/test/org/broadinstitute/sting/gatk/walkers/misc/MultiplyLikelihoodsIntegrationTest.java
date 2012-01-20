package org.broadinstitute.sting.gatk.walkers.misc;

import org.broadinstitute.sting.WalkerTest;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/20/12
 * Time: 1:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class MultiplyLikelihoodsIntegrationTest extends WalkerTest {

    public static final String baseTestString = " -T MultiplyLikelihoods -o %s -R "+b37KGReference;

    @Test
    public void testLikelihoodsAddedAndRenormalized() {
        WalkerTestSpec spec = new WalkerTestSpec(baseTestString +
        " -V:ex,vcf "+validationDataLocation+"multLik.test.exome.vcf"+
        " -V:chip,vcf "+validationDataLocation+"multLik.test.chip.vcf"+
        " -V:wg,vcf "+validationDataLocation+"multLik.test.genome.vcf",
                1, Arrays.asList("3c59956a3252cbad246ceee2dc7a40f6"));
    }
}
