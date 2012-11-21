package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MongoGenotypeUnitTest extends BaseTest {
    final VariantContext vc = new VariantContextBuilder("x", "1", 1, 1, Arrays.asList(Allele.create("C", true), Allele.create("T"))).make();

    @DataProvider(name = "MGBasic")
    public Object[][] makeMGBasic() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( int i = 0; i < 2; i++ ) {
            for ( int j = 0; j < 2; j++ ) {
                for ( final int DP : Arrays.asList(1, 10, 100))
                    for ( final int GQ : Arrays.asList(1, 10, 100))
                        tests.add(new Object[]{vc.getAlleles(), i, j, DP, GQ});
            }
        }

        tests.add(new Object[]{vc.getAlleles(), -1, -1, -1, -1});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MGBasic")
    public void testVCToMVC(final List<Allele> alleles, final int allele1, final int allele2, final int GQ, final int DP) {
        final MongoGenotype mg = new MongoGenotype(allele1, allele2, GQ, DP);
        Assert.assertEquals(mg.getAllele1(), allele1);
        Assert.assertEquals(mg.getAllele2(), allele2);
        Assert.assertEquals(mg.getGQ(), GQ);
        Assert.assertEquals(mg.getDP(), DP);

        mg.setDP(1);
        Assert.assertEquals(mg.getDP(), 1);

        mg.setGQ(1);
        Assert.assertEquals(mg.getGQ(), 1);

        final Genotype gt = mg.toGenotype(alleles);
        Assert.assertEquals(gt.getSampleName(), "NA12878");
        Assert.assertEquals(mg.getAllele1(), alleles.indexOf(gt.getAllele(0)));
        Assert.assertEquals(mg.getAllele2(), alleles.indexOf(gt.getAllele(1)));
        Assert.assertEquals(mg.getDP() != -1, gt.hasDP());
        Assert.assertEquals(mg.getGQ() != -1, gt.hasGQ());

        Assert.assertNotNull(mg.toString());
    }

    @DataProvider(name = "MGPolyStatus")
    public Object[][] makeMGPolyStatus() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new MongoGenotype(-1, -1), PolymorphicStatus.UNKNOWN});
        tests.add(new Object[]{new MongoGenotype(0, 0), PolymorphicStatus.MONOMORPHIC});
        tests.add(new Object[]{new MongoGenotype(0, 1), PolymorphicStatus.POLYMORPHIC});
        tests.add(new Object[]{new MongoGenotype(1, 0), PolymorphicStatus.POLYMORPHIC});
        tests.add(new Object[]{new MongoGenotype(1, 1), PolymorphicStatus.POLYMORPHIC});
        tests.add(new Object[]{new MongoGenotype(1, 1, MongoGenotype.DISCORDANT_GQ, 1), PolymorphicStatus.DISCORDANT});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MGPolyStatus")
    public void testMGPolyStatus(final MongoGenotype mg, final PolymorphicStatus expected) {
        Assert.assertEquals(mg.getPolymorphicStatus(), expected);

        Assert.assertEquals(mg.isMonomorphic(), expected == PolymorphicStatus.MONOMORPHIC);
        Assert.assertEquals(mg.isPolymorphic(), expected == PolymorphicStatus.POLYMORPHIC);
        Assert.assertEquals(mg.isUnknown(), expected == PolymorphicStatus.UNKNOWN);
        Assert.assertEquals(mg.isDiscordant(), expected == PolymorphicStatus.DISCORDANT);

    }
}