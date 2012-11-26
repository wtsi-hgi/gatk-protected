package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class PolymorphicStatusUnitTest extends BaseTest {
    @DataProvider(name = "PSTest")
    public Object[][] makeGLsWithNonInformative() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final PolymorphicStatus x : PolymorphicStatus.values() )
            tests.add(new Object[]{x, PolymorphicStatus.UNKNOWN, x});

        for ( final PolymorphicStatus x : PolymorphicStatus.values() )
            tests.add(new Object[]{x, x, x});

        tests.add(new Object[]{PolymorphicStatus.POLYMORPHIC, PolymorphicStatus.MONOMORPHIC, PolymorphicStatus.DISCORDANT});
        tests.add(new Object[]{PolymorphicStatus.MONOMORPHIC, PolymorphicStatus.POLYMORPHIC, PolymorphicStatus.DISCORDANT});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PSTest")
    public void testMakeConsensus(final PolymorphicStatus ps1, final PolymorphicStatus ps2, final PolymorphicStatus expected) {
        Assert.assertEquals(ps1.makeConsensus(ps2), expected, "Polymorphic consensus of " + ps1 + " + " + ps2 + " was not expected " + expected);
    }
}