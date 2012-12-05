package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class TruthStatusUnitTest extends BaseTest {
    @DataProvider(name = "TSTest")
    public Object[][] makeGLsWithNonInformative() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final TruthStatus x : TruthStatus.values() )
            tests.add(new Object[]{x, TruthStatus.UNKNOWN, x});

        for ( final TruthStatus x : TruthStatus.values() )
            tests.add(new Object[]{x, x, x});

        tests.add(new Object[]{TruthStatus.FALSE_POSITIVE, TruthStatus.TRUE_POSITIVE, TruthStatus.FALSE_POSITIVE});
        tests.add(new Object[]{TruthStatus.TRUE_POSITIVE, TruthStatus.FALSE_POSITIVE, TruthStatus.FALSE_POSITIVE});

        tests.add(new Object[]{TruthStatus.FALSE_POSITIVE, TruthStatus.SUSPECT, TruthStatus.FALSE_POSITIVE});
        tests.add(new Object[]{TruthStatus.TRUE_POSITIVE, TruthStatus.SUSPECT, TruthStatus.SUSPECT});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TSTest")
    public void testMakeConsensus(final TruthStatus ps1, final TruthStatus ps2, final TruthStatus expected) {
        Assert.assertEquals(ps1.makeConsensus(ps2), expected, "Truth status consensus of " + ps1 + " + " + ps2 + " was not expected " + expected);
    }
}