package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class SiteSelectorUnitTest extends NA12878KBUnitTestBase {
    final MongoVariantContext mvc19_2 = MongoVariantContext.create("y", "19", 2, "A", "C", true);
    final MongoVariantContext mvc20_1 = MongoVariantContext.create("x", "20", 1, "A", "C", false);
    final MongoVariantContext mvc20_3 = MongoVariantContext.create("x", "20", 3, "A", "C", false);
    final MongoVariantContext mvc20_4 = MongoVariantContext.create("y", "20", 4, "A", "C", false);
    final MongoVariantContext mvc20_5 = MongoVariantContext.create("y", "20", 5, "A", "C", true);

    final List<MongoVariantContext> allMVCs = Arrays.asList(mvc19_2, mvc20_1, mvc20_3, mvc20_4, mvc20_5);

    @BeforeMethod
    public void setup() {
        setupBeforeMethod();

        // note out of order adding
        for ( final MongoVariantContext mvc : allMVCs )
            db.addCall(mvc);
    }

    @AfterMethod
    public void teardown() {
        teardownMethod();
    }

    private class CompareMVCs implements Comparator<MongoVariantContext> {
        @Override
        public int compare(MongoVariantContext o1, MongoVariantContext o2) {
            return o1.getLocation(parser).compareTo(o2.getLocation(parser));
        }
    }

    @DataProvider(name = "SiteSelectorQueryTest")
    public Object[][] makeAt() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new SiteSelector(parser), allMVCs});

        for ( final MongoVariantContext mvc : allMVCs )
            tests.add(new Object[]{new SiteSelector(parser).addInterval(mvc.getChr(), mvc.getStart(), mvc.getStart()), Arrays.asList(mvc)});

        tests.add(new Object[]{new SiteSelector(parser).addInterval("20", 1, 1), Arrays.asList(mvc20_1)});
        tests.add(new Object[]{new SiteSelector(parser).addInterval("20", 1, 3), Arrays.asList(mvc20_1, mvc20_3)});
        tests.add(new Object[]{new SiteSelector(parser).addInterval("20", 1, 4), Arrays.asList(mvc20_1, mvc20_3, mvc20_4)});
        tests.add(new Object[]{new SiteSelector(parser).addInterval("20", 3, 4), Arrays.asList(mvc20_3, mvc20_4)});
        tests.add(new Object[]{new SiteSelector(parser).addInterval("20", 1, 5), Arrays.asList(mvc20_1, mvc20_3, mvc20_4, mvc20_5)});
        tests.add(new Object[]{new SiteSelector(parser).addInterval("20", 1, 100), Arrays.asList(mvc20_1, mvc20_3, mvc20_4, mvc20_5)});

        for ( final List<MongoVariantContext> mvcs : Utils.makePermutations(allMVCs, 3, false) ) {
            Collections.sort(mvcs, new CompareMVCs());
            final SiteSelector selector = new SiteSelector(parser);
            for ( final MongoVariantContext mvc : mvcs ) {
                selector.addInterval(mvc.getChr(), mvc.getStart(), mvc.getStart());
            }
            tests.add(new Object[]{selector, mvcs});
        }

        tests.add(new Object[]{new SiteSelector(parser).onlyReviewed(), Arrays.asList(mvc19_2, mvc20_5)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SiteSelectorQueryTest")
    public void testSiteSelectorQuery(final SiteSelector selector, final List<MongoVariantContext> expectMVCs) {
        final SiteIterator<MongoVariantContext> it = db.getCalls(selector);
        final List<MongoVariantContext> actualMVCs = it.toList();
        Assert.assertEquals(actualMVCs, expectMVCs, "Expected " + expectMVCs);
    }
}