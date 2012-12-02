package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.InvalidRecordsLogError;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.InvalidRecordsRemove;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.InvalidRecordsThrowError;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.MongoVariantContextException;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SiteIteratorUnitTest extends NA12878KBUnitTestBase {
    SiteIterator<MongoVariantContext> it;

    @BeforeMethod
    public void setup() {
        setupBeforeMethod();

        // note out of order adding
        db.addCall(MongoVariantContext.create("x", "20", 1, "A", "C", false));
        db.addCall(MongoVariantContext.create("x", "20", 2, "A", "C", false));
        db.addCall(MongoVariantContext.create("x", "20", 4, "A", "C", false));
        db.addCall(MongoVariantContext.create("x", "20", 3, "A", "C", false));
        db.addCall(MongoVariantContext.create("y", "20", 3, "A", "G", false));
        db.addCall(MongoVariantContext.create("y", "20", 4, "A", "C", false));
        db.addCall(MongoVariantContext.create("y", "20", 5, "A", "C", false));

        // adding duplicate record to ensure that dups are filtered out on the fly
        db.addCall(MongoVariantContext.create("y", "20", 4, "A", "C", false));

        it = db.getCalls();
        it.setErrorHandler(new InvalidRecordsThrowError<MongoVariantContext>());
    }

    @AfterMethod
    public void teardown() {
        teardownMethod();
        it.close();
    }

    @Test(enabled = true)
    public void testBasic() {
        Assert.assertTrue(it.hasNext());
        final List<MongoVariantContext> l = it.toList();
        Assert.assertEquals(l.size(), 7);
        Assert.assertFalse(it.hasNext());
    }

    @Test(enabled = true)
    public void testOrder() {
        int lastStart = -1;
        for ( final MongoVariantContext vc : it ) {
            Assert.assertTrue(vc.getStart() >= lastStart);
            lastStart = vc.getStart();
        }
    }

    @Test(enabled = true)
    public void testNextEquivalents() {
        Assert.assertEquals(it.getNextEquivalents().size(), 1); // only 1 at 1
        Assert.assertEquals(it.getNextEquivalents().size(), 1); // only 1 at 2
        Assert.assertEquals(it.getNextEquivalents().size(), 1); // 2 sites at 3, but not equivalent
        Assert.assertEquals(it.getNextEquivalents().size(), 1); // 2 sites at 3, but not equivalent
        Assert.assertEquals(it.getNextEquivalents().size(), 2); // 2 sites at 4, and they are equivalent
        Assert.assertEquals(it.getNextEquivalents().size(), 1); // 1 sites at 5
        Assert.assertEquals(it.hasNext(), false); // no more variants
    }

    @DataProvider(name = "Before")
    public Object[][] makeBefore() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{1, 0});
        tests.add(new Object[]{2, 1});
        tests.add(new Object[]{3, 2});
        tests.add(new Object[]{4, 4});
        tests.add(new Object[]{5, 6});
        tests.add(new Object[]{6, 7});
        tests.add(new Object[]{8, 7});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "Before")
    public void testBefore(final int startThres, final int expectedCount) {
        final List<MongoVariantContext> l = it.getSitesBefore(parser.createGenomeLoc("20", startThres, startThres));
        Assert.assertEquals(l.size(), expectedCount, "Query returned more results than expected");
        for ( final MongoVariantContext mvc : l )
            Assert.assertTrue(mvc.getStart() <= startThres, "MVC " + mvc + " has start > threshold " + startThres);
    }

    @DataProvider(name = "BeforeIndels")
    public Object[][] makeBeforeIndels() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final MongoVariantContext mvc1 = MongoVariantContext.create("x", "20", 1, "A", "C", false);
        final MongoVariantContext mvc2 = MongoVariantContext.create("x", "20", 2, "ACT", "C", false);
        final MongoVariantContext mvc2_1 = MongoVariantContext.create("x", "20", 2, "ACTGT", "C", false);
        final MongoVariantContext mvc3 = MongoVariantContext.create("x", "20", 3, "A", "G", false);
        final MongoVariantContext mvc4 = MongoVariantContext.create("x", "20", 4, "A", "C", false);
        final List<MongoVariantContext> mvcs = Arrays.asList(mvc1, mvc2, mvc2_1, mvc3, mvc4);

        tests.add(new Object[]{mvcs, 1, 0});
        tests.add(new Object[]{mvcs, 2, 1});
        tests.add(new Object[]{mvcs, 3, 3});
        tests.add(new Object[]{mvcs, 4, 4});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "BeforeIndels")
    public void testBeforeIndels(final List<MongoVariantContext> mvcs, final int startThres, final int expectedCount) {
        db.reset();
        db.addCalls(mvcs);
        it = db.getCalls();
        final List<MongoVariantContext> l = it.getSitesBefore(parser.createGenomeLoc("20", startThres, startThres));
        Assert.assertEquals(l.size(), expectedCount, "Query returned more results than expected");
        for ( final MongoVariantContext mvc : l )
            Assert.assertTrue(mvc.getStart() <= startThres, "MVC " + mvc + " has start > threshold " + startThres);
    }

    @DataProvider(name = "At")
    public Object[][] makeAt() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{1, 1});
        tests.add(new Object[]{2, 1});
        tests.add(new Object[]{3, 2});
        tests.add(new Object[]{4, 2});
        tests.add(new Object[]{5, 1});
        tests.add(new Object[]{6, 0});
        tests.add(new Object[]{8, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true,dataProvider = "At")
    public void testAt(final int start, final int expectedCount) {
        final List<MongoVariantContext> l = it.getSitesAtLocation(parser.createGenomeLoc("20", start, start));
        Assert.assertEquals(l.size(), expectedCount, "Query returned more results than expected");
        for ( final MongoVariantContext mvc : l )
            Assert.assertTrue(mvc.getStart() == start, "MVC " + mvc + " has start != threshold " + start);
    }

    @DataProvider(name = "BasicIteration")
    public Object[][] makeBasicIteration() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final MongoVariantContext mvc1 = MongoVariantContext.create("x", "20", 1, "A", "C", false);
        final MongoVariantContext mvc2 = MongoVariantContext.create("x", "20", 2, "A", "C", false);
        final MongoVariantContext mvc3 = MongoVariantContext.create("x", "20", 2, "A", "G", false);
        final MongoVariantContext mvc4 = MongoVariantContext.create("x", "20", 3, "A", "C", false);

        tests.add(new Object[]{Arrays.asList()});
        tests.add(new Object[]{Arrays.asList(mvc1)});
        tests.add(new Object[]{Arrays.asList(mvc2)});
        tests.add(new Object[]{Arrays.asList(mvc3)});
        tests.add(new Object[]{Arrays.asList(mvc1, mvc2)});
        tests.add(new Object[]{Arrays.asList(mvc1, mvc2, mvc3)});
        tests.add(new Object[]{Arrays.asList(mvc1, mvc2, mvc3, mvc4)});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "BasicIteration")
    public void testBasicIteration(final List<MongoVariantContext> mvcs) {
        db.reset();
        for ( final MongoVariantContext mvc : mvcs ) db.addCall(mvc);
        final SiteIterator<MongoVariantContext> it = db.getCalls();

        for ( int n = 0; n < mvcs.size(); n++ ) {
            Assert.assertTrue(it.hasNext());
            final MongoVariantContext read1 = it.next();
            Assert.assertEquals(read1, mvcs.get(n));
        }

        Assert.assertFalse(it.hasNext());
    }

    @DataProvider(name = "BadMVCs")
    public Object[][] makeBadMVCsProvider() throws CloneNotSupportedException {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final MongoVariantContext bad : MongoVariantContextUnitTest.makeBadMVCs() )
            tests.add(new Object[]{bad});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "BadMVCs", expectedExceptions = MongoVariantContextException.class)
    public void testIteratorWithBadMVCsErrors(final MongoVariantContext mvc) {
        it.close();
        db.addCall(mvc);
        it = db.getCalls();
        it.setErrorHandler(new InvalidRecordsThrowError<MongoVariantContext>());
        it.toList();
    }

    @Test(enabled = true, dataProvider = "BadMVCs")
    public void testIteratorWithBadMVCsLogError(final MongoVariantContext bad) {
        final List<MongoVariantContext> expected = it.toList();
        it.close();
        db.addCall(bad);
        it = db.getCalls();
        final InvalidRecordsLogError<MongoVariantContext> handler = new InvalidRecordsLogError<MongoVariantContext>();
        it.setErrorHandler(handler);
        final List<MongoVariantContext> withoutBad = it.toList();
        Assert.assertEquals(withoutBad, expected);
        Assert.assertEquals(handler.getnBad(), 1);
    }

    @Test(enabled = true, dataProvider = "BadMVCs")
    public void testIteratorWithBadMVCsRemoving(final MongoVariantContext bad) {
        final List<MongoVariantContext> expected = it.toList();
        it.close();
        db.addCall(bad);
        it = db.getCalls();
        final InvalidRecordsRemove<MongoVariantContext> handler = new InvalidRecordsRemove<MongoVariantContext>(db.sites);
        it.setErrorHandler(handler);
        final List<MongoVariantContext> withoutBad = it.toList();
        Assert.assertEquals(withoutBad, expected);
        Assert.assertEquals(handler.getnBad(), 1);

        // now that we've removed the record, we should be able to reread the db without protection and get the right answer
        it = db.getCalls();
        final InvalidRecordsThrowError<MongoVariantContext> errorThrower = new InvalidRecordsThrowError<MongoVariantContext>();
        it.setErrorHandler(errorThrower);
        final List<MongoVariantContext> withoutBadFilter = it.toList();
        Assert.assertEquals(withoutBadFilter, expected);
    }
}