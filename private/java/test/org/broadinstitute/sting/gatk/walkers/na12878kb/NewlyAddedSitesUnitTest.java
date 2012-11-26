package org.broadinstitute.sting.gatk.walkers.na12878kb;

import com.mongodb.DBCursor;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class NewlyAddedSitesUnitTest extends NA12878KBUnitTestBase {
    private static Logger logger = Logger.getLogger(NewlyAddedSitesUnitTest.class);

    private List<MongoVariantContext> makeAllMVCs() {
        final MongoVariantContext mvc19_2 = MongoVariantContext.create("y", "19", 2, "A", "C", true);
        final MongoVariantContext mvc20_1 = MongoVariantContext.create("x", "20", 1, "A", "C", false);
        final MongoVariantContext mvc20_3 = MongoVariantContext.create("x", "20", 3, "A", "C", false);
        final MongoVariantContext mvc20_4 = MongoVariantContext.create("y", "20", 4, "A", "C", false);
        final MongoVariantContext mvc20_4_2 = MongoVariantContext.create("z", "20", 4, "A", "C", false);
        final MongoVariantContext mvc20_5 = MongoVariantContext.create("y", "20", 5, "A", "C", true);
        final MongoVariantContext mvc20_5_2 = MongoVariantContext.create("z", "20", 5, "A", "G", false);

        final List<MongoVariantContext> allMVCs =
                Arrays.asList(mvc19_2,
                        mvc20_1,
                        mvc20_3,
                        mvc20_4,
                        mvc20_4_2,
                        mvc20_5,
                        mvc20_5_2);

        return allMVCs;
    }

    @BeforeMethod
    public void setup() {
        setupBeforeMethod();
    }

    @AfterMethod
    public void teardown() {
        teardownMethod();
    }

    @DataProvider(name = "NewlyAddedTest")
    public Object[][] makeAt() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final List<MongoVariantContext> empty = Collections.emptyList();
        final List<MongoVariantContext> sites = new ArrayList<MongoVariantContext>(makeAllMVCs());
        for ( int i = 0; i < sites.size(); i++ ) {
            final int n = sites.size();
            for ( int split = 0; split <= sites.size(); split++ ) {
                final List<MongoVariantContext> newlyAllocated = new ArrayList<MongoVariantContext>(makeAllMVCs());
                Collections.rotate(newlyAllocated, i);
                final List<MongoVariantContext> befores = newlyAllocated.subList(0, split);
                final List<MongoVariantContext> afters = split == n ? empty : newlyAllocated.subList(split + 1, n);
                tests.add(new Object[]{befores, afters, -1});
            }
        }

        // specific tests that the by location query is working
        tests.add(new Object[]{
                Arrays.asList(MongoVariantContext.create("y", "20", 4, "A", "C", false)),
                Arrays.asList(MongoVariantContext.create("z", "20", 4, "A", "C", false)),
                2});

        tests.add(new Object[]{
                Arrays.asList(MongoVariantContext.create("y", "20", 3, "A", "C", false), MongoVariantContext.create("y", "20", 4, "A", "C", false)),
                Arrays.asList(MongoVariantContext.create("z", "20", 4, "A", "C", false)),
                2});

        tests.add(new Object[]{
                Arrays.asList(MongoVariantContext.create("y", "20", 3, "A", "C", false), MongoVariantContext.create("y", "20", 4, "A", "C", false), MongoVariantContext.create("y", "20", 5, "A", "C", false)),
                Arrays.asList(MongoVariantContext.create("z", "20", 4, "A", "C", false)),
                2});

        tests.add(new Object[]{
                Arrays.asList(MongoVariantContext.create("y", "20", 3, "A", "C", false), MongoVariantContext.create("y", "20", 4, "A", "C", false), MongoVariantContext.create("y", "20", 5, "A", "C", false)),
                Arrays.asList(MongoVariantContext.create("z", "20", 4, "A", "C", false), MongoVariantContext.create("z", "20", 5, "A", "C", false)),
                4});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "NewlyAddedTest")
    public void testQueryForRecords(final List<MongoVariantContext> befores, final List<MongoVariantContext> newlyAdded, int expecteByLocation) {
        db.addCalls(befores);

        final NewlyAddedSites newGetter = new NewlyAddedSites(db);

        db.addCalls(newlyAdded);

        final DBCursor cursor = newGetter.getNewlyAddedRecords();
        Assert.assertEquals(cursor.size(), newlyAdded.size());

        final SiteIterator<MongoVariantContext> it = new SiteIterator<MongoVariantContext>(parser, cursor);
        for ( final MongoVariantContext mvc : it ) {
            Assert.assertTrue(newlyAdded.contains(mvc));
        }
    }

    @Test(enabled = true, dataProvider = "NewlyAddedTest", dependsOnMethods =  "testQueryForRecords")
    public void testQueryForLocations(final List<MongoVariantContext> befores, final List<MongoVariantContext> newlyAdded, int expecteByLocation) {
//        logger.warn("testQueryForLocations " + befores + " " + newlyAdded);
        db.addCalls(befores);

        final int maxToItemize = 100;
        final NewlyAddedSites newGetter = new NewlyAddedSites(db);

        db.addCalls(newlyAdded);

        final SiteSelector selector = newGetter.getNewlyAddedLocations(parser, maxToItemize); // TODO -- fixme

        if ( newlyAdded.isEmpty() ) {
            Assert.assertNull(selector);
        } else {
            Assert.assertNotNull(selector);
            final SiteIterator<MongoVariantContext> it = db.getCalls(selector);
            final List<MongoVariantContext> l = it.toList();

            if ( expecteByLocation != -1 )
                Assert.assertEquals(l.size(), expecteByLocation);

            final List<MongoVariantContext> allAdded = new LinkedList<MongoVariantContext>();
            allAdded.addAll(befores);
            allAdded.addAll(newlyAdded);

            final Collection<MongoVariantContext> expected = new HashSet<MongoVariantContext>();
            for ( final MongoVariantContext possible : allAdded ) {
                final GenomeLoc possibleLoc = possible.getLocation(parser);
                for ( final MongoVariantContext added : newlyAdded ) {
                    if ( possibleLoc.startsAt(added.getLocation(parser)) ) {
                        expected.add(possible);
                    }
                }
            }

//            logger.warn("Added " + newlyAdded.size() + " expected " + expected.size());
            // the only contract is that l must contain at least all of the records in expected, but potentially more
            for ( final MongoVariantContext mvc : expected ) {
//                logger.warn("  Expected mvc " + mvc);
                Assert.assertTrue(l.contains(mvc));
            }
        }
    }

    @DataProvider(name = "NewlyAddedTestManyValues")
    public Object[][] makeNewlyAddedTestManyValues() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final int nInDB : Arrays.asList(2, 10, 20) ) {
            for ( int maxQueries = 1; maxQueries < 2 * nInDB; maxQueries++ ) {
                final List<MongoVariantContext> toAdd = new LinkedList<MongoVariantContext>();
                for ( int j = 1; j <= nInDB; j++ )
                    toAdd.add(MongoVariantContext.create("y", "20", j, "A", "C", false));
                tests.add(new Object[]{toAdd, maxQueries});
                tests.add(new Object[]{toAdd, -1});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "NewlyAddedTestManyValues")
    public void testQueryForLocationsManyResults(final List<MongoVariantContext> toAdd, final int maxQueries) {
        final int n = toAdd.size();
        final List<MongoVariantContext> pre = toAdd.subList(0, n / 2);
        final List<MongoVariantContext> post = toAdd.subList(n / 2, n);
        db.addCalls(pre);
        final NewlyAddedSites newGetter = new NewlyAddedSites(db);
        db.addCalls(post);

        final SiteSelector selector = newGetter.getNewlyAddedLocations(parser, maxQueries);

        final SiteIterator<MongoVariantContext> it = db.getCalls(selector);
        final List<MongoVariantContext> l = it.toList();

        if ( post.size() > maxQueries && maxQueries != -1 ) {
            // we get everything when there are more maxResults to
            Assert.assertEquals(l.size(), toAdd.size());
        } else {
            Assert.assertEquals(l.size(), post.size());
        }
    }

}