/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

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