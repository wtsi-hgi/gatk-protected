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

package org.broadinstitute.sting.utils.distributedutils;


// the imports for unit testing.


import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.distributedutils.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.*;

/**
 * Basic unit test for GenomeLoc
 */
public class GenomeLocProcessingTrackerUnitTest extends BaseTest {
    IndexedFastaSequenceFile fasta = null;
    GenomeLocParser genomeLocParser = null;
    String chr1 = null;
    private final static String FILE_ROOT = "public/testdata/GLPTFile";

    @BeforeTest
    public void before() {
        File referenceFile = new File(hg18Reference);
        try {
            fasta = new IndexedFastaSequenceFile(referenceFile);
            chr1 = fasta.getSequenceDictionary().getSequence(1).getSequenceName();
            genomeLocParser = new GenomeLocParser(fasta);

        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }
    }

    @BeforeMethod
    public void beforeMethod(Object[] data) {
        if ( data.length > 0 )
            ((TestTarget)data[0]).init();
    }

    @AfterMethod
    public void afterMethod(Object[] data) {
        if ( data.length > 0 ) {
            ((TestTarget)data[0]).getTracker().close();
            ((TestTarget)data[0]).cleanup();
        }
    }

    abstract private class TestTarget {
        String name;
        int nShards;
        int shardSize;
        File file;

        public void init() { cleanup(); }

        public void cleanup() {
            if ( file != null && file.exists() )
                file.delete();
        }

        public boolean isThreadSafe() { return true; }

        protected TestTarget(String name, int nShards, int shardSize, File file) {
            this.name = name;
            this.nShards = nShards;
            this.shardSize = shardSize;
            this.file = file;
        }

        public abstract GenomeLocProcessingTracker getTracker();

        public List<GenomeLoc> getShards() {
            List<GenomeLoc> shards = new ArrayList<GenomeLoc>();
            for ( int i = 0; i < nShards; i++ ) {
                int start = shardSize * i;
                int stop = start + shardSize;
                shards.add(genomeLocParser.createGenomeLoc(chr1, start, stop));
            }
            return shards;
        }

        public String toString() {
            return String.format("TestTarget %s: nShards=%d shardSize=%d", name, nShards, shardSize);
        }
    }

    @DataProvider(name = "threadData")
    public Object[][] createThreadData() {
        // gotta keep the tests small...
        return createData(Arrays.asList(10, 100), Arrays.asList(10));
        //return createData(Arrays.asList(10, 100, 1000, 10000), Arrays.asList(10));
    }

    public Object[][] createData(List<Integer> nShards, List<Integer> shardSizes) {
        List<TestTarget> params = new ArrayList<TestTarget>();

        int counter = 0;
        String name = null;
        for ( int nShard : nShards ) {
            for ( int shardSize : shardSizes ) {
                // shared mem -- canonical implementation
                params.add(new TestTarget("ThreadSafeSharedMemory", nShard, shardSize, null) {
                    GenomeLocProcessingTracker tracker = new SharedMemoryGenomeLocProcessingTracker(new ClosableReentrantLock());
                    public GenomeLocProcessingTracker getTracker() { return tracker; }
                });

                final File file1 = new File(String.format("%s_ThreadSafeFileBacked_%d_%d", FILE_ROOT, counter++, nShard, shardSize));
                params.add(new TestTarget("ThreadSafeFileBacked", nShard, shardSize, file1) {
                    GenomeLocProcessingTracker tracker = new FileBackedGenomeLocProcessingTracker(file1, genomeLocParser, new ClosableReentrantLock(), null);
                    public GenomeLocProcessingTracker getTracker() { return tracker; }
                });

                name = "FileBackedSharedFileThreadSafe";
                final File file2 = new File(String.format("%s_%s_%d_%d", FILE_ROOT, name, counter++, nShard, shardSize));
                params.add(new TestTarget(name, nShard, shardSize, file2) {
                    GenomeLocProcessingTracker tracker = new FileBackedGenomeLocProcessingTracker(file2, genomeLocParser, new SharedFileThreadSafeLock(file2, -1), null);
                    public GenomeLocProcessingTracker getTracker() { return tracker; }
                });

                name = "FileBackedSharedFile";
                final File file3 = new File(String.format("%s_%s_%d_%d", FILE_ROOT, name, counter++, nShard, shardSize));
                params.add(new TestTarget(name, nShard, shardSize, file3) {
                    GenomeLocProcessingTracker tracker = new FileBackedGenomeLocProcessingTracker(file3, genomeLocParser, new SharedFileLock(file3, -1), null);
                    public GenomeLocProcessingTracker getTracker() { return tracker; }
                    public boolean isThreadSafe() { return false; }
                });
            }
        }

        List<Object[]> params2 = new ArrayList<Object[]>();
        for ( TestTarget x : params ) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }

    @DataProvider(name = "simpleData")
    public Object[][] createSimpleData() {
        return createData(Arrays.asList(1000), Arrays.asList(100));
    }

    private static final String NAME_ONE   = "name1";
    private static final String NAME_TWO   = "name2";

    @Test(enabled = true)
    public void testNoop() {
        GenomeLocProcessingTracker tracker = new NoOpGenomeLocProcessingTracker();
        for ( int start = 1; start < 100; start++ ) {
            for ( int n = 0; n < 2; n++ ) {
                GenomeLoc loc = genomeLocParser.createGenomeLoc(chr1, start, start +1);
                ProcessingLoc ploc = tracker.claimOwnership(loc, NAME_ONE);
                Assert.assertTrue(ploc.isOwnedBy(NAME_ONE));
                Assert.assertEquals(tracker.updateAndGetProcessingLocs(NAME_ONE).size(), 0);
            }
        }
    }

    @Test(dataProvider = "simpleData", enabled = true)
    public void testSingleProcessTracker(TestTarget test) {
        GenomeLocProcessingTracker tracker = test.getTracker();
        List<GenomeLoc> shards = test.getShards();
        logger.warn("testSingleProcessTracker " + test);

        int counter = 0;
        for ( GenomeLoc shard : shards ) {
            counter++;

            Assert.assertNull(tracker.findOwner(shard, NAME_ONE));
            Assert.assertFalse(tracker.locIsOwned(shard, NAME_ONE));

            ProcessingLoc proc = tracker.claimOwnership(shard,NAME_ONE);
            Assert.assertNotNull(proc);
            Assert.assertNotNull(proc.getLocation());
            Assert.assertNotNull(proc.getOwner());
            Assert.assertEquals(proc.getLocation(), shard);
            Assert.assertEquals(proc.getOwner(), NAME_ONE);
            Assert.assertEquals(tracker.findOwner(shard, NAME_ONE), proc);
            Assert.assertTrue(tracker.locIsOwned(shard, NAME_ONE));
            Assert.assertNotNull(tracker.updateAndGetProcessingLocs(NAME_ONE));
            Assert.assertEquals(tracker.updateAndGetProcessingLocs(NAME_ONE).size(), counter);

            ProcessingLoc badClaimAttempt = tracker.claimOwnership(shard,NAME_TWO);
            Assert.assertFalse(badClaimAttempt.getOwner().equals(NAME_TWO));
            Assert.assertEquals(badClaimAttempt.getOwner(), NAME_ONE);
        }
    }

    @Test(dataProvider = "simpleData", enabled = true)
    public void testIterator(TestTarget test) {
        GenomeLocProcessingTracker tracker = test.getTracker();
        List<GenomeLoc> shards = test.getShards();
        logger.warn("testIterator " + test);

        List<GenomeLoc> markedShards = new ArrayList<GenomeLoc>();
        List<GenomeLoc> toFind = new ArrayList<GenomeLoc>();

        for ( int i = 0; i < shards.size(); i++ ) {
            if ( ! (i % 10 == 0) ) {
                markedShards.add(shards.get(i));
                tracker.claimOwnership(shards.get(i), NAME_TWO);
            } else {
                toFind.add(shards.get(i));
            }
        }

        int nFound = 0;
        Iterator<GenomeLoc> it = shards.iterator();
        while ( it.hasNext() ) {
            GenomeLoc shard = tracker.claimOwnershipOfNextAvailable(it, NAME_ONE);

            if ( shard == null ) { // everything to get is done
                Assert.assertEquals(nFound, toFind.size(), "Didn't find all of the available shards");
            } else {
                nFound++;
                ProcessingLoc proc = tracker.findOwner(shard, NAME_ONE);

                Assert.assertTrue(proc.isOwnedBy(NAME_ONE));
                Assert.assertTrue(! markedShards.contains(shard), "Ran process was already marked!");
                Assert.assertTrue(toFind.contains(shard), "Claimed shard wasn't one of the unmarked!");
            }
        }
    }

    @Test(dataProvider = "simpleData", enabled = true)
    public void testMarkedProcesses(TestTarget test) {
        GenomeLocProcessingTracker tracker = test.getTracker();
        List<GenomeLoc> shards = test.getShards();
        logger.warn("testMarkedProcesses " + test);

        List<GenomeLoc> markedShards = new ArrayList<GenomeLoc>();

        for ( int i = 0; i < shards.size(); i++ ) {
            if ( i % 2 == 0 ) {
                markedShards.add(shards.get(i));
                tracker.claimOwnership(shards.get(i), NAME_TWO);
            }
        }

        for ( GenomeLoc shard : shards ) {
            ProcessingLoc proc = tracker.claimOwnership(shard,NAME_ONE);

            Assert.assertTrue(proc.isOwnedBy(NAME_ONE) || proc.isOwnedBy(NAME_TWO));

            if ( proc.isOwnedBy(NAME_ONE) )
                Assert.assertTrue(! markedShards.contains(shard), "Ran process was already marked!");
            else
                Assert.assertTrue(markedShards.contains(shard), "Unran process wasn't marked");

            if ( ! markedShards.contains(shard) ) {
                Assert.assertEquals(tracker.findOwner(shard, NAME_ONE), proc);
            }
        }
    }

    public class TestThread implements Callable<Integer> {
        public TestTarget test;
        public String name;
        public List<GenomeLoc> ran, toRun;
        boolean useIterator;

        public TestThread(TestTarget test, int count, List<GenomeLoc> toRun, boolean useIterator) {
            this.test = test;
            this.toRun = toRun;
            this.name = "thread" + count;
            this.ran = new ArrayList<GenomeLoc>();
            this.useIterator = useIterator;
        }

        public Integer call() {
            //logger.warn(String.format("Call() Thread %s", name));
            if ( useIterator ) {
                for ( GenomeLoc shard : test.getTracker().onlyOwned(toRun.iterator(), name) ) {
                    if ( shard != null ) { // ignore the unclaimable end of the stream
                        ran.add(shard);
                        // do some work here
                        for ( int sum =0, i = 0; i < 100000; i++) sum += i;
                    }
                }

            } else {
                for ( GenomeLoc shard : toRun ) {
                    //System.out.printf("Claiming ownership in %s on %s%n", name, shard);
                    ProcessingLoc proc = test.getTracker().claimOwnership(shard,name);
                    //System.out.printf("  => ownership of %s is %s (I own? %b)%n", shard, proc.getOwner(), proc.isOwnedBy(name));
                    if ( proc.isOwnedBy(name) ) {
                        ran.add(proc.getLocation());
                        // do some work here
                        for ( int sum =0, i = 0; i < 100000; i++) sum += i;
                    }
                    //logger.warn(String.format("Thread %s on %s -> owned by %s", name, shard, proc.getOwner()));
                }
            }

            return 1;
        }
    }

    private static TestThread findOwner(String name, List<TestThread> threads) {
        for ( TestThread thread : threads ) {
            if ( thread.name.equals(name) )
                return thread;
        }
        return null;
    }

    private static final <T> void assertAllThreadsFinished(List<Future<T>> futures) {
        try {
            for ( Future f : futures ) {
                Assert.assertTrue(f.isDone(), "Thread never finished running");
                Assert.assertTrue(f.get() != null, "Finished successfully");
            }
        } catch (InterruptedException e) {
            Assert.fail("Thread failed to run to completion", e);
        } catch (ExecutionException e) {
            Assert.fail("Thread generated an exception", e);
        }
    }

    private static final List<GenomeLoc> subList(List<GenomeLoc> l, int i) {
        List<GenomeLoc> r = new ArrayList<GenomeLoc>();
        for ( int j = 0; j < l.size(); j++ ) {
            if ( j % i == 0 )
                r.add(l.get(j));
        }

        return r;
    }

    @Test(dataProvider = "threadData", enabled = true)
    public void testThreadedProcessesLowLevelFunctions(TestTarget test) {
        testThreading(test, false);
    }

    @Test(dataProvider = "threadData", enabled = true)
    public void testThreadedProcessesIterator(TestTarget test) {
        testThreading(test, true);
    }

    private void testThreading(TestTarget test, boolean useIterator) {
        if ( ! test.isThreadSafe() )
            // skip tests that aren't thread safe
            return;

        // start up 3 threads
        logger.warn("ThreadedTesting " + test + " using iterator " + useIterator);
        List<TestThread> threads = new ArrayList<TestThread>();
        for ( int i = 0; i < 4; i++) {
            List<GenomeLoc> toRun = subList(test.getShards(), i+1);
            TestThread thread = new TestThread(test, i, toRun, useIterator);
            threads.add(thread);
        }
        ExecutorService exec = java.util.concurrent.Executors.newFixedThreadPool(threads.size());

        try {
            List<Future<Integer>> results = exec.invokeAll(threads, 300, TimeUnit.SECONDS);
            GenomeLocProcessingTracker tracker = test.getTracker();
            List<GenomeLoc> shards = test.getShards();

            for ( TestThread thread : threads )
                logger.warn(String.format("TestThread %s ran %d jobs of %d to run", thread.name, thread.ran.size(), thread.toRun.size()));

            assertAllThreadsFinished(results);

            // we ran everything
            Assert.assertEquals(tracker.updateAndGetProcessingLocs(NAME_ONE).size(), shards.size(), "Not all shards were run");

            for ( GenomeLoc shard : shards ) {
                Assert.assertTrue(tracker.locIsOwned(shard, NAME_ONE), "Unowned shard");

                ProcessingLoc proc = tracker.findOwner(shard, NAME_ONE);
                Assert.assertNotNull(proc, "Proc was null");

                Assert.assertNotNull(proc.getOwner(), "Owner was null");
                Assert.assertEquals(proc.getLocation(), shard, "Shard loc doesn't make ProcessingLoc");

                TestThread owner = findOwner(proc.getOwner(), threads);
                Assert.assertNotNull(owner, "Couldn't find owner");

                Assert.assertTrue(owner.ran.contains(shard), "Owner doesn't contain ran shard");

                for ( TestThread thread : threads )
                    if ( ! proc.isOwnedBy(thread.name) && thread.ran.contains(shard) )
                        Assert.fail("Shard appears in another run list: proc=" + proc + " shard=" + shard + " also in jobs of " + thread.name + " obj=" + thread.ran.get(thread.ran.indexOf(shard)));

            }
        }  catch (InterruptedException e) {
            Assert.fail("Thread failure", e);
        }
    }
}