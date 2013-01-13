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

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * Abstract base class to coordinating data processing by a collecting for processes / threads.
 *
 * Conceptually, the genome is viewed as a collection of non-overlapping genome location:
 *
 * chr1:1-10
 * chr1:11-20
 * chr1:21-30
 * etc.
 *
 * This class, and it's concrete derived classes, provide the ability to claim individual locations
 * as "mine", and exclude other processes / threads from processing them.  At the lowest-level this
 * is implemented by the claimOwnership(loc, name) function, that returns true if loc free (unclaimed)
 * and makes name the owner of loc.  High-level, and more efficient operations provide claiming
 * iterators over streams of objects implementing the HasGenomeLocation interface, so that you can
 * write code that looks like:
 *
 * for ( GenomeLoc ownedLoc : onlyOwned(allLocsToProcess.iterator) ) {
 *   doSomeWork(ownedLoc)
 *
 * Much of the code in this class is actually surrounding debugging and performance metrics code.
 * The actual synchronization code is separated out into the ClosableReentrantLock() system
 * and the two abstract functions:
 *
 *  protected abstract void registerNewLocs(Collection<ProcessingLoc> plocs);
 *  protected abstract Collection<ProcessingLoc> readNewLocs();
 *
 * That maintain the state of the tracker.
 *
 * That is, the ProcessingTracker is made of two components: a thread / process locking system and
 * a subclass that implements the methods to record new claimed state changes and to read out updates
 * that may have occurred by another thread or process.
 *
 * NOTE: this class assumes that all threads / processes are working with the same set of potential
 * GenomeLocs to own.  Claiming chr1:1-10 and then chr1:5-6 is allowed by the system.  Basically,
 * you only can stake claim to GenomeLocs that are .equal().
 */
public abstract class GenomeLocProcessingTracker {
    private final static Logger logger = Logger.getLogger(FileBackedGenomeLocProcessingTracker.class);
    private final static SimpleDateFormat STATUS_FORMAT = new SimpleDateFormat("HH:mm:ss,SSS");
    private final static int DEFAULT_OWNERSHIP_ITERATOR_SIZE = 1;

    /**
     * Useful state strings for printing status
     */
    private final static String GOING_FOR_LOCK = "going_for_lock";
    private final static String RELEASING_LOCK = "releasing_lock";
    private final static String HAVE_LOCK = "have_lock";
    private final static String RUNNING = "running";

    /**
     * A map, for efficiency, that allows quick lookup of the processing loc for a
     * given GenomeLoc.  The map points from loc -> loc / owner as a ProcessingLoc
     */
    private final Map<GenomeLoc, ProcessingLoc> processingLocs;

    /**
     * The locking object used to protect data from simulatanous access by multiple
     * threads or processes.
     */
    private final ClosableReentrantLock lock;

    /** A stream for writing status messages.  Can be null if we aren't writing status */
    private final PrintStream status;

    //
    // Timers for recording performance information
    // Note -- these cannot be used because this class isn't thread safe, and neither are the
    // timers, so they result in invalid operations w.r.t. the SimpleTimer contract
    //
//    protected final SimpleTimer writeTimer = new SimpleTimer("writeTimer");
//    protected final SimpleTimer readTimer = new SimpleTimer("readTimer");
//    protected final SimpleTimer lockWaitTimer = new SimpleTimer("lockWaitTimer");
    protected final SimpleTimer timer = new SimpleTimer();
    protected long nLocks = 0, nWrites = 0, nReads = 0;

    // --------------------------------------------------------------------------------
    //
    // Creating ProcessingTrackers
    //
    // --------------------------------------------------------------------------------
    public GenomeLocProcessingTracker(ClosableReentrantLock lock, PrintStream status) {
        this.processingLocs = new HashMap<GenomeLoc, ProcessingLoc>();
        this.status = status;
        this.lock = lock;
        printStatusHeader();
    }

    // --------------------------------------------------------------------------------
    //
    // Code to override to change the dynamics of the the GenomeLocProcessingTracker
    //
    // --------------------------------------------------------------------------------

    protected void close() {
        lock.close();
        if ( status != null ) status.close();
    }

    /**
     * Takes a collection of newly claimed (i.e., previous unclaimed) genome locs
     * and the name of their owner and "registers" this data in some persistent way that's
     * visible to all threads / processes communicating via this GenomeLocProcessingTracker.
     *
     * Could be a in-memory data structure (a list) if we are restricting ourselves to intra-memory
     * parallelism, a locked file on a shared file system, or a server we communicate with.
     *
     * @param plocs
     */
    protected abstract void registerNewLocs(Collection<ProcessingLoc> plocs);

    /**
     * The inverse of the registerNewLocs() function.  Looks at the persistent data store
     * shared by all threads / processes and returns the ones that have appeared since the last
     * call to readNewLocs().  Note that we expect the pair of registerNewLocs and readNewLocs to
     * include everything, even locs registered by this thread / process.  For example:
     *
     * readNewLocs() => List()
     * registerNewLocs(List(x, y,)) => void
     * readNewLocs() => List(x,y))
     *
     * even for this thread or process.
     * @return
     */
    protected abstract Collection<ProcessingLoc> readNewLocs();


    // --------------------------------------------------------------------------------
    //
    // Code to claim intervals for processing and query for their ownership
    //
    // --------------------------------------------------------------------------------

    /**
     * Queries the current database if a location is owned.  Does not guarantee that the
     * loc can be owned in a future call, though.
     *
     * @param loc
     * @return
     */
    public final boolean locIsOwned(GenomeLoc loc, String id) {
        return findOwner(loc, id) != null;
    }

    /**
     * The workhorse routine.  Attempt to claim processing ownership of loc, with my name.
     * This is an atomic operation -- other threads / processes will wait until this function
     * returns.  The return result is the ProcessingLoc object describing who owns this
     * location.  If the location isn't already claimed and we now own the location, the pl owner
     * will be myName.  Otherwise, the name of the owner can found in the pl.
     *
     * @param loc
     * @param myName
     * @return
     */
    public final ProcessingLoc claimOwnership(final GenomeLoc loc, final String myName) {
        // processingLocs is a shared memory synchronized object, and this
        // method is synchronized, so we can just do our processing
        return new WithLock<ProcessingLoc>(myName) {
            public ProcessingLoc doBody() {
                ProcessingLoc owner = findOwner(loc, myName);
                if ( owner == null ) { // we are unowned
                    owner = new ProcessingLoc(loc, myName);
                    registerNewLocsWithTimers(Arrays.asList(owner), myName);
                }
                return owner;
            }
        }.run();
    }


    // --------------------------------------------------------------------------------
    //
    // High-level iterator-style interface to claiming ownership
    //
    // --------------------------------------------------------------------------------

    /**
     * A higher-level, and more efficient, interface to obtain the next location we own.  Takes an
     * iterator producing objects that support the getLocation() interface, and returns the next
     * object in that stream that we can claim ownership of.  Returns null if we run out of elements
     * during the iteration.
     *
     * Can be more efficiently implemented in subclasses to avoid multiple unlocking
     *
     * @param iterator
     * @param myName
     * @return
     */
    public final <T extends HasGenomeLocation> T claimOwnershipOfNextAvailable(Iterator<T> iterator, String myName) {
        OwnershipIterator<T> myIt = new OwnershipIterator<T>(iterator, myName, 1);
        return myIt.next();
    }

    public final <T extends HasGenomeLocation> Iterable<T> onlyOwned(Iterator<T> iterator, String myName) {
        return new OwnershipIterator<T>(iterator, myName);
    }

    private final class OwnershipIterator<T extends HasGenomeLocation> implements Iterator<T>, Iterable<T> {
        private final Iterator<T> subit;
        private final String myName;
        private final Queue<T> cache;
        private final int cacheSize;

        public OwnershipIterator(Iterator<T> subit, String myName) {
            this(subit, myName, DEFAULT_OWNERSHIP_ITERATOR_SIZE);
        }

        public OwnershipIterator(Iterator<T> subit, String myName, int cacheSize) {
            this.subit = subit;
            this.myName = myName;
            cache = new LinkedList<T>();
            this.cacheSize = cacheSize;
        }

        /**
         * Will return true for all elements of subit, even if we can't get ownership of some of the future
         * elements and so will return null there
         * @return
         */
        public final boolean hasNext() {
            return cache.peek() != null || subit.hasNext();
        }

        /**
         * High performance iterator that only locks and unlocks once per claimed object found.  Avoids
         * locking / unlocking for each query
         *
         * @return an object of type T owned by this thread, or null if none of the remaining object could be claimed
         */
        public final T next() {
            if ( cache.peek() != null)
                return cache.poll();
            else {
                // cache is empty, we need to fill up the cache and return the first element of the queue
                return new WithLock<T>(myName) {
                    public T doBody() {
                        // read once the database of owners at the start
                        updateAndGetProcessingLocs(myName);

                        boolean done = false;
                        Queue<ProcessingLoc> pwns = new LinkedList<ProcessingLoc>(); // ;-)
                        while ( !done && cache.size() < cacheSize && subit.hasNext() ) {
                            final T elt = subit.next();
                            GenomeLoc loc = elt.getLocation();

                            ProcessingLoc owner = processingLocs.get(loc);

                            if ( owner == null ) { // we are unowned
                                owner = new ProcessingLoc(loc, myName);
                                pwns.offer(owner);
                                if ( ! cache.offer(elt) ) throw new ReviewedStingException("Cache offer unexpectedly failed");
                                if ( GenomeLoc.isUnmapped(loc) ) done = true;
                            }
                            // if not, we continue our search
                        }

                        registerNewLocsWithTimers(pwns, myName);

                        // we've either filled up the cache or run out of elements.  Either way we return
                        // the first element of the cache. If the cache is empty, we return null here.
                        return cache.poll();
                    }
                }.run();
            }
        }

        public final void remove() {
            throw new UnsupportedOperationException();
        }

        public final Iterator<T> iterator() {
            return this;
        }
    }

    // --------------------------------------------------------------------------------
    //
    // private / protected low-level accessors / manipulators and utility functions
    //
    // --------------------------------------------------------------------------------

    /**
     * Useful debugging function that returns the ProcessingLoc who owns loc.  ID
     * is provided for debugging purposes
     * @param loc
     * @param id
     * @return
     */
    protected final ProcessingLoc findOwner(GenomeLoc loc, String id) {
        // fast path to check if we already have the existing genome loc in memory for ownership claims
        // getProcessingLocs() may be expensive [reading from disk, for example] so we shouldn't call it
        // unless necessary
        ProcessingLoc x = processingLocs.get(loc);
        return x == null ? updateAndGetProcessingLocs(id).get(loc) : x;
    }

    /**
     * Returns the list of currently owned locations, updating the database as necessary.
     * DO NOT MODIFY THIS MAP! As with all parallelizing data structures, the list may be
     * out of date immediately after the call returns, or may be updating on the fly.
     * @return
     */
    protected final Map<GenomeLoc, ProcessingLoc> updateAndGetProcessingLocs(String myName) {
        return new WithLock<Map<GenomeLoc, ProcessingLoc>>(myName) {
            public Map<GenomeLoc, ProcessingLoc> doBody() {
//                readTimer.restart();
                for ( ProcessingLoc p : readNewLocs() )
                    processingLocs.put(p.getLocation(), p);
//                readTimer.stop();
                nReads++;
                return processingLocs;
            }
        }.run();
    }

    /**
     * Wrapper around registerNewLocs that also times the operation
     *
     * @param plocs
     * @param myName
     */
    protected final void registerNewLocsWithTimers(Collection<ProcessingLoc> plocs, String myName) {
//        writeTimer.restart();
        registerNewLocs(plocs);
        nWrites++;
//        writeTimer.stop();
    }

    private final void printStatusHeader() {
        if ( status != null ) status.printf("process.id\thr.time\ttime\tstate%n");
    }

    private final void printStatus(String id, long machineTime, String state) {
        // prints a line like processID human-readable-time machine-time state
        if ( status != null  ) {
            status.printf("%s\t%s\t%d\t%s%n", id, STATUS_FORMAT.format(machineTime), machineTime, state);
            status.flush();
        }
    }


    /**
     * Lock the data structure, preventing other threads / processes from reading and writing to the
     * common store
     * @param id the name of the process doing the locking
     */
    private final void lock(String id) {
        //lockWaitTimer.restart();
        boolean hadLock = lock.ownsLock();
        if ( ! hadLock ) {
            nLocks++;
            //printStatus(id, lockWaitTimer.currentTime(), GOING_FOR_LOCK);
        }
        lock.lock();
        //lockWaitTimer.stop();
        //if ( ! hadLock ) printStatus(id, lockWaitTimer.currentTime(), HAVE_LOCK);
    }

    /**
     * Unlock the data structure, allowing other threads / processes to read and write to the common store
     * @param id the name of the process doing the unlocking
     */
    private final void unlock(String id) {
        if ( lock.getHoldCount() == 1 ) printStatus(id, timer.currentTime(), RELEASING_LOCK);
        lock.unlock();
        if ( ! lock.ownsLock() ) printStatus(id, timer.currentTime(), RUNNING);
    }

    // useful code for getting
    public final long getNLocks() { return nLocks; }
    public final long getNReads() { return nReads; }
    public final long getNWrites() { return nWrites; }
//    public final double getTimePerLock() { return lockWaitTimer.getElapsedTime() / Math.max(nLocks, 1); }
//    public final double getTimePerRead() { return readTimer.getElapsedTime() / Math.max(nReads,1); }
//    public final double getTimePerWrite() { return writeTimer.getElapsedTime() / Math.max(nWrites,1); }

    // --------------------------------------------------------------------------------
    //
    // Java-style functional form for with lock do { x };
    //
    // --------------------------------------------------------------------------------

    /**
     * Private utility class that executes doBody() method with the lock() acquired and
     * handles property unlock()ing the system, even if an error occurs.  Allows one to write
     * clean code like:
     *
     * new WithLock<Integer>(name) {
     *   public Integer doBody() { doSomething(); return 1; }
     * }.run()
     *
     * @param <T> the return type of the doBody() method
     */
    private abstract class WithLock<T> {
        private final String myName;

        public WithLock(String myName) {
            this.myName = myName;
        }

        protected abstract T doBody();

        public T run() {
            boolean locked = false;
            try {
                lock(myName);
                locked = true;
                return doBody();
            } finally {
                if (locked) unlock(myName);
            }
        }
    }

    // --------------------------------------------------------------------------------
    //
    // main function for testing performance
    //
    // --------------------------------------------------------------------------------
    public static void main(String[] args) {
        //BasicConfigurator.configure();

        final String ref = args[0];
        final File file = new File(args[1]);
        final int cycles = Integer.valueOf(args[2]);

        File referenceFile = new File(ref);
        try {
            final IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(referenceFile);
            final String chr1 = fasta.getSequenceDictionary().getSequence(1).getSequenceName();
            final GenomeLocParser genomeLocParser = new GenomeLocParser(fasta);

            final class MyTest {
                String name;
                GenomeLocProcessingTracker tracker;

                MyTest(String name, GenomeLocProcessingTracker tracker) {
                    this.name = name;
                    this.tracker = tracker;
                }

                public void execute(int cycles) {
                    SimpleTimer delta = new SimpleTimer("delta");
                    SimpleTimer timer = new SimpleTimer("none");

                    if ( file.exists() ) file.delete();
                    timer.start();
                    delta.start();
                    for ( int i = 1; i < cycles; i++ ) {
                        tracker.claimOwnership(genomeLocParser.createGenomeLoc(chr1, i, i+1), "ABCDEFGHIJKL");
                        if ( i % 1000 == 0 ) {
                            System.out.printf("%s\t%d\t%d\t%.4f\t%.4f%n", name, i, timer.currentTime(), timer.getElapsedTime(), delta.getElapsedTime() );
                            delta.restart();
                        }
                    }
                }
            }

            System.out.printf("name\tcycle\tcurrent.time\telapsed.time\tdelta%n");
            new MyTest("in-memory", new SharedMemoryGenomeLocProcessingTracker(new ClosableReentrantLock())).execute(cycles);
            new MyTest("nio", new FileBackedGenomeLocProcessingTracker(file, genomeLocParser, new ClosableReentrantLock(), null)).execute(cycles);
            new MyTest("nio-file-lock", new FileBackedGenomeLocProcessingTracker(file, genomeLocParser, new SharedFileThreadSafeLock(file,1), null)).execute(cycles);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }
    }
}
