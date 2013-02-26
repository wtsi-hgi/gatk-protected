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

package org.broadinstitute.sting.gatk.executive;

import com.google.java.contract.Ensures;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.datasources.reads.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.iterators.NullSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.resourcemanagement.ThreadAllocation;
import org.broadinstitute.sting.gatk.traversals.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.AutoFormattingTime;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.activeregion.ExperimentalActiveRegionShardType;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.progressmeter.ProgressMeter;
import org.broadinstitute.sting.utils.threading.ThreadEfficiencyMonitor;

import javax.management.JMException;
import javax.management.MBeanServer;
import javax.management.ObjectName;
import java.io.File;
import java.lang.management.ManagementFactory;
import java.util.*;
import java.util.concurrent.TimeUnit;


/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 26, 2009
 * Time: 12:37:23 PM
 *
 * General base class for all scheduling algorithms
 * Shards and schedules data in manageable chunks.
 *
 * Creates N TraversalEngines for each data thread for the MicroScheduler.  This is necessary
 * because in the HMS case you have multiple threads executing a traversal engine independently, and
 * these engines may need to create separate resources for efficiency or implementation reasons.  For example,
 * the nanoScheduler creates threads to implement the traversal, and this creation is instance specific.
 * So each HMS thread needs to have it's own distinct copy of the traversal engine if it wants to have
 * N data threads x M nano threads => N * M threads total.  These are borrowed from this microscheduler
 * and returned when done.  Also allows us to tracks all created traversal engines so this microscheduler
 * can properly shut them all down when the scheduling is done.
 *
 */
public abstract class MicroScheduler implements MicroSchedulerMBean {
    protected static final Logger logger = Logger.getLogger(MicroScheduler.class);

    /**
     * The list of all Traversal engines we've created in this micro scheduler
     */
    final List<TraversalEngine> allCreatedTraversalEngines = new LinkedList<TraversalEngine>();

    /**
     * All available engines.  Engines are borrowed and returned when a subclass is actually
     * going to execute the engine on some data.  This allows us to have N copies for
     * N data parallel executions, but without the dangerous code of having local
     * ThreadLocal variables.
     */
    final LinkedList<TraversalEngine> availableTraversalEngines = new LinkedList<TraversalEngine>();

    /**
     * Engines that have been allocated to a key already.
     */
    final HashMap<Object, TraversalEngine> allocatedTraversalEngines = new HashMap<Object, TraversalEngine>();

    /**
     * Counts the number of instances of the class that are currently alive.
     */
    private static int instanceNumber = 0;

    /**
     * The engine invoking this scheduler.
     */
    protected final GenomeAnalysisEngine engine;

    protected final IndexedFastaSequenceFile reference;

    private final SAMDataSource reads;
    protected final Collection<ReferenceOrderedDataSource> rods;

    private final MBeanServer mBeanServer;
    private final ObjectName mBeanName;

    /**
     * Threading efficiency monitor for tracking the resource utilization of the GATK
     *
     * may be null
     */
    ThreadEfficiencyMonitor threadEfficiencyMonitor = null;

    final ProgressMeter progressMeter;

    /**
     * MicroScheduler factory function.  Create a microscheduler appropriate for reducing the
     * selected walker.
     *
     * @param walker        Which walker to use.
     * @param reads         the informations associated with the reads
     * @param reference     the reference file
     * @param rods          the rods to include in the traversal
     * @param threadAllocation Number of threads to utilize.
     *
     * @return The best-fit microscheduler.
     */
    public static MicroScheduler create(GenomeAnalysisEngine engine, Walker walker, SAMDataSource reads, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods, ThreadAllocation threadAllocation) {
        if ( threadAllocation.isRunningInParallelMode() ) {
            logger.info(String.format("Running the GATK in parallel mode with %d total threads, " +
                    "%d CPU thread(s) for each of %d data thread(s), of %d processors available on this machine",
                    threadAllocation.getTotalNumThreads(),
                    threadAllocation.getNumCPUThreadsPerDataThread(),
                    threadAllocation.getNumDataThreads(),
                    Runtime.getRuntime().availableProcessors()));
            if ( threadAllocation.getTotalNumThreads() > Runtime.getRuntime().availableProcessors() )
                logger.warn(String.format("Number of requested GATK threads %d is more than the number of " +
                        "available processors on this machine %d", threadAllocation.getTotalNumThreads(),
                        Runtime.getRuntime().availableProcessors()));
//            if ( threadAllocation.getNumDataThreads() > 1 && threadAllocation.getNumCPUThreadsPerDataThread() > 1)
//                throw new UserException("The GATK currently doesn't support running with both -nt > 1 and -nct > 1");
        }

        if ( threadAllocation.getNumDataThreads() > 1 ) {
            if (walker.isReduceByInterval())
                throw new UserException.BadArgumentValue("nt", String.format("The analysis %s aggregates results by interval.  Due to a current limitation of the GATK, analyses of this type do not currently support parallel execution.  Please run your analysis without the -nt option.", engine.getWalkerName(walker.getClass())));

            if ( ! (walker instanceof TreeReducible) ) {
                throw badNT("nt", engine, walker);
            }
        }

        if ( threadAllocation.getNumCPUThreadsPerDataThread() > 1 && ! (walker instanceof NanoSchedulable) ) {
            throw badNT("nct", engine, walker);
        }

        if ( threadAllocation.getNumDataThreads() > 1 ) {
            return new HierarchicalMicroScheduler(engine, walker, reads, reference, rods, threadAllocation);
        } else {
            return new LinearMicroScheduler(engine, walker, reads, reference, rods, threadAllocation);
        }
    }

    private static UserException badNT(final String parallelArg, final GenomeAnalysisEngine engine, final Walker walker) {
        throw new UserException.BadArgumentValue(parallelArg,
                String.format("The analysis %s currently does not support parallel execution with %s.  " +
                        "Please run your analysis without the %s option.", engine.getWalkerName(walker.getClass()), parallelArg, parallelArg));
    }

    /**
     * Create a microscheduler given the reads and reference.
     *
     * @param walker  the walker to execute with
     * @param reads   The reads.
     * @param reference The reference.
     * @param rods    the rods to include in the traversal
     * @param threadAllocation the allocation of threads to use in the underlying traversal
     */
    protected MicroScheduler(final GenomeAnalysisEngine engine,
                             final Walker walker,
                             final SAMDataSource reads,
                             final IndexedFastaSequenceFile reference,
                             final Collection<ReferenceOrderedDataSource> rods,
                             final ThreadAllocation threadAllocation) {
        this.engine = engine;
        this.reads = reads;
        this.reference = reference;
        this.rods = rods;

        final File progressLogFile = engine.getArguments() == null ? null : engine.getArguments().performanceLog;

        // Creates uninitialized TraversalEngines appropriate for walker and threadAllocation,
        // and adds it to the list of created engines for later shutdown.
        for ( int i = 0; i < threadAllocation.getNumDataThreads(); i++ ) {
            final TraversalEngine traversalEngine = createTraversalEngine(walker, threadAllocation);
            allCreatedTraversalEngines.add(traversalEngine);
            availableTraversalEngines.add(traversalEngine);
        }

        // Create our progress meter
        this.progressMeter = new ProgressMeter(progressLogFile,
                availableTraversalEngines.peek().getTraversalUnits(),
                engine.getRegionsOfGenomeBeingProcessed());

        // Now that we have a progress meter, go through and initialize the traversal engines
        for ( final TraversalEngine traversalEngine : allCreatedTraversalEngines )
            traversalEngine.initialize(engine, progressMeter);

        // JMX does not allow multiple instances with the same ObjectName to be registered with the same platform MXBean.
        // To get around this limitation and since we have no job identifier at this point, register a simple counter that
        // will count the number of instances of this object that have been created in this JVM.
        int thisInstance = instanceNumber++;
        mBeanServer = ManagementFactory.getPlatformMBeanServer();
        try {
            mBeanName = new ObjectName("org.broadinstitute.sting.gatk.executive:type=MicroScheduler,instanceNumber="+thisInstance);
            mBeanServer.registerMBean(this, mBeanName);
        }
        catch (JMException ex) {
            throw new ReviewedStingException("Unable to register microscheduler with JMX", ex);
        }
    }

    /**
     * Really make us a traversal engine of the appropriate type for walker and thread allocation
     *
     * @return a non-null uninitialized traversal engine
     */
    @Ensures("result != null")
    private TraversalEngine createTraversalEngine(final Walker walker, final ThreadAllocation threadAllocation) {
        if (walker instanceof ReadWalker) {
            return new TraverseReadsNano(threadAllocation.getNumCPUThreadsPerDataThread());
        } else if (walker instanceof LocusWalker) {
            return new TraverseLociNano(threadAllocation.getNumCPUThreadsPerDataThread());
        } else if (walker instanceof DuplicateWalker) {
            return new TraverseDuplicates();
        } else if (walker instanceof ReadPairWalker) {
            return new TraverseReadPairs();
        } else if (walker instanceof ActiveRegionWalker) {
            switch (engine.getArguments().activeRegionShardType) {
                case LOCUSSHARD: return new TraverseActiveRegions();
                case READSHARD: return new ExperimentalReadShardTraverseActiveRegions();
                case ACTIVEREGIONSHARD: return new ExperimentalActiveRegionShardTraverseActiveRegions();
                default: throw new UnsupportedOperationException("Unable to determine traversal type, the walker is an unknown type of ActiveRegionWalker.");
            }
        } else {
            throw new UnsupportedOperationException("Unable to determine traversal type, the walker is an unknown type.");
        }
    }


    /**
     * Return the ThreadEfficiencyMonitor we are using to track our resource utilization, if there is one
     *
     * @return the monitor, or null if none is active
     */
    public ThreadEfficiencyMonitor getThreadEfficiencyMonitor() {
        return threadEfficiencyMonitor;
    }

    /**
     * Inform this Microscheduler to use the efficiency monitor used to create threads in subclasses
     *
     * @param threadEfficiencyMonitor
     */
    public void setThreadEfficiencyMonitor(final ThreadEfficiencyMonitor threadEfficiencyMonitor) {
        this.threadEfficiencyMonitor = threadEfficiencyMonitor;
    }

    /**
     * Should we stop all execution work and exit gracefully?
     *
     * Returns true in the case where some external signal or time limit has been received, indicating
     * that this GATK shouldn't continue executing.  This isn't a kill signal, it is really a "shutdown
     * gracefully at the next opportunity" signal.  Concrete implementations of the MicroScheduler
     * examine this value as often as reasonable and, if it returns true, stop what they are doing
     * at the next available opportunity, shutdown their resources, call notify done, and return.
     *
     * @return true if we should abort execution, or false otherwise
     */
    protected boolean abortExecution() {
        final boolean abort = engine.exceedsRuntimeLimit(progressMeter.getRuntimeInNanoseconds(), TimeUnit.NANOSECONDS);
        if ( abort ) {
            final AutoFormattingTime aft = new AutoFormattingTime(engine.getRuntimeLimitInNanoseconds(), -1, 4);
            logger.info("Aborting execution (cleanly) because the runtime has exceeded the requested maximum " + aft);
        }
        return abort;
    }

    /**
     * Walks a walker over the given list of intervals.
     *
     * @param walker        Computation to perform over dataset.
     * @param shardStrategy A strategy for sharding the data.
     *
     * @return the return type of the walker
     */
    public abstract Object execute(Walker walker, Iterable<Shard> shardStrategy);

    /**
     * Retrieves the object responsible for tracking and managing output.
     * @return An output tracker, for loading data in and extracting results.  Will not be null.
     */
    public abstract OutputTracker getOutputTracker();

    /**
     * Gets the an iterator over the given reads, which will iterate over the reads in the given shard.
     * @param shard the shard to use when querying reads.
     * @return an iterator over the reads specified in the shard.
     */
    protected StingSAMIterator getReadIterator(Shard shard) {
        return (!reads.isEmpty()) ? reads.seek(shard) : new NullSAMIterator();
    }

    /**
     * Must be called by subclasses when execute is done
     */
    protected void executionIsDone() {
        progressMeter.notifyDone(engine.getCumulativeMetrics().getNumIterations());
        printReadFilteringStats();
        shutdownTraversalEngines();

        // Print out the threading efficiency of this HMS, if state monitoring is enabled
        if ( threadEfficiencyMonitor != null ) {
            // include the master thread information
            threadEfficiencyMonitor.threadIsDone(Thread.currentThread());
            threadEfficiencyMonitor.printUsageInformation(logger);
        }
    }

    /**
     * Shutdown all of the created engines, and clear the list of created engines, dropping
     * pointers to the traversal engines
     */
    public synchronized void shutdownTraversalEngines() {
        // no longer applicable because engines are allocated to keys now
//        if ( availableTraversalEngines.size() != allCreatedTraversalEngines.size() )
//            throw new IllegalStateException("Shutting down TraversalEngineCreator but not all engines " +
//                    "have been returned.  Expected " + allCreatedTraversalEngines.size() + " but only " + availableTraversalEngines.size()
//                    + " have been returned");

        for ( final TraversalEngine te : allCreatedTraversalEngines)
            te.shutdown();

        allCreatedTraversalEngines.clear();
        availableTraversalEngines.clear();
    }

    /**
     * Prints out information about number of reads observed and filtering, if any reads were used in the traversal
     *
     * Looks like:
     *
     * INFO  10:40:47,370 MicroScheduler - 22 reads were filtered out during traversal out of 101 total (21.78%)
     * INFO  10:40:47,370 MicroScheduler -   -> 1 reads (0.99% of total) failing BadMateFilter
     * INFO  10:40:47,370 MicroScheduler -   -> 20 reads (19.80% of total) failing DuplicateReadFilter
     * INFO  10:40:47,370 MicroScheduler -   -> 1 reads (0.99% of total) failing FailsVendorQualityCheckFilter
     */
    private void printReadFilteringStats() {
        final ReadMetrics cumulativeMetrics = engine.getCumulativeMetrics();
        if ( cumulativeMetrics.getNumReadsSeen() > 0 ) {
            // count up the number of skipped reads by summing over all filters
            long nSkippedReads = 0L;
            for ( final long countsByFilter : cumulativeMetrics.getCountsByFilter().values())
                nSkippedReads += countsByFilter;

            logger.info(String.format("%d reads were filtered out during traversal out of %d total (%.2f%%)",
                    nSkippedReads,
                    cumulativeMetrics.getNumReadsSeen(),
                    100.0 * MathUtils.ratio(nSkippedReads, cumulativeMetrics.getNumReadsSeen())));

            for ( final Map.Entry<String, Long> filterCounts : cumulativeMetrics.getCountsByFilter().entrySet() ) {
                long count = filterCounts.getValue();
                logger.info(String.format("  -> %d reads (%.2f%% of total) failing %s",
                        count, 100.0 * MathUtils.ratio(count,cumulativeMetrics.getNumReadsSeen()), filterCounts.getKey()));
            }
        }
    }

    /**
     * Gets the engine that created this microscheduler.
     * @return The engine owning this microscheduler.
     */
    public GenomeAnalysisEngine getEngine() { return engine; }

    /**
     * Returns data source maintained by this scheduler
     * @return
     */
    public SAMDataSource getSAMDataSource() { return reads; }

    /**
     * Returns the reference maintained by this scheduler.
     * @return The reference maintained by this scheduler.
     */
    public IndexedFastaSequenceFile getReference() { return reference; }

    protected void cleanup() {
        try {
            mBeanServer.unregisterMBean(mBeanName);
        }
        catch (JMException ex) {
            throw new ReviewedStingException("Unable to unregister microscheduler with JMX", ex);
        }
    }

    /**
     * Returns a traversal engine suitable for use, associated with key
     *
     * Key is an arbitrary object that is used to retrieve the same traversal
     * engine over and over.  This can be important in the case where the
     * traversal engine has data associated with it in some other context,
     * and we need to ensure that the context always sees the same traversal
     * engine.  This happens in the HierarchicalMicroScheduler, where you want
     * the a thread executing traversals to retrieve the same engine each time,
     * as outputs are tracked w.r.t. that engine.
     *
     * If no engine is associated with key yet, pops the next available engine
     * from the available ones maintained by this
     * microscheduler.  Note that it's a runtime error to pop a traversal engine
     * from this scheduler if there are none available.  Callers that
     * once pop'd an engine for use must return it with returnTraversalEngine
     *
     * @param key the key to associate with this engine
     * @return a non-null TraversalEngine suitable for execution in this scheduler
     */
    @Ensures("result != null")
    protected synchronized TraversalEngine borrowTraversalEngine(final Object key) {
        if ( key == null ) throw new IllegalArgumentException("key cannot be null");

        final TraversalEngine engine = allocatedTraversalEngines.get(key);
        if ( engine == null ) {
            if ( availableTraversalEngines.isEmpty() )
                throw new IllegalStateException("no traversal engines were available");
            allocatedTraversalEngines.put(key, availableTraversalEngines.pop());
            return allocatedTraversalEngines.get(key);
        } else {
            return engine;
        }
    }

    /**
     * Return a borrowed traversal engine to this MicroScheduler, for later use
     * in another traversal execution
     *
     * @param key the key used to id the engine, provided to the borrowTraversalEngine function
     * @param traversalEngine the borrowed traversal engine.  Must have been previously borrowed.
     */
    protected synchronized void returnTraversalEngine(final Object key, final TraversalEngine traversalEngine) {
        if ( traversalEngine == null )
            throw new IllegalArgumentException("Attempting to push a null traversal engine");
        if ( ! allCreatedTraversalEngines.contains(traversalEngine) )
            throw new IllegalArgumentException("Attempting to push a traversal engine not created by this MicroScheduler" + engine);
        if ( ! allocatedTraversalEngines.containsKey(key) )
            throw new IllegalArgumentException("No traversal engine was never checked out with key " + key);

        // note there's nothing to actually do here, but a function implementation
        // might want to do something
    }
}
