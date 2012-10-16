/*
 * Copyright (c) 2012, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.tools;

import net.sf.picard.util.PeekableIterator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.collections.ExperimentalNestedIntegerArray;
import org.broadinstitute.sting.utils.collections.LoggingNestedIntegerArray;
import org.broadinstitute.sting.utils.collections.LoggingNestedIntegerArray.NestedIntegerArrayOperation;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.recalibration.RecalDatum;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Semaphore;
import java.util.concurrent.TimeUnit;

/**
 * Utility to profile the scalability by # of threads of the NestedIntegerArray data structure (specifically
 * a NestedIntegerArray of RecalDatum as used by the BQSR).
 *
 * Allows toggling between the original NestedIntegerArray implementation and an ExperimentalNestedIntegerArray
 * implementation via the --useExperimentalArrays option.
 *
 * Requires a log file generated from an actual BQSR run with the --recal_table_update_log option giving
 * the sequence of GET and PUT operations to perform on the data structure.
 *
 * Given such a log file, dispatches the various GET/PUT operations to a thread pool of size specified
 * by the --threadPoolSize option, and outputs the elapsed wall clock time at the end.
 *
 * @author David Roazen
 */
public class ProfileNestedIntegerArray extends CommandLineProgram {

    private static Logger logger = Logger.getLogger(ProfileNestedIntegerArray.class);

    @Argument(fullName = "operationLog", shortName = "operationLog", doc = "File containing output from a LoggingNestedIntegerArray containing the update operations to perform during testing", required = true)
    private File operationLog = null;

    @Argument(fullName = "threadPoolSize", shortName = "threadPoolSize", doc = "Size of the thread pool to use for this test", required = false)
    private int threadPoolSize = 1;

    @Argument(fullName = "useExperimentalArrays", shortName = "useExperimentalArrays", doc = "Use the experimental NestedIntegerArray implementation?", required = false)
    private boolean useExperimentalArrays = false;

    @Argument(fullName = "maxEnqueuedTasks", shortName = "maxEnqueuedTasks", doc = "Maximum number of tasks that can be submitted to the thread pool at once", required = false)
    private int maxEnqueuedTasks = 10000;

    @Argument(fullName = "operationsPerThread", shortName = "operationsPerThread", doc = "Number of array operations to execute within each thread", required = false)
    private int operationsPerThread = 100000;

    @Argument(fullName = "operationBufferSize", shortName = "operationBufferSize", doc = "Number of operations to load at a time from the log file before dispatching them for execution", required = false)
    private int operationBufferSize = 20000000;

    @Argument(fullName = "debug", shortName = "debug", doc = "Output debugging information", required = false)
    private boolean debug = false;

    private PeekableIterator<String> operationLogIterator;
    private Iterator<ArrayOperation> operationIterator;
    private List<BulkArrayOperationRunner> bulkArrayOperationBuffer;
    private Map<String, NestedIntegerArray<RecalDatum>> arrays;
    private ExecutorService threadPool;
    private Semaphore threadPoolSlot;

    protected int execute() throws Exception {
        operationLogIterator = new PeekableIterator<String>(new XReadLines(operationLog, false));

        arrays = new HashMap<String, NestedIntegerArray<RecalDatum>>();
        initializeArrays();

        operationIterator = new ArrayOperationIterator();
        bulkArrayOperationBuffer = new ArrayList<BulkArrayOperationRunner>(operationBufferSize / operationsPerThread + 1);

        threadPool = Executors.newFixedThreadPool(threadPoolSize);
        threadPoolSlot = new Semaphore(maxEnqueuedTasks);

        logger.info("Running test with settings:");
        logger.info(String.format("operationLog=%s threadPoolSize=%d maxEnqueuedTasks=%d operationsPerThread=%d operationBufferSize=%d useExperimentalArrays=%b",
                                  operationLog, threadPoolSize, maxEnqueuedTasks, operationsPerThread, operationBufferSize, useExperimentalArrays));

        SimpleTimer wallClockTimer = new SimpleTimer("wallClock");
        wallClockTimer.start();
        runTest();
        wallClockTimer.stop();

        logger.info(String.format("Total wall clock time: %.2f seconds (%.2f minutes)", wallClockTimer.getElapsedTime(), wallClockTimer.getElapsedTime() / 60.0));

        return 0;
    }

    private void initializeArrays() {
        while ( operationLogIterator.hasNext() && operationLogIterator.peek().startsWith(LoggingNestedIntegerArray.HEADER_LINE_PREFIX) ) {
            String headerLine = operationLogIterator.next().substring(LoggingNestedIntegerArray.HEADER_LINE_PREFIX.length());

            String[] tokens = headerLine.split("\t", -1);
            if ( tokens.length < 2 ) {
                throw new UserException.MalformedFile(operationLog, "Found a header line with too few tokens (no dimensions specified)");
            }

            String arrayLabel = tokens[0];
            int[] dimensions = new int[tokens.length - 1];

            for ( int i = 1; i < tokens.length; i++ ) {
                try {
                    dimensions[i - 1] = Integer.parseInt(tokens[i]);
                }
                catch ( NumberFormatException e ) {
                    throw new UserException.MalformedFile(operationLog, "Error parsing a numerical header field", e);
                }
            }

            arrays.put(arrayLabel, useExperimentalArrays ? new ExperimentalNestedIntegerArray<RecalDatum>(dimensions) :
                                                           new NestedIntegerArray<RecalDatum>(dimensions));

            if ( debug ) {
                logger.info(String.format("Created NestedIntegerArray %s with dimensions %s", arrayLabel, Arrays.toString(dimensions)));
            }
        }

        if ( arrays.size() == 0 ) {
            throw new UserException.MalformedFile(operationLog, "No array metadata was found in the header");
        }
    }

    private void runTest() {
        try {
            int maxObservedEnqueuedOperations = 0;
            long enqueuedOperationsRunningTotal = 0;
            int numObservations = 0;

            do {
                bufferUpcomingOperations();

                for  ( BulkArrayOperationRunner bulkOperation : bulkArrayOperationBuffer ) {
                    // Record stats so that we can make sure the thread pool is always using all its threads
                    // and that we're not I/O bound
                    int currentlyEnqueuedOperations = maxEnqueuedTasks - threadPoolSlot.availablePermits();
                    maxObservedEnqueuedOperations = Math.max(maxObservedEnqueuedOperations, currentlyEnqueuedOperations);
                    enqueuedOperationsRunningTotal += currentlyEnqueuedOperations;
                    numObservations++;

                    threadPoolSlot.acquire();
                    threadPool.execute(bulkOperation);
                }
            } while ( bulkArrayOperationBuffer.size() > 0 );

            threadPool.shutdown();
            if ( ! threadPool.awaitTermination(300, TimeUnit.SECONDS) ) {
                throw new ReviewedStingException("Final tasks in thread pool did not complete within a reasonable amount of time");
            }

            logger.info(String.format("Max observed enqueued operations: %d\tAverage # of enqueued operations: %.2f",
                                      maxObservedEnqueuedOperations, (double)enqueuedOperationsRunningTotal / numObservations));
        }
        catch ( InterruptedException e ) {
            threadPool.shutdownNow();
            throw new ReviewedStingException("Thread interrupted during execution");
        }
    }

    private void bufferUpcomingOperations() {
        bulkArrayOperationBuffer.clear();
        int totalOperationsLoaded = 0;

        while ( operationIterator.hasNext() && totalOperationsLoaded < operationBufferSize ) {
            List<ArrayOperation> operations = new ArrayList<ArrayOperation>(operationsPerThread);
            while ( operationIterator.hasNext() && operations.size() < operationsPerThread && totalOperationsLoaded < operationBufferSize ) {
                operations.add(operationIterator.next());
                totalOperationsLoaded++;
            }

            bulkArrayOperationBuffer.add(new BulkArrayOperationRunner(operations));
        }
    }

    private class BulkArrayOperationRunner implements Runnable {
        private List<ArrayOperation> operations;

        public BulkArrayOperationRunner( List<ArrayOperation> operations ) {
            this.operations = operations;
        }

        public void run() {
            try {
                for ( ArrayOperation operation : operations ) {
                    operation.run();
                }
            }
            finally {
                threadPoolSlot.release();
            }
        }
    }

    private class ArrayOperation implements Runnable {
        private NestedIntegerArrayOperation operationType;
        private String arrayLabel;
        private RecalDatum datum;
        private int[] keys;

        public ArrayOperation( NestedIntegerArrayOperation operationType, String arrayLabel, RecalDatum datum, int[] keys ) {
            this.operationType = operationType;
            this.arrayLabel = arrayLabel;
            this.datum = datum;
            this.keys = keys.clone();
        }

        public void run() {
            if ( debug ) {
                logger.info(String.format("Running task %s", toString()));
            }

            NestedIntegerArray<RecalDatum> array = arrays.get(arrayLabel);
            if ( array == null ) {
                throw new ReviewedStingException(String.format("Attempted to access non-existent array %s", arrayLabel));
            }

            switch ( operationType ) {
                case GET:
                    RecalDatum existingDatum = array.get(keys);
                    if ( existingDatum != null ) {
                        existingDatum.increment(true);  // Assume that isError is always true for profiling purposes --
                                                        // will not affect performance results, as all that matters is
                                                        // that we call a synchronized method on this RecalDatum
                    }
                    // If existingDatum was null, we'll assume that there is an upcoming PUT operation
                    break;
                case PUT:
                    array.put(datum, keys);
                    break;
            }
        }

        public String toString() {
            return String.format("[%s in table %s, datum: %s keys: %s]", operationType,
                                                                         arrayLabel,
                                                                         datum != null ? datum : "(none)",
                                                                         Arrays.toString(keys));
        }
    }

    private class ArrayOperationIterator implements Iterator<ArrayOperation>, Iterable<ArrayOperation> {

        private ArrayOperation nextOperation = null;

        public ArrayOperationIterator() {
            advance();
        }

        public boolean hasNext() {
            return nextOperation != null;
        }

        public ArrayOperation next() {
            if ( nextOperation == null ) {
                throw new NoSuchElementException("next() called when there are no more items");
            }

            ArrayOperation toReturn = nextOperation;
            advance();
            return toReturn;
        }

        private void advance() {
            if ( ! operationLogIterator.hasNext() ) {
                nextOperation = null;
                return;
            }

            String nextOperationLogEntry = operationLogIterator.next();
            String[] tokens = nextOperationLogEntry.split("\t", -1);

            if ( tokens.length < 4 ) {
                throw new UserException.MalformedFile(operationLog, String.format("Found an array operation log entry with less than 4 fields: %s", nextOperationLogEntry));
            }

            String arrayLabel = tokens[0];

            NestedIntegerArrayOperation operationType;
            try {
                operationType = NestedIntegerArrayOperation.valueOf(tokens[1]);
            }
            catch ( IllegalArgumentException e ) {
                throw new UserException.MalformedFile(operationLog, String.format("Illegal operation type %s in log entry %s", tokens[1], nextOperationLogEntry));
            }

            RecalDatum datum = tokens[2].length() > 0 ? parseRecalDatumString(tokens[2]) : null;

            int[] keys = new int[tokens.length - 3];
            for ( int i = 3; i < tokens.length; i++ ) {
                try {
                    keys[i - 3] = Integer.parseInt(tokens[i]);
                }
                catch ( NumberFormatException e ) {
                    throw new UserException.MalformedFile(operationLog, String.format("Found an array operation log entry with non-integer key: %s", nextOperationLogEntry), e);
                }
            }

            nextOperation = new ArrayOperation(operationType, arrayLabel, datum, keys);
        }

        private RecalDatum parseRecalDatumString( String recalDatumString ) {
            String[] tokens = recalDatumString.split(",", -1);
            double numObservations, numMismatches;
            byte quality;

            try {
                numObservations = Double.parseDouble(tokens[0]);
                numMismatches = Double.parseDouble(tokens[1]);
                quality = (byte)Math.round(Double.parseDouble(tokens[2]));   // using empirical quality as substitute for reported quality for profiling purposes
            }
            catch ( NumberFormatException e ) {
                throw new UserException.MalformedFile(operationLog, String.format("Failed to parse RecalDatum string %s", recalDatumString), e);
            }

            return new RecalDatum(numObservations, numMismatches, quality);
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }

        public Iterator<ArrayOperation> iterator() {
            return this;
        }
    }

    public static void main ( String[] args ) {
        try {
            ProfileNestedIntegerArray instance = new ProfileNestedIntegerArray();
            start(instance, args);
            System.exit(CommandLineProgram.result);
        } catch ( UserException e ) {
            exitSystemWithUserError(e);
        } catch ( Exception e ) {
            exitSystemWithError(e);
        }
    }
}
