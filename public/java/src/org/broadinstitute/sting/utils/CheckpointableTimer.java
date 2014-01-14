/*
* Copyright (c) 2012 The Broad Institute
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.sting.utils;


import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

import java.util.concurrent.TimeUnit;

import static java.lang.Math.abs;

/**
 * A timer using actual system time and thus supporting checkpoint/restart
 * semantics without breaking due to differing JVM clocks.
 *
 * Note that this code is not thread-safe.  If you have a single timer
 * being started and stopped by multiple threads you will need to protect the
 * calls to avoid meaningless results of having multiple starts and stops
 * called sequentially.
 *
 * User: depristo
 * Date: Dec 10, 2010
 * Time: 9:07:44 AM
 */
public class CheckpointableTimer implements Timer {

    protected static final double NANO_TO_SECOND_DOUBLE = 1.0 / TimeUnit.SECONDS.toNanos(1);

    /**
     * Allowable clock drift in nanoseconds.
     */
    private static final double CLOCK_DRIFT = 5000000;

    private final String name;

    /**
     * The difference between system time and nano time at construction.
     * This is used to detect checkpoint/restart events, and should be 
     * reset when a checkpoint/restart is detected.
     */
    private long nanoTimeOffset;

    /**
     * The elapsedTimeNano time in nanoSeconds of this timer.  The elapsedTimeNano time is the
     * sum of times between starts/restrats and stops.
     */
    private long elapsedTimeNano = 0l;

    /**
     * The start time of the last start/restart in nanoSeconds
     */
    private long startTimeNano = 0l;

    /**
     * Is this timer currently running (i.e., the last call was start/restart)
     */
    private boolean running = false;

    /**
     * Creates an anonymous simple timer
     */
    public CheckpointableTimer() {
        this("Anonymous");
    }

    /**
     * Creates a timer named name
     * @param name of the timer, must not be null
     */
    public CheckpointableTimer(final String name) {
        if ( name == null ) throw new IllegalArgumentException("CheckpointableTimer name cannot be null");
        this.name = name;

        this.nanoTimeOffset = getNanoOffset();
    }

    /**
     * @return the name associated with this timer
     */
    public synchronized String getName() {
        return name;
    }

    /**
     * Starts the timer running, and sets the elapsedTimeNano time to 0.  This is equivalent to
     * resetting the time to have no history at all.
     *
     * @return this object, for programming convenience
     */
    @Ensures("elapsedTimeNano == 0l")
    public synchronized CheckpointableTimer start() {
        elapsedTimeNano = 0l;
        return restart();
    }

    /**
     * Starts the timer running, without resetting the elapsedTimeNano time.  This function may be
     * called without first calling start().  The only difference between start and restart
     * is that start resets the elapsedTimeNano time, while restart does not.
     *
     * @return this object, for programming convenience
     */
    public synchronized CheckpointableTimer restart() {
        running = true;
        startTimeNano = currentTimeNano();
        return this;
    }

    /**
     * @return is this timer running?
     */
    public synchronized boolean isRunning() {
        return running;
    }

    /**
     * @return A convenience function to obtain the current time in milliseconds from this timer
     */
    public long currentTime() {
        return System.currentTimeMillis();
    }

    /**
     * @return A convenience function to obtain the current time in nanoSeconds from this timer
     */
    public long currentTimeNano() {
        return System.nanoTime();
    }

    /**
     * Stops the timer.  Increases the elapsedTimeNano time by difference between start and now.
     * This also checks that there hasn't been a massive change between the system and nano time
     * since the initiation - this would be indicative of a C/R. If there has been a change, the
     * elapsed time is discarded.
     *
     * It's ok to call stop on a timer that's not running.  It has no effect on the timer.
     *
     * @return this object, for programming convenience
     */
    @Requires("startTimeNano != 0l")
    public synchronized CheckpointableTimer stop() {
        if ( running ) {
            running = false;
            long currentOffset = getNanoOffset();
            if (abs(currentOffset - nanoTimeOffset) <= CLOCK_DRIFT) {
                elapsedTimeNano += currentTimeNano() - startTimeNano;
            }
            // Reset the drift meter to stay in sync.
            this.nanoTimeOffset = currentOffset;
        }
        return this;
    }

    /**
     * Returns the total elapsedTimeNano time of all start/stops of this timer.  If the timer is currently
     * running, includes the difference from currentTime() and the start as well
     *
     * @return this time, in seconds
     */
    public synchronized double getElapsedTime() {
        return nanoToSecondsAsDouble(getElapsedTimeNano());
    }

    protected static double nanoToSecondsAsDouble(final long nano) {
        return nano * NANO_TO_SECOND_DOUBLE;
    }

    /**
     * @see #getElapsedTime() but returns the result in nanoseconds
     *
     * @return the elapsed time in nanoseconds
     */
    public synchronized long getElapsedTimeNano() {
        return running ? (currentTimeNano() - startTimeNano + elapsedTimeNano) : elapsedTimeNano;
    }

    /**
     * Add the elapsed time from toAdd to this elapsed time
     *
     * @param toAdd the timer whose elapsed time we want to add to this timer
     */
    public synchronized void addElapsed(final Timer toAdd) {
        elapsedTimeNano += toAdd.getElapsedTimeNano();
    }

    /**
     * Get the current offset of nano time from system time.
     */
    private static long getNanoOffset() {
        return System.nanoTime() - (System.currentTimeMillis() * 1000);
    }
}
