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

import java.util.concurrent.TimeUnit;

/**
 * Defines a class of simple 'timing' facilities, which can be used
 * to measure durations.
 *
 * Note that this code is not thread-safe.  If you have a single timer
 * being started and stopped by multiple threads you will need to protect the
 * calls to avoid meaningless results of having multiple starts and stops
 * called sequentially.
 *
 * User: nc6@sanger.ac.uk
 * Date: Jan 13th, 2014
 * Time: 13:45
 */
public interface Timer {

    /**
     * @return the name associated with this timer
     */
    public String getName();

    /**
     * Starts the timer running, and sets the elapsedTimeNano time to 0.  This is equivalent to
     * resetting the time to have no history at all.
     *
     * @return this object, for programming convenience
     */
    public Timer start();

    /**
     * Starts the timer running, without resetting the elapsedTimeNano time.  This function may be
     * called without first calling start().  The only difference between start and restart
     * is that start resets the elapsedTimeNano time, while restart does not.
     *
     * @return this object, for programming convenience
     */
    public Timer restart();

    /**
     * @return is this timer running?
     */
    public boolean isRunning();

    /**
     * @return A convenience function to obtain the current time in milliseconds from this timer
     */
    public long currentTime();

    /**
     * Stops the timer.  Increases the elapsedTimeNano time by difference between start and now.
     *
     * It's ok to call stop on a timer that's not running.  It has no effect on the timer.
     *
     * @return this object, for programming convenience
     */
    public Timer stop();

    /**
     * Returns the total elapsedTimeNano time of all start/stops of this timer.  If the timer is currently
     * running, includes the difference from currentTime() and the start as well
     *
     * @return this time, in seconds
     */
    public double getElapsedTime();

    /**
     * @see #getElapsedTime() but returns the result in nanoseconds
     *
     * @return the elapsed time in nanoseconds
     */
    public long getElapsedTimeNano();

    /**
     * Add the elapsed time from toAdd to this elapsed time
     *
     * @param toAdd the timer whose elapsed time we want to add to this timer
     */
    public void addElapsed(final Timer toAdd);
}