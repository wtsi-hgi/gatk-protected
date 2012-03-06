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

package org.broadinstitute.sting.utils.recalibration;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * A general algorithm for quantizing quality score distributions to use a specific number of levels
 *
 * @author Mark Depristo
 * @since 3/2/12
 */
public class QualQuantizer {
    private static Logger logger = Logger.getLogger(QualQuantizer.class);

    final int nLevels, originalSize;
    final List<Long> nObservationsPerQual;

    /** Map from original qual (e.g., Q30) to new quantized qual (e.g., Q28) */
    final List<Byte> originalToQuantizedMap;
    final List<Double> penaltyPerQual;
    final TreeSet<QualInterval> quantizedIntervals;
    double overallPenalty;

    public QualQuantizer(final List<Long> nObservationsPerQual, final int nLevels) {
        this.nObservationsPerQual = nObservationsPerQual;
        this.nLevels = nLevels;

        this.originalSize = nObservationsPerQual.size();
        this.penaltyPerQual = new ArrayList<Double>(originalSize);

        // for QC.  These values should never appear after a correct run
        Collections.fill(penaltyPerQual, -1.0);

        this.quantizedIntervals = quantize();

        this.originalToQuantizedMap = intervalsToMap(quantizedIntervals);
    }

    public static class QualInterval implements Comparable<QualInterval> {
        final int qStart, qEnd;
        final long nObservations;
        final long nErrors;
        final int fixedQual;
        Set<QualInterval> subIntervals = new HashSet<QualInterval>();

        public QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors) {
            this(qStart, qEnd, nObservations, nErrors, -1);
        }

        public QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors, final int fixedQual) {
            this.qStart = qStart;
            this.qEnd = qEnd;
            this.nObservations = nObservations;
            this.nErrors = nErrors;
            this.fixedQual = fixedQual;
        }

        public double getErrorRate() {
            if ( nObservations == 0 )
                return 0.0;
            else
                return nErrors / (1.0 * nObservations);
            //return (nErrors+1) / (1.0*(nObservations + 1));
        }

        public byte getQual() {
            if ( fixedQual == -1 )
                return QualityUtils.probToQual(1-getErrorRate(), 0);
            else
                return (byte)fixedQual;
        }

        @Override
        public int compareTo(final QualInterval qualInterval) {
            return new Integer(this.qStart).compareTo(qualInterval.qStart);
        }

        public QualInterval merge(final QualInterval toMerge) {
            final QualInterval left = this.compareTo(toMerge) < 0 ? this : toMerge;
            final QualInterval right = this.compareTo(toMerge) < 0 ? toMerge : this;

            final long nCombinedObs = left.nObservations + right.nObservations;
            final long nCombinedErr = left.nErrors + right.nErrors;

            QualInterval merged = new QualInterval(left.qStart, right.qEnd, nCombinedObs, nCombinedErr );
            merged.subIntervals.add(left);
            merged.subIntervals.add(right);

            return merged;
        }

        public double getPenalty() {
            return calcPenalty(getErrorRate());
        }

        private double calcPenalty(final double globalErrorRate) {
            // the penalty is simply the sum of miscalibrated errors down the tree.
            // Suppose this interval has error e*.
            // the penalty is the sum of (ei - e*) * Ni for all bins covered by this interval

            if ( subIntervals.isEmpty() ) {
                return (Math.abs(getErrorRate() - globalErrorRate)) * nObservations;
            } else {
                double sum = 0;
                for ( QualInterval interval : subIntervals )
                    sum += interval.calcPenalty(globalErrorRate);
                return sum;
            }
        }
    }

    private TreeSet<QualInterval> quantize() {
        // create intervals for each qual individually
        final TreeSet<QualInterval> intervals = new TreeSet<QualInterval>();
        for ( int qStart = 0; qStart < originalSize; qStart++ ) {
            final long nObs = nObservationsPerQual.get(qStart);
            final double errorRate = QualityUtils.qualToErrorProb((byte)qStart);
            final double nErrors = nObs * errorRate;
            final QualInterval qi = new QualInterval(qStart, qStart, nObs, (int)Math.floor(nErrors), (byte)qStart);
            intervals.add(qi);
        }

        // merge intervals for Q0-Q5
        // TODO -- is this even necessary?

        // greedy algorithm:
        // while ( n intervals >= nLevels ):
        //   find intervals to merge with least penalty
        //   merge it
        while ( intervals.size() > nLevels ) {
            mergeLowestPenaltyIntervals(intervals);
        }

        // calculate the penalty for merging
        for ( final QualInterval interval : intervals )
            this.overallPenalty += interval.getPenalty();

        return intervals;
    }

    private void mergeLowestPenaltyIntervals(final TreeSet<QualInterval> intervals) {
        // setup the iterators
        final Iterator<QualInterval> it1 = intervals.iterator();
        final Iterator<QualInterval> it1p = intervals.iterator();
        it1p.next(); // skip one

        QualInterval minMerge = null;
        logger.info("mergeLowestPenaltyIntervals: " + intervals.size());
        while ( it1p.hasNext() ) {
            final QualInterval left = it1.next();
            final QualInterval right = it1p.next();
            final QualInterval merged = left.merge(right);
            if ( minMerge == null || (merged.getPenalty() < minMerge.getPenalty() ) ) {
                logger.info("  Updating merge " + minMerge);
                minMerge = merged;
            }
        }
        logger.info("  => final min merge " + minMerge);
        intervals.removeAll(minMerge.subIntervals);
        intervals.add(minMerge);
        logger.info("updated intervals: " + intervals);
    }

    private List<Byte> intervalsToMap(final TreeSet<QualInterval> intervals) {
        List<Byte> map = new ArrayList<Byte>(originalSize);
        map.addAll(Collections.nCopies(originalSize, Byte.MIN_VALUE));
        for ( QualInterval interval : intervals ) {
            for ( int q = interval.qStart; q <= interval.qEnd; q++ ) {
                map.set(q, interval.getQual());
            }
        }

        if ( Collections.min(map) == Byte.MIN_VALUE )
            throw new ReviewedStingException("quantized quality score map contains an un-initialized value");

        return map;
    }
}
