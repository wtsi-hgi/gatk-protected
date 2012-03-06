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

import net.sf.samtools.SAMReadGroupRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.PrintStream;
import java.util.*;

/**
 * A general algorithm for quantizing quality score distributions to use a specific number of levels
 *
 * @author Mark Depristo
 * @since 3/2/12
 */
public class QualQuantizer {
    private static Logger logger = Logger.getLogger(QualQuantizer.class);

    final int nLevels, originalSize, minInterestingQual;
    final List<Long> nObservationsPerQual;

    /** Map from original qual (e.g., Q30) to new quantized qual (e.g., Q28) */
    final List<Byte> originalToQuantizedMap;
    final List<Double> penaltyPerQual;
    final TreeSet<QualInterval> quantizedIntervals;
    double overallPenalty;

    protected QualQuantizer() {
        this(Collections.<Long>emptyList(), 0, 0);
        // for testing purposes only
    }

    public QualQuantizer(final List<Long> nObservationsPerQual, final int nLevels, final int minInterestingQual) {
        this.nObservationsPerQual = nObservationsPerQual;
        this.nLevels = nLevels;
        this.minInterestingQual = minInterestingQual;

        this.originalSize = nObservationsPerQual.size();
        this.penaltyPerQual = new ArrayList<Double>(originalSize);

        // for QC.  These values should never appear after a correct run
        Collections.fill(penaltyPerQual, -1.0);

        this.quantizedIntervals = quantize();

        this.originalToQuantizedMap = intervalsToMap(quantizedIntervals);
    }

    public class QualInterval implements Comparable<QualInterval> {
        final int qStart, qEnd;
        final long nObservations;
        final long nErrors;
        final int fixedQual;
        final int level;
        int mergeOrder;
        Set<QualInterval> subIntervals = new HashSet<QualInterval>();

        public QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors, final int level) {
            this(qStart, qEnd, nObservations, nErrors, level, -1);
        }

        public QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors, final int level, final int fixedQual) {
            this.qStart = qStart;
            this.qEnd = qEnd;
            this.nObservations = nObservations;
            this.nErrors = nErrors;
            this.fixedQual = fixedQual;
            this.level = level;
            this.mergeOrder = 0;
        }

        public String getName() {
            return qStart + "-" + qEnd;
        }

        @Override
        public String toString() {
            return "QQ:" + getName();
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

            final int level = Math.max(left.level, right.level) + 1;
            QualInterval merged = new QualInterval(left.qStart, right.qEnd, nCombinedObs, nCombinedErr, level );
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
                if ( this.qEnd <= minInterestingQual )
                    // It's free to merge up quality scores below the smallest interesting one
                    return 0;
                else {
                    if ( nErrors == 0 )
                        return 0;
                    else {
                        return (Math.abs(Math.log10(getErrorRate()) - Math.log10(globalErrorRate))) * nObservations;
                    }
                }
                // this is the linear error rate penalty
                //return (Math.abs(getErrorRate() - globalErrorRate)) * nObservations;
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
            final QualInterval qi = new QualInterval(qStart, qStart, nObs, (int)Math.floor(nErrors), 0, (byte)qStart);
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
        int lastMergeOrder = 0;
        while ( it1p.hasNext() ) {
            final QualInterval left = it1.next();
            final QualInterval right = it1p.next();
            final QualInterval merged = left.merge(right);
            lastMergeOrder = Math.max(Math.max(lastMergeOrder, left.mergeOrder), right.mergeOrder);
            if ( minMerge == null || (merged.getPenalty() < minMerge.getPenalty() ) ) {
                logger.info("  Updating merge " + minMerge);
                minMerge = merged;
            }
        }
        logger.info("  => final min merge " + minMerge);
        intervals.removeAll(minMerge.subIntervals);
        intervals.add(minMerge);
        minMerge.mergeOrder = lastMergeOrder + 1;
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

    public void writeReport(PrintStream out) {
        final GATKReport report = new GATKReport();

        addQualHistogramToReport(report);
        addIntervalsToReport(report);

        report.print(out);
    }

    private final void addQualHistogramToReport(final GATKReport report) {
        report.addTable("QualHistogram", "Quality score histogram provided to report");
        GATKReportTable table = report.getTable("QualHistogram");

        table.addPrimaryKey("qual");
        table.addColumn("count", "NA");

        for ( int q = 0; q < nObservationsPerQual.size(); q++ ) {
            table.set(q, "count", nObservationsPerQual.get(q));
        }
    }


    private final void addIntervalsToReport(final GATKReport report) {
        report.addTable("QualQuantizerIntervals", "Table of QualQuantizer quantization intervals");
        GATKReportTable table = report.getTable("QualQuantizerIntervals");

        table.addPrimaryKey("name");
        table.addColumn("qStart", "NA");
        table.addColumn("qEnd", "NA");
        table.addColumn("level", "NA");
        table.addColumn("merge.order", "NA");
        table.addColumn("nErrors", "NA");
        table.addColumn("nObservations", "NA");
        table.addColumn("qual", "NA");
        table.addColumn("penalty", "NA");
        table.addColumn("root.node", "NA");
        //table.addColumn("subintervals", "NA");

        for ( QualInterval interval : quantizedIntervals)
            addIntervalToReport(table, interval, true);
    }

    private final void addIntervalToReport(final GATKReportTable table, QualInterval interval, final boolean atRootP) {
        final String name = interval.getName();
        table.set(name, "qStart", interval.qStart);
        table.set(name, "qEnd", interval.qEnd);
        table.set(name, "level", interval.level);
        table.set(name, "merge.order", interval.mergeOrder);
        table.set(name, "nErrors", interval.nErrors);
        table.set(name, "nObservations", interval.nObservations);
        table.set(name, "qual", interval.getQual());
        table.set(name, "penalty", String.format("%.1f", interval.getPenalty()));
        table.set(name, "root.node", atRootP);

        for ( final QualInterval sub : interval.subIntervals )
            addIntervalToReport(table, sub, false);
    }

    public List<Byte> getOriginalToQuantizedMap() {
        return originalToQuantizedMap;
    }

    public TreeSet<QualInterval> getQuantizedIntervals() {
        return quantizedIntervals;
    }

    public double getOverallPenalty() {
        return overallPenalty;
    }
}
