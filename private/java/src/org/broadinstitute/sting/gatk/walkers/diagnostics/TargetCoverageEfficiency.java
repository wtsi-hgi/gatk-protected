package org.broadinstitute.sting.gatk.walkers.diagnostics;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Reports the ratio of coverage between the target (middle of the interval) and the average distribution of coverages in the interval.
 * <p/>
 * <p>
 * For each interval, in an interval file, it computes the target coverage (target is defined as the middle of the interval), the average
 * coverage of the interval and the ratio between the target and the interval coverage.
 *
 * This tool is designed to evaluate the penalty payed to achieve the target coverage
 * </p>
 * <p/>
 * <p/>
 * <h2>Input</h2>
 * <p>
 * A BAM file and an interval list with the targeted sequencing
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * A table with the values per interval for : target coverage, average interval coverage and the ratio between the two.
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T TargetCoverageEfficiency
 *      -R reference.fasta
 *      -I input.bam
 *      -L intervals.list
 *  </pre>
 *
 * @author Mauricio Carneiro
 * @since 11/10/11
 */
@By(DataSource.REFERENCE)
public class TargetCoverageEfficiency extends LocusWalker<Long, ArrayList<Long>> {
    @Output
    PrintStream out;

    GATKReportTable reportTable;

    public void initialize () {
        reportTable = new GATKReportTable("TargetCoverageEfficiency", "A table with the values per interval for: target coverage, average interval coverage and the ratio between the two.", true);
        reportTable.addPrimaryKey("Interval", true);
        reportTable.addColumn("Target", 0, true);
        reportTable.addColumn("Average", 0, true);
        reportTable.addColumn("Ratio", 0, true);
    }

    public boolean isReduceByInterval () {
        return true;
    }

    @Override
    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return (long) context.size();
    }

    @Override
    public ArrayList<Long> reduceInit() {
        return new ArrayList<Long>();
    }

    @Override
    public ArrayList<Long> reduce(Long value, ArrayList<Long> sum) {
        sum.add(value);
        return sum;
    }

    @Override
    public void onTraversalDone(List<Pair<GenomeLoc, ArrayList<Long>>> results) {
        for (Pair<GenomeLoc, ArrayList<Long>> intervalPair : results) {
            GenomeLoc interval = intervalPair.getFirst();
            List<Long> distribution = intervalPair.getSecond();

            int targetIndex = (int) Math.floor((interval.getStart() + interval.getStop()) / 2) - interval.getStart();    // Our target will be in the middle of the interval
            long targetCoverage = distribution.get(targetIndex);
            double averageCoverage = MathUtils.average(distribution);
            double ratio = (averageCoverage == 0) ? 0 : targetCoverage / averageCoverage;

            reportTable.set(interval.toString(), "Target", targetCoverage);
            reportTable.set(interval.toString(), "Average", averageCoverage);
            reportTable.set(interval.toString(), "Ratio", ratio);
        }
        reportTable.write(out);
    }
}
