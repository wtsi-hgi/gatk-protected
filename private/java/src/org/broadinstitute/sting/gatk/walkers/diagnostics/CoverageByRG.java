package org.broadinstitute.sting.gatk.walkers.diagnostics;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 12/26/11
 * Time: 11:01 PM
 * To change this template use File | Settings | File Templates.
 */
public class CoverageByRG extends LocusWalker<List<Pair<String, Long>>, List<Pair<String, Long>>> {

    @Output
    PrintStream out;

    Long currentIntervalSize;
    GATKReportTable reportTable;

    public void initialize() {
        reportTable = new GATKReportTable("CoverageByRG", "A table with the values per interval for each read group", true);
        reportTable.addPrimaryKey("Interval", true);
    }

    public boolean isReduceByInterval() {
        return true;
    }

    @Override
    public List<Pair<String, Long>> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        ReadBackedPileup pileup = context.getBasePileup();

        // If the table has no columns, add them
        if (reportTable.getColumns().size() < 1) {
            for (String RG : pileup.getReadGroups()) {
                reportTable.addColumn(RG, 0, true);
            }
        }

        List<Pair<String, Long>> output = new LinkedList<Pair<String, Long>>();

        for (String RG : pileup.getReadGroups()) {
            output.add(new Pair<String, Long>(RG, (long) pileup.getPileupForReadGroup(RG).depthOfCoverage()));
        }
        return output;
    }

    @Override
    public List<Pair<String, Long>> reduceInit() {
        currentIntervalSize = 0L;
        return new LinkedList<Pair<String, Long>>();

    }

    @Override
    public List<Pair<String, Long>> reduce(List<Pair<String, Long>> value, List<Pair<String, Long>> sum) {
        if (sum.isEmpty()) {
            currentIntervalSize = 1L;
            return value;
        }

        for (Pair<String, Long> valuePair : value) {
            // Find the sum with the same RG
            Iterator<Pair<String, Long>> i = sum.iterator();
            Pair<String, Long> sumPair = i.next();

            while (i.hasNext()) {
                if (valuePair.getFirst().equals(sumPair.getFirst())) {
                    sumPair.second += valuePair.getSecond();
                    break;
                }
                sumPair = i.next();
            }
        }

        currentIntervalSize++;
        return sum;
    }

    @Override
    public void onTraversalDone(List<Pair<GenomeLoc, List<Pair<String, Long>>>> results) {
        for (Pair<GenomeLoc, List<Pair<String, Long>>> intervalPair : results) {
            GenomeLoc interval = intervalPair.getFirst();
            List<Pair<String, Long>> counts = intervalPair.getSecond();

            for (Pair<String, Long> sumPair : counts) {
                reportTable.set(interval.toString(), sumPair.getFirst(), (double) sumPair.getSecond() / (double) interval.size());
            }
            currentIntervalSize = 0L;
        }

        reportTable.write(out);

    }
}
