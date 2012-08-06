package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.utils.recalibration.ContextDatum;
import org.broadinstitute.sting.utils.recalibration.RecalDatumNode;

import java.io.PrintStream;

/**
 * Allows one to create and write out analyses of adaptive context trees
 *
 * User: depristo
 * Date: 8/3/12
 * Time: 9:30 AM
 * To change this template use File | Settings | File Templates.
 */
public final class AdaptiveTreeAnalysis {
    private final GATKReport report;

    public AdaptiveTreeAnalysis() {
        report = GATKReport.newSimpleReport("AdaptiveContextAnalysis", "QualityScore", "eventType", "operation",
                "pruneTarget", "size", "numLeafs", "minDepth",
                "maxDepth", "penalty", "maxPenalty", "minPenalty");
    }

    /**
     * Add information about pruned tree to this analysis report
     *
     * @param contextType the context type (I, M, D) corresponding to the tree
     * @param operation the operation we applied to the tree
     * @param pruneTarget the target value used during pruning
     * @param prunedTree the pruned tree itself
     */
    public void add(final int qual,
                    final String contextType,
                    final VisualizeContextTree.Operation operation,
                    final int pruneTarget,
                    final RecalDatumNode<ContextDatum> prunedTree) {
        report.addRow(qual, contextType, operation, pruneTarget, prunedTree.size(), prunedTree.numLeaves(),
                prunedTree.minDepth(), prunedTree.maxDepth(),
                prunedTree.totalPenalty(), prunedTree.maxPenalty(), prunedTree.minPenalty());
    }


    public void print(final PrintStream out) {
        report.print(out);
    }
}
