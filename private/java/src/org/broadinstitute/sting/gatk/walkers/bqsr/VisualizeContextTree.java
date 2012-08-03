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

package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.recalibration.AdaptiveContext;
import org.broadinstitute.sting.utils.recalibration.ContextDatum;
import org.broadinstitute.sting.utils.recalibration.RecalDatumNode;
import org.jgrapht.DirectedGraph;
import org.jgrapht.ext.*;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleDirectedGraph;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.*;

/**
 * Consumes a BQSR table and creates a graphviz visualization of the context covariate tree or produces
 * an adaptive context determined version of the same BQSR table
 *
 * TODO -- split into base functionality and determine the best place to put it, when / if we
 * TODO -- decide that adaptive contexts are the way to go in the future.
 */
public class VisualizeContextTree extends RefWalker<Integer, Integer> {
    private final static Logger logger = Logger.getLogger(VisualizeContextTree.class);

    @Output(doc="Write analysis output to this file")
    PrintStream out;

    @Argument(fullName = "treeOutPrefix", shortName = "treeOutPrefix", doc="Write tree output to files with this prefix", required = false)
    public String treeOutputPrefix;

    @Input(fullName = "recalFile", doc="", required = true)
    public File RECAL_FILE;

    @Argument(fullName = "recalFileOut", doc="", required = false)
    public File RECAL_FILE_OUT = null;

    @Argument(fullName = "maxDepth", shortName = "maxDepth", doc="Quality scores to visualize", required = false)
    public int maxDepth = 4;

    @Argument(fullName = "prefix", shortName = "prefix", doc="", required = false)
    public String prefix = "";

    @Argument(fullName = "operations", shortName = "ops", doc="", required = false)
    public Set<Operation> operations = EnumSet.allOf(Operation.class);

    @Argument(fullName = "contextType", shortName = "C", doc="", required = false)
    public List<String> contextTypes = Arrays.asList("I", "M", "D");

    @Argument(fullName = "pruneTarget", shortName = "pt", doc="", required = false)
    public List<Integer> pruneTargets = Arrays.asList(4);

    @Argument(fullName = "mode", shortName = "mode", doc="", required = false)
    public Mode MODE = Mode.ANALYZE;

    @Argument(fullName = "test", shortName = "test", doc="", required = false)
    public boolean TEST = false;

    @Argument(fullName = "qualForMVisualization", shortName = "qualForMVisualization", doc="", required = false)
    public int QUAL_FOR_M_VISUALIZATION = 30;


    @Argument(fullName = "writePartialContexts", shortName = "writePartialContexts", doc="", required = false)
    public boolean WRITE_PARTIAL_CONTEXTS = false;

    public enum Mode {
        ANALYZE,
        UPDATE
    }

    public enum Operation {
        KEEP_ORIGINAL,
        PRUNE_BY_DEPTH,
        PRUNE_BY_PENALTY
    }

    @Override
    public void initialize() {
        if ( RECAL_FILE_OUT == null )
            RECAL_FILE_OUT = new File(RECAL_FILE.getAbsolutePath() + ".adaptive.grp");
    }

    @Override
    public boolean isDone() {
        return true;
    }

    @Override
    public Integer map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.return null;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce( Integer read, Integer output ) {
        return 0;
    }

    /**
     * Walk over the BQSR covariates table and return a set of all of the distinct RecalDataSubsets
     * we may want to build adaptive trees for.
     *
     * @param optionalCovariates
     * @return
     */
    private Set<RecalDataSubset> getQualAndEventTypesForContext(final GATKReportTable optionalCovariates) {
        final Set<RecalDataSubset> all = new HashSet<RecalDataSubset>();
        for ( int i = 0; i < optionalCovariates.getNumRows(); i++ ) {
            final String rg = (String)optionalCovariates.get(i, "ReadGroup");
            final int qual = Integer.valueOf((String)optionalCovariates.get(i, "QualityScore"));
            final String eventType = (String)optionalCovariates.get(i, "EventType");
            all.add(new RecalDataSubset(rg, qual, eventType));
        }

        for ( final RecalDataSubset qe : all ) {
            logger.info("Found QE " + qe);
        }

        return TEST ? Collections.singleton(new RecalDataSubset(all.iterator().next().rg, 45, "I")) : all;
    }

    /**
     * Pull out the subset of contexts matches the requested RecalDataSubset information
     * @param optionalCovariates
     * @param selectInfo
     * @return
     */
    private List<ContextDatum> subsetToOurContexts(final GATKReportTable optionalCovariates, final RecalDataSubset selectInfo) {
        final List<ContextDatum> toKeep = new ArrayList<ContextDatum>();

        for ( Object rowKey : optionalCovariates.getRowIDs() ) {
            if ( keepRow(rowKey, optionalCovariates, selectInfo) ) {
                final String context = (String)optionalCovariates.get(rowKey, "CovariateValue");
                final long observations = (Long)optionalCovariates.get(rowKey, "Observations");
                final long errors = (Long)optionalCovariates.get(rowKey, "Errors");
                toKeep.add( new ContextDatum(context, observations, errors) );
            }
        }

        return toKeep;
    }

    private RecalDatumNode<ContextDatum> filterContexts(final RecalDatumNode<ContextDatum> node) {
        return node;

        // TODO -- fixme
//        final ContextDatum datum = node.getRecalDatum();
//        if ( datum.size() <= maxDepth && (prefix == null || datum.context.startsWith(prefix)) ) {
//            final Set<RecalDatumNode<ContextDatum>> filteredSubs = new HashSet<RecalDatumNode<ContextDatum>>();
//            for ( final RecalDatumNode<ContextDatum> sub : node.getSubnodes() ) {
//                final RecalDatumNode<ContextDatum> filteredSub = filterContexts(sub);
//                if ( filteredSub != null )
//                    filteredSubs.add(filteredSub);
//            }
//            return new RecalDatumNode<ContextDatum>(datum, filteredSubs);
//        }
//        else
//            return null;
    }


    private boolean keepRow( final Object rowKey, final GATKReportTable optionalCovariates, final RecalDataSubset selectInfo) {
        final List<String> columnsToCheck = Arrays.asList("ReadGroup", "QualityScore", "CovariateName", "EventType");
        final List<Object> valuesToMatch = Arrays.asList(selectInfo.rg, (Object)Integer.toString(selectInfo.qual), "Context", selectInfo.eventType);
        for ( int i = 0; i < columnsToCheck.size(); i++ ) {
            final Object actual = optionalCovariates.get(rowKey, columnsToCheck.get(i));
            final Object expected = valuesToMatch.get(i);
            if ( ! expected.equals("ignored") && ! expected.equals(actual) )
                return false;
        }
        return true;
    }

    @Override
    public void onTraversalDone(final Integer result) {
        // TODO -- split into two walkers

        final GATKReport recalibrationReport = new GATKReport(RECAL_FILE);
        final GATKReportTable optionalCovariates = recalibrationReport.getTable("RecalTable2");

        final AdaptiveTreeAnalysis analysisReport = new AdaptiveTreeAnalysis();
        if ( MODE == Mode.ANALYZE) {
            for ( final String contextType : contextTypes ) {
                final RecalDataSubset selectInfo = new RecalDataSubset(contextType);
                final List<ContextDatum> ourContexts = subsetToOurContexts(optionalCovariates, selectInfo);
                final RecalDatumNode<ContextDatum> root = AdaptiveContext.createTreeFromFlatContexts(ourContexts);

                RecalDatumNode<ContextDatum> initialTree = filterContexts(root);

                analyzeAdaptiveTrees(selectInfo, analysisReport, initialTree);
            }
        } else if ( MODE == Mode.UPDATE ) {
            if ( pruneTargets.size() != 1 )
                throw new UserException.BadArgumentValue("pruneTarget", "Only one pruneTarget allowed with UPDATE mode");

            final GATKReport updatedReport = startUpdatedReport(recalibrationReport);
            final int pruneTarget = pruneTargets.get(0);
            for ( final RecalDataSubset selectInfo : getQualAndEventTypesForContext(optionalCovariates) ) {
                logger.info("Processing recal data " + selectInfo);
                final List<ContextDatum> ourContext = subsetToOurContexts(optionalCovariates, selectInfo);
                final RecalDatumNode<ContextDatum> root = AdaptiveContext.createTreeFromFlatContexts(ourContext);
                final RecalDatumNode<ContextDatum> pruned = root.pruneByPenalty(pruneTarget+1);

                final RecalDatumNode<ContextDatum> fullTree =
                        WRITE_PARTIAL_CONTEXTS ? pruned : AdaptiveContext.fillToDepth(pruned, root.maxDepth() - 1);

                analysisReport.add(selectInfo.qual, selectInfo.eventType, Operation.PRUNE_BY_PENALTY, pruneTarget, pruned);
                addContextsToReport(updatedReport, selectInfo, fullTree);
            }

            try {
                final PrintStream recalOut = new PrintStream(RECAL_FILE_OUT);
                updatedReport.print(recalOut);
                recalOut.close();
            } catch ( FileNotFoundException e ) {
                throw new UserException.CouldNotCreateOutputFile(RECAL_FILE_OUT, e);
            }

        } else {
            throw new UserException("Unexpected mode " + MODE);
        }

        analysisReport.print(out);
    }

    /**
     * All of the other information in the BSQR context table associated with the context covariates
     * we've turned into a tree.  Needed when writing the adaptive tree back out again.
     */
    private class RecalDataSubset {
        final String rg;
        final int qual;
        final String eventType;

        private RecalDataSubset(final String eventType) {
            rg = "ignored";
            qual = eventType.equals("M") ? QUAL_FOR_M_VISUALIZATION : 45; // TODO fixme!
            this.eventType = eventType;
        }

        private RecalDataSubset(final String readGroup, int qual, String eventType) {
            this.rg = readGroup;
            this.qual = qual;
            this.eventType = eventType;
        }


        @Override
        public String toString() {
            return String.format("qual=%d eventType=%s", qual, eventType);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof RecalDataSubset)) return false;

            RecalDataSubset that = (RecalDataSubset) o;

            if (qual != that.qual) return false;
            if (eventType != null ? !eventType.equals(that.eventType) : that.eventType != null) return false;
            if (rg != null ? !rg.equals(that.rg) : that.rg != null) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = rg != null ? rg.hashCode() : 0;
            result = 31 * result + qual;
            result = 31 * result + (eventType != null ? eventType.hashCode() : 0);
            return result;
        }
    }

    /**
     * Take the original input BQSR GATK report and return a newly allocated
     * minus the context covariates (as these will be added later after the adaptive context algorithms
     * have run)
     *
     * @param originalReport GATKReport from BQSR output
     * @return
     */
    private GATKReport startUpdatedReport(final GATKReport originalReport) {
        final GATKReport updated = new GATKReport();

        // most of the tables in the recal report we just want to copy verbatim
        final List<String> tablesToCopy = Arrays.asList("Arguments", "Quantized", "RecalTable0", "RecalTable1");
        for ( final String tableToCopy : tablesToCopy )
            updated.addTable(originalReport.getTable(tableToCopy));

        // we need to create a new table for RecalTable2 containing all Covariate data except the contexts
        final GATKReportTable originalRecalTable2 = originalReport.getTable("RecalTable2");
        final GATKReportTable updatedRecalTable2 = new GATKReportTable(originalRecalTable2, false);

        int newI = 0;
        for ( int i = 0; i < originalRecalTable2.getNumRows(); i++ ) {
            if ( ! originalRecalTable2.get(i, "CovariateName").equals("Context") ) {
                for ( int j = 0; j < originalRecalTable2.getNumColumns(); j++ )
                    updatedRecalTable2.set(newI, j, originalRecalTable2.get(i, j));
                newI++;
            }
        }
        updated.addTable(updatedRecalTable2);

        return updated;
    }

    private void addContextsToReport(final GATKReport updatedReport,
                                     final RecalDataSubset info,
                                     final RecalDatumNode<ContextDatum> fullTree) {
        final GATKReportTable covariatesTable = updatedReport.getTable("RecalTable2");
        for ( final ContextDatum context : fullTree.getAllLeaves() ) {
            final int rowIndex = covariatesTable.getNumRows();
            covariatesTable.set(rowIndex, "ReadGroup", info.rg);
            covariatesTable.set(rowIndex, "QualityScore", String.valueOf(info.qual));
            covariatesTable.set(rowIndex, "CovariateValue", context.context);
            covariatesTable.set(rowIndex, "CovariateName", "Context");
            covariatesTable.set(rowIndex, "EventType", info.eventType);
            covariatesTable.set(rowIndex, "EmpiricalQuality", String.valueOf(context.getEmpiricalQuality()));
            covariatesTable.set(rowIndex, "Observations", String.valueOf(context.getNumObservations()));
            covariatesTable.set(rowIndex, "Errors", String.valueOf(context.getNumMismatches()));
        }
    }

    /**
     * loop over the prune targets, collecting information in analysisReport
     */
    private void analyzeAdaptiveTrees(final RecalDataSubset selectInfo,
                                      final AdaptiveTreeAnalysis analysisReport,
                                      final RecalDatumNode<ContextDatum> initialTree) {

        for ( final int pruneTarget : pruneTargets ) {
            for ( final Operation operation : operations ) {
                RecalDatumNode<ContextDatum> prunedTree = null;
                switch (operation) {
                    case KEEP_ORIGINAL:
                        prunedTree = initialTree;
                        break;
                    case PRUNE_BY_DEPTH:
                        final double potentialDepth = Math.log10(pruneTarget) / Math.log10(4);
                        if ( Math.floor(potentialDepth) == potentialDepth ) { // we are divisable by 4
                            final int depth = (int)Math.floor(potentialDepth) + 1;
                            prunedTree = initialTree.pruneToDepth(depth);
                        }
                        break;
                    case PRUNE_BY_PENALTY:
                        prunedTree = initialTree.pruneByPenalty(pruneTarget+1);
                        break;
                    default:
                        throw new IllegalArgumentException("Unexpected operation");
                }

                if ( prunedTree != null ) {
                    analysisReport.add(selectInfo.qual, selectInfo.eventType, operation, pruneTarget, prunedTree);
                    final String name = String.format("contextType_%s.qual_%d.prune_op_%s.prune_target_%d",
                            selectInfo.eventType, selectInfo.qual, operation, pruneTarget);
                    visualizeTree(name, prunedTree);
                }
            }
        }
    }

    // ---------------------------------------------------------------------------
    //
    // Functions to visualize trees
    //
    // ---------------------------------------------------------------------------

    /**
     * Write out a DOT graph representation of tree root
     *
     * @param outputName
     * @param root
     */
    private void visualizeTree(final String outputName, final RecalDatumNode<ContextDatum> root) {
        final DirectedGraph<RecalDatumNode<ContextDatum>, DefaultEdge> contextGraph = buildContextTree(root);

        // create the visualizer, and write out the DOT graph for us
        final VertexNameProvider<RecalDatumNode<ContextDatum>> vertexIDProvider = new StringNameProvider<RecalDatumNode<ContextDatum>>();
        final VertexNameProvider<RecalDatumNode<ContextDatum>> vertexLabelProvider = new ContextDataLabelProvider();
        final EdgeNameProvider<RecalDatumNode<ContextDatum>> edgeNameProvider = null; // new StringEdgeNameProvider();
        final ComponentAttributeProvider<RecalDatumNode<ContextDatum>> vertexAttributeProvider = new ContextDataAttributeProvider();
        final DOTExporter exporter = new DOTExporter(vertexIDProvider, vertexLabelProvider, edgeNameProvider, vertexAttributeProvider, null);

        final File dest = new File(treeOutputPrefix + "." + outputName + ".dot");
        try {
            final PrintWriter out = new PrintWriter(new PrintStream(dest));
            exporter.export(out, contextGraph);
            out.close();
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(dest, e);
        }
    }

    final static class ContextDataLabelProvider implements VertexNameProvider<RecalDatumNode<ContextDatum>> {
        @Override
        public String getVertexName(final RecalDatumNode<ContextDatum> contextDatum) {
            return String.format("%s:Q%d:N%d:P%.2e",
                    contextDatum.getRecalDatum().context,
                    (int)contextDatum.getRecalDatum().getEmpiricalQuality(),
                    (int)(Math.log10(contextDatum.getRecalDatum().getNumObservations()) * 10),
                    contextDatum.getPenalty());
        }
    }

    final static class ContextDataAttributeProvider implements ComponentAttributeProvider<RecalDatumNode<ContextDatum>> {
        @Override
        public Map<String, String> getComponentAttributes(final RecalDatumNode<ContextDatum> contextDatum) {
            final Map<String, String> components = new HashMap<String, String>();

            final double Q = contextDatum.getRecalDatum().getEmpiricalQuality();
            final double Qscale = Math.min(Math.max(Q - 10, 0), 60) / 60.0; // from 0.0 -> 1.0
            final String color = heatmapColors(Qscale);
            components.put("color", color);
            components.put("fontcolor", color);
            components.put("penwidth", "0.5");
            components.put("fontsize", "10.0");
            components.put("shape", "none");
            return components;
        }
    }

    private static String heatmapColors(final double value) {
        final double H = value;
        final double S = 0.75;
        final double V = 1.0;
        return String.format("%.3f %.3f %.3f", H, S, V);
        //final int colorValue = (int)(Qscale * 255);
        //return String.format("#%2x%2x%2x", (int)R, (int)G, (int)B);
    }

    /**
     * Create a graph representing the tree root
     *
     * @param root
     * @return
     */
    private DirectedGraph<RecalDatumNode<ContextDatum>, DefaultEdge> buildContextTree(final RecalDatumNode<ContextDatum> root) {
        final DirectedGraph<RecalDatumNode<ContextDatum>, DefaultEdge> g = new SimpleDirectedGraph<RecalDatumNode<ContextDatum>, org.jgrapht.graph.DefaultEdge>(DefaultEdge.class);
        buildContextTreeRec(g, root);
        return g;
    }

    /**
     * Recursive helper routine for #buildContextTree
     * @param g
     * @param root
     */
    private void buildContextTreeRec(final DirectedGraph<RecalDatumNode<ContextDatum>, DefaultEdge> g, final RecalDatumNode<ContextDatum> root) {
        // add the vertices
        g.addVertex(root);

        for ( final RecalDatumNode<ContextDatum> sub : root.getSubnodes() ) {
            g.addVertex(sub);
            g.addEdge(root, sub);
            buildContextTreeRec(g, sub);
        }
    }
}
