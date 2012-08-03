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
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.recalibration.RecalDatum;
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
 * Consumes a BQSR table and creates a graphviz visualization of the context covariate tree
 *
 * <p>
 * x
 *
 * <h2>Input</h2>
 * <p>
 * One BQSR recalibration GATK report.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A dot file visualizing the BQSR context coverage.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T ContextTree \
 *   -o output.dot \
 *   -recalData bqsr.dat.gatkreport.txt
 * </pre>
 *
 */
public class VisualizeContextTree extends RefWalker<Integer, Integer> {
    private final static Logger logger = Logger.getLogger(VisualizeContextTree.class);

    @Output(doc="Write analysis output to this file")
    PrintStream out;

    @Argument(fullName = "treeOutPrefix", shortName = "treeOutPrefix", doc="Write tree output to files with this prefix", required = true)
    public String treeOutputPrefix;

    @Argument(fullName = "recalFile", doc="", required = true)
    public File RECAL_FILE;

    @Argument(fullName = "maxDepth", shortName = "maxDepth", doc="Quality scores to visualize", required = false)
    public int maxDepth = 4;

    @Argument(fullName = "prefix", shortName = "prefix", doc="", required = false)
    public String prefix = null;

    @Argument(fullName = "operations", shortName = "ops", doc="", required = false)
    public Set<Operation> operations = EnumSet.allOf(Operation.class);

    @Argument(fullName = "contextType", shortName = "C", doc="", required = false)
    public List<String> contextTypes = Arrays.asList("I", "M", "D");

    @Argument(fullName = "pruneTarget", shortName = "pt", doc="", required = false)
    public List<Integer> pruneTargets = Arrays.asList(4);

    public enum Operation {
        KEEP_ORIGINAL,
        PRUNE_BY_DEPTH,
        PRUNE_BY_PENALTY
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

    public final static String ROOT_CONTEXT = "x";
    final class ContextDatum extends RecalDatum {
        final String context;

        public ContextDatum(final String context, final long observations, final long errors ) {
            super(observations, errors, (byte)30); // TODO -- should use default value?
            this.context = context;
        }

        @Override
        public String toString() {
            return context;
        }

        public String getParentContext() {
            return size() == 1 ? ROOT_CONTEXT : context.substring(0, size() - 1);
        }

        public int size() { return context.length(); }
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

    private final static String heatmapColors(final double value) {
        final double H = value;
        final double S = 0.75;
        final double V = 1.0;
        return String.format("%.3f %.3f %.3f", H, S, V);
        //final int colorValue = (int)(Qscale * 255);
        //return String.format("#%2x%2x%2x", (int)R, (int)G, (int)B);
    }

    private final List<ContextDatum> subsetToOurContexts(final GATKReportTable optionalCovariates, final String contextType) {
        final List<ContextDatum> toKeep = new ArrayList<ContextDatum>();

        for ( Object rowKey : optionalCovariates.getRowIDs() ) {
            if ( keepRow(rowKey, optionalCovariates, contextType) ) {
                final String context = (String)optionalCovariates.get(rowKey, "CovariateValue");
                final long observations = (Long)optionalCovariates.get(rowKey, "Observations");
                final long errors = (Long)optionalCovariates.get(rowKey, "Errors");
                toKeep.add( new ContextDatum(context, observations, errors) );
            }
        }

        return toKeep;
    }

    private final boolean keepRow( final Object rowKey, final GATKReportTable optionalCovariates, final String contextType) {
        final List<String> columnsToCheck = Arrays.asList("QualityScore", "CovariateName", "EventType");
        final int qualToTake = contextType.equals("M") ? 30 : 45; // TODO fixme!
        final List<Object> valuesToMatch = Arrays.asList((Object)Integer.toString(qualToTake), "Context", contextType);
        for ( int i = 0; i < columnsToCheck.size(); i++ ) {
            final Object actual = optionalCovariates.get(rowKey, columnsToCheck.get(i));
            final Object expected = valuesToMatch.get(i);
            if ( ! expected.equals(actual) )
                return false;
        }
        return true;
    }

    @Override
    public void onTraversalDone(final Integer result) {
        final GATKReport recalibrationReport = new GATKReport(RECAL_FILE);
        final GATKReportTable optionalCovariates = recalibrationReport.getTable("RecalTable2");

        for ( final String contextType : contextTypes ) {
            final List<ContextDatum> ourContexts = subsetToOurContexts(optionalCovariates, contextType);
            final RecalDatumNode<ContextDatum> root = createAllSubcontexts(ourContexts);

            RecalDatumNode<ContextDatum> initialTree = filterContexts(root);

            analyzeAdaptiveTrees(initialTree, contextType);
        }
    }

    /**
     * loop over the prune targets, collecting information and writing out the analysis output
     *
     * @param initialTree
     */
    private void analyzeAdaptiveTrees(final RecalDatumNode<ContextDatum> initialTree, final String contextType) {
        GATKReport report = GATKReport.newSimpleReport("AdaptiveContextAnalysis", "contextType", "operation", "pruneTarget", "size", "numLeafs", "minDepth", "maxDepth", "penalty");

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
                    report.addRow(contextType, operation, pruneTarget, prunedTree.size(), prunedTree.numLeaves(), prunedTree.minDepth(), prunedTree.maxDepth(), prunedTree.totalPenalty());
                    final String name = String.format("contextType_%s.prune_op_%s.prune_target_%d", contextType, operation, pruneTarget);
                    visualizeTree(name, prunedTree);
                }
            }
        }

        report.print(out);
    }

    private final void visualizeTree(final String outputName, final RecalDatumNode<ContextDatum> root) {
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

    private final RecalDatumNode<ContextDatum> filterContexts(final RecalDatumNode<ContextDatum> node) {
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

    // ---------------------------------------------------------------------------
    //
    // Go from a flat list of contexts to the full tree
    //
    // ---------------------------------------------------------------------------

    private final RecalDatumNode<ContextDatum> createAllSubcontexts(final List<ContextDatum> nContexts) {
        final Queue<RecalDatumNode<ContextDatum>> remaining = new LinkedList<RecalDatumNode<ContextDatum>>();
        final Map<String, RecalDatumNode<ContextDatum>> contextToNodes = new HashMap<String, RecalDatumNode<ContextDatum>>();
        RecalDatumNode<ContextDatum> root = null;

        // initialize -- start with all of the contexts
        for ( final ContextDatum cd : nContexts )
            remaining.add(new RecalDatumNode<ContextDatum>(cd));

        while ( remaining.peek() != null ) {
            final RecalDatumNode<ContextDatum> add = remaining.poll();
            final ContextDatum cd = add.getRecalDatum();

            final String parentContext = cd.getParentContext();
            RecalDatumNode<ContextDatum> parent = contextToNodes.get(parentContext);
            if ( parent == null ) {
                // haven't yet found parent, so make one, and enqueue it for processing
                parent = new RecalDatumNode<ContextDatum>(new ContextDatum(parentContext, 0, 0));
                contextToNodes.put(parentContext, parent);

                if ( parentContext != ROOT_CONTEXT )
                    remaining.add(parent);
                else
                    root = parent;
            }

            parent.getRecalDatum().incrementNumObservations(cd.getNumObservations());
            parent.getRecalDatum().incrementNumMismatches(cd.getNumMismatches());
            parent.addSubnode(add);
        }

        if ( root == null )
            throw new RuntimeException("root is unexpectedly null");

        // set the fixed penalty everywhere in the tree, so that future modifications
        // don't
        root.calcAndSetFixedPenalty(true);

        return root;
    }


    // ---------------------------------------------------------------------------
    //
    // Functions to visualize trees
    //
    // ---------------------------------------------------------------------------

    private final DirectedGraph<RecalDatumNode<ContextDatum>, DefaultEdge> buildContextTree(final RecalDatumNode<ContextDatum> root) {
        final DirectedGraph<RecalDatumNode<ContextDatum>, DefaultEdge> g = new SimpleDirectedGraph<RecalDatumNode<ContextDatum>, org.jgrapht.graph.DefaultEdge>(DefaultEdge.class);
        buildContextTreeRec(g, root);
        return g;
    }

    private final void buildContextTreeRec(final DirectedGraph<RecalDatumNode<ContextDatum>, DefaultEdge> g, final RecalDatumNode<ContextDatum> root) {
        // add the vertices
        g.addVertex(root);

        for ( final RecalDatumNode<ContextDatum> sub : root.getSubnodes() ) {
            g.addVertex(sub);
            g.addEdge(root, sub);
            buildContextTreeRec(g, sub);
        }
    }
}
