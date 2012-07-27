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
import org.broadinstitute.sting.utils.recalibration.RecalDatum;
import org.jgrapht.DirectedGraph;
import org.jgrapht.ext.*;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleDirectedGraph;

import java.io.File;
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

    @Output(doc="Write output to this BAM filename instead of STDOUT")
    PrintStream out;

    @Argument(fullName = "recalFile", doc="", required = true)
    public File RECAL_FILE;

    @Argument(fullName = "Qual", shortName = "Q", doc="Quality scores to visualize", required = false)
    public int qualToTake = 45;

    @Argument(fullName = "maxDepth", shortName = "maxDepth", doc="Quality scores to visualize", required = false)
    public int maxDepth = 4;

    @Argument(fullName = "prefix", shortName = "prefix", doc="", required = false)
    public String prefix = null;

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

    final class ContextDatum extends RecalDatum {
        final String context;

        ContextDatum(final String context, final long observations, final long errors ) {
            super(observations, errors, (byte)qualToTake);
            this.context = context;
        }

        @Override
        public String toString() {
            return context;
        }

        /**
         * Are the context bases of this object a prefix of the those in cd2?
         *
         * @param cd2
         * @return
         */
        public boolean isParentContext(final ContextDatum cd2) {
            return (cd2.context.startsWith(context) && cd2.size() - 1 == size()) ||
                    (cd2.size() == 1 && size() == 0);
        }

        public String getParentContext() {
            return context.substring(0, size() - 1);
        }

        public int size() { return context.length(); }
    }

    final static class ContextDataLabelProvider implements VertexNameProvider<ContextDatum> {
        @Override
        public String getVertexName(final ContextDatum contextDatum) {
            return String.format("%s:Q%d:N%d",
                    contextDatum.context,
                    (int)contextDatum.getEmpiricalQuality(),
                    (int)(Math.log10(contextDatum.getNumObservations()) * 10));
        }
    }

    final static class ContextDataAttributeProvider implements ComponentAttributeProvider<ContextDatum> {
        @Override
        public Map<String, String> getComponentAttributes(final ContextDatum contextDatum) {
            final Map<String, String> components = new HashMap<String, String>();

            final double Q = contextDatum.getEmpiricalQuality();
            final double Qscale = Math.min(Math.max(Q - 10, 0), 60) / 60.0; // from 0.0 -> 1.0
            final String color = heatmapColors(Qscale);
            components.put("color", color);
            components.put("fontcolor", color);
            components.put("penwidth", "2.0");
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

    private final List<ContextDatum> subsetToOurContexts(final GATKReportTable optionalCovariates) {
        final List<ContextDatum> toKeep = new ArrayList<ContextDatum>();

        for ( Object rowKey : optionalCovariates.getRowIDs() ) {
            if ( keepRow(rowKey, optionalCovariates) ) {
                final String context = (String)optionalCovariates.get(rowKey, "CovariateValue");
                final long observations = (Long)optionalCovariates.get(rowKey, "Observations");
                final long errors = (Long)optionalCovariates.get(rowKey, "Errors");
                toKeep.add( new ContextDatum(context, observations, errors) );
            }
        }

        return toKeep;
    }

    private final boolean keepRow( final Object rowKey, final GATKReportTable optionalCovariates) {
        final List<String> columnsToCheck = Arrays.asList("QualityScore", "CovariateName", "EventType");
        final List<Object> valuesToMatch = Arrays.asList((Object)Integer.toString(this.qualToTake), "Context", "I");
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
        final List<ContextDatum> ourContexts = subsetToOurContexts(optionalCovariates);
        final List<ContextDatum> allCombinations = createAllSubcontexts(ourContexts);
        final List<ContextDatum> filteredCombinations = filterContexts(allCombinations);
        final DirectedGraph<ContextDatum, DefaultEdge> contextTree = buildContextTree(filteredCombinations);

        // create the visualizer, and write out the DOT graph for us
        final VertexNameProvider<ContextDatum> vertexIDProvider = new StringNameProvider<ContextDatum>();
        final VertexNameProvider<ContextDatum> vertexLabelProvider = new ContextDataLabelProvider();
        final EdgeNameProvider<ContextDatum> edgeNameProvider = null; // new StringEdgeNameProvider();
        final ComponentAttributeProvider<ContextDatum> vertexAttributeProvider = new ContextDataAttributeProvider();
        final DOTExporter exporter = new DOTExporter(vertexIDProvider, vertexLabelProvider, edgeNameProvider, vertexAttributeProvider, null);

        exporter.export(new PrintWriter(out), contextTree);
    }

    private final List<ContextDatum> filterContexts(final List<ContextDatum> contexts) {
        final List<ContextDatum> filtered = new ArrayList<ContextDatum>(contexts.size());

        for ( final ContextDatum cd : contexts )
            if ( cd.size() <= maxDepth &&
                 (prefix == null || cd.context.startsWith(prefix)) )
                filtered.add(cd);

        return filtered;
    }

    private final List<ContextDatum> createAllSubcontexts(final List<ContextDatum> nContexts) {
        final int n = nContexts.get(0).size();

        if ( n == 0 ) // recursive case -- we are done
            return Arrays.asList(new ContextDatum("x", 0, 0));
        else {
            final List<ContextDatum> oneLayerUp = createSubcontextsOnNextLayer(nContexts, n);
            final List<ContextDatum> allUp = createAllSubcontexts(oneLayerUp);

            final List<ContextDatum> all = new ArrayList<ContextDatum>(nContexts);
            all.addAll(allUp);
            return all;
        }
    }

    private final List<ContextDatum> createSubcontextsOnNextLayer(final List<ContextDatum> contexts, int contextSize) {
        Map<String, ContextDatum> up = new HashMap<String, ContextDatum>(contexts.size());

        for ( final ContextDatum cd : contexts ) {
            final String parentContext = cd.getParentContext();
            ContextDatum parent = up.get(parentContext);
            if ( parent == null ) {
                parent = new ContextDatum(parentContext, 0, 0);
                up.put(parentContext, parent);
            }

            parent.incrementNumObservations(cd.getNumObservations());
            parent.incrementNumMismatches(cd.getNumMismatches());
        }

        return new ArrayList<ContextDatum>(up.values());
    }

    private final DirectedGraph<ContextDatum, DefaultEdge> buildContextTree(final List<ContextDatum> ourContexts) {
        DirectedGraph<ContextDatum, DefaultEdge> g = new SimpleDirectedGraph<ContextDatum, org.jgrapht.graph.DefaultEdge>(DefaultEdge.class);

        for ( final ContextDatum cd : ourContexts ) {
            // add the vertices
            g.addVertex(cd);

            for ( final ContextDatum cd2 : ourContexts ) {
                if ( cd.isParentContext(cd2) )
                    g.addEdge(cd, cd2);
            }
        }

        return g;
    }
}
