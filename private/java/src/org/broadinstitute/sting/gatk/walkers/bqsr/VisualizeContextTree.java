/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
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
    public String treeOutputPrefix = "";

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

    @Argument(fullName = "eventType", shortName = "eventType", doc="", required = false)
    public List<String> eventTypes = Arrays.asList("I", "M", "D");

    @Argument(fullName = "maxPTarget", shortName = "mpt", doc="The maximum P-value to allow when pruning to Pvalue, as phred-scaled value (i.e., 30 means 10^-3)", required = false)
    public List<Integer> maxPTargets = Arrays.asList(30);

    @Argument(fullName = "numNodesTarget", shortName = "nnt", doc="The maximum number of nodes to allow in the tree when pruning to num nodes", required = false)
    public List<Integer> numNodesTargets = Arrays.asList(64);

    @Argument(fullName = "maxDepthTarget", shortName = "mdt", doc="The maximum depth allowed in the tree when pruning to depth", required = false)
    public List<Integer> maxDepthTargets = Arrays.asList(3);

    @Argument(fullName = "mode", shortName = "mode", doc="", required = false)
    public Mode MODE = Mode.ANALYZE;

    @Argument(fullName = "test", shortName = "test", doc="", required = false)
    public boolean TEST = false;

    @Argument(fullName = "qualForMVisualization", shortName = "qualForMVisualization", doc="", required = false)
    public int QUAL_FOR_M_VISUALIZATION = 30;

    @Argument(fullName = "writePartialContexts", shortName = "writePartialContexts", doc="", required = false)
    public boolean WRITE_PARTIAL_CONTEXTS = false;

    @Argument(fullName = "nContext", shortName = "nContext", doc="", required = false)
    public boolean DEBUG_CONTEXT_WITH_N = false;

    public enum Mode {
        ANALYZE,
        UPDATE
    }

    public enum Operation {
        KEEP_ORIGINAL,
        PRUNE_TO_DEPTH,
        PRUNE_TO_NUM_NODES_BY_PENALTY,
        PRUNE_TO_PENALTY
    }

    private List<Integer> getPruneTargets(final Operation operation) {
        switch ( operation ) {
            case KEEP_ORIGINAL: return Arrays.asList(0); // doesn't matter
            case PRUNE_TO_DEPTH: return maxDepthTargets;
            case PRUNE_TO_NUM_NODES_BY_PENALTY: return numNodesTargets;
            case PRUNE_TO_PENALTY: return maxPTargets;
            default: throw new IllegalArgumentException("Expected operation " + operation);
        }
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
        Set<RecalDataSubset> all = new HashSet<RecalDataSubset>();
        for ( int i = 0; i < optionalCovariates.getNumRows(); i++ ) {
            if ( optionalCovariates.get(i, "CovariateName").equals("Context") ) {
                final String rg = (String)optionalCovariates.get(i, "ReadGroup");
                final int qual = Integer.valueOf((String)optionalCovariates.get(i, "QualityScore"));
                final String eventType = (String)optionalCovariates.get(i, "EventType");
                all.add(new RecalDataSubset(rg, qual, eventType));
            }
        }

        // sort all by transforming it into a treeset
        all = new TreeSet<RecalDataSubset>(all);

        for ( final RecalDataSubset qe : all ) {
            logger.info("Found QE " + qe);
        }

        if ( TEST )
            return Collections.singleton(new RecalDataSubset(all.iterator().next().rg, 45, "I"));
        else
            return all;
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
                toKeep.add( new ContextDatum(context, (int)observations, (double)errors) );
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

        logger.info("Reading recalibration report " + RECAL_FILE);
        final GATKReport recalibrationReport = new GATKReport(RECAL_FILE);
        final GATKReportTable optionalCovariates = recalibrationReport.getTable("RecalTable2");
        logger.info(" ... done");

        final AdaptiveTreeAnalysis analysisReport = new AdaptiveTreeAnalysis();
        if ( MODE == Mode.ANALYZE) {
            for ( final String eventType : eventTypes) {
                final RecalDataSubset selectInfo = new RecalDataSubset(eventType);
                final List<ContextDatum> ourContexts = subsetToOurContexts(optionalCovariates, selectInfo);
                final RecalDatumNode<ContextDatum> root = AdaptiveContext.createTreeFromFlatContexts(ourContexts);

                RecalDatumNode<ContextDatum> initialTree = filterContexts(root);

                analyzeAdaptiveTrees(selectInfo, analysisReport, initialTree);
            }
        } else if ( MODE == Mode.UPDATE ) {
            if ( operations.size() != 1 )
                throw new UserException.BadArgumentValue("operations", "Only one operation allowed with UPDATE mode");

            final Operation operation = operations.iterator().next();
            final List<Integer> pruneTargets = getPruneTargets(operation);
            if ( pruneTargets.size() != 1 )
                throw new UserException.BadArgumentValue("pruneTarget", "Only one pruneTarget allowed with UPDATE mode");
            final int pruneTarget = pruneTargets.get(0);

            final GATKReport updatedReport = startUpdatedReport(recalibrationReport);
            for ( final RecalDataSubset selectInfo : getQualAndEventTypesForContext(optionalCovariates) ) {
                logger.info("Processing recal data " + selectInfo);
                final List<ContextDatum> ourContext = subsetToOurContexts(optionalCovariates, selectInfo);
                final RecalDatumNode<ContextDatum> root = AdaptiveContext.createTreeFromFlatContexts(ourContext);
                final RecalDatumNode<ContextDatum> pruned = applyOpToTree(operation, pruneTarget, root);

                // write out the analysis of the pruning as requested
                analysisReport.add(selectInfo.qual, selectInfo.eventType, operation, pruneTarget, pruned);

                // and for comparison purposes the pruning by depth
                final RecalDatumNode<ContextDatum> prunedByDepth = applyOpToTree(Operation.PRUNE_TO_DEPTH, pruneTarget, root);
                if ( prunedByDepth != null )
                    analysisReport.add(selectInfo.qual, selectInfo.eventType, Operation.PRUNE_TO_DEPTH, pruneTarget, prunedByDepth);

                // write out the full tree to the new BQSR report
                final RecalDatumNode<ContextDatum> fullTree =
                        WRITE_PARTIAL_CONTEXTS ? pruned : AdaptiveContext.fillToDepth(pruned, root.maxDepth() - 1, DEBUG_CONTEXT_WITH_N);
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
    private class RecalDataSubset implements Comparable<RecalDataSubset> {
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
        public int compareTo(RecalDataSubset recalDataSubset) {
            int cmp = Integer.valueOf(qual).compareTo(recalDataSubset.qual);

            if ( cmp == 0 )
                cmp = eventType.compareTo(recalDataSubset.eventType);

            if ( cmp == 0 )
                cmp = rg.compareTo(recalDataSubset.eventType);

            return cmp;
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

    private RecalDatumNode<ContextDatum> applyOpToTree(final Operation operation,
                                                       final int pruneTarget,
                                                       final RecalDatumNode<ContextDatum> initialTree) {
        switch (operation) {
            case KEEP_ORIGINAL:
                return initialTree;
            case PRUNE_TO_DEPTH:
                return initialTree.pruneToDepth(pruneTarget + 1); // add one for the root
            case PRUNE_TO_NUM_NODES_BY_PENALTY:
                return initialTree.pruneByPenalty(pruneTarget+1); // add one for the root
            case PRUNE_TO_PENALTY:
                return initialTree.pruneToNoMoreThanPenalty(pruneTarget, true);
            default:
                throw new IllegalArgumentException("Unexpected operation");
        }
    }

    /**
     * loop over the prune targets, collecting information in analysisReport
     */
    private void analyzeAdaptiveTrees(final RecalDataSubset selectInfo,
                                      final AdaptiveTreeAnalysis analysisReport,
                                      final RecalDatumNode<ContextDatum> initialTree) {

        for ( final Operation operation : operations ) {
            for ( final int pruneTarget : getPruneTargets(operation) ) {
                logger.info("Analyzing " + selectInfo + " with pruneType " + pruneTarget + " and operation " + operation);
                RecalDatumNode<ContextDatum> prunedTree = applyOpToTree(operation, pruneTarget, initialTree);

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

        final File dest = new File(treeOutputPrefix + outputName + ".dot");
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
