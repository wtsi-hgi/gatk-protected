package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.jgrapht.graph.DefaultDirectedGraph;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 14, 2011
 */
public class SimpleDeBruijnAssembler extends LocalAssemblyEngine {

    private static final int KMER_OVERLAP = 6; // the additional size of a valid chunk of sequence, used to string together k-mers
    private static final int PRUNE_FACTOR = 1;
    private static final int NUM_BEST_PATHS_PER_KMER = 9;
    
    private final boolean DEBUG;
    private final PrintStream GRAPH_WRITER;

    // the deBruijn graph object
    private final ArrayList<DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>> graphs = new ArrayList<DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>>();

    public SimpleDeBruijnAssembler( final boolean debug, final PrintStream graphWriter ) {
        super();
        DEBUG = debug;
        GRAPH_WRITER = graphWriter;
    }

    public ArrayList<Haplotype> runLocalAssembly(final ArrayList<GATKSAMRecord> reads, final Haplotype refHaplotype) {
        // create the graphs
        createDeBruijnGraphs( reads, refHaplotype );

        /*
        if( GRAPH_WRITER != null ) {
            printGraphs( false );
        }
        */

        // clean up the graph by merging and pruning iteratively
        mergeNodes();
        pruneGraphs();
        mergeNodes();
        pruneGraphs();
        mergeNodes();

        if( GRAPH_WRITER != null ) {
            printGraphs( true );
        }

        // find the best paths in the graphs
        return findBestPaths( refHaplotype );
    }

    private void createDeBruijnGraphs( final ArrayList<GATKSAMRecord> reads, final Haplotype refHaplotype ) {
        graphs.clear();

        // create the graph
        for( int kmer = 7; kmer <= 60; kmer += 8 ) {
            final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
            createGraphFromSequences( graph, reads, kmer, refHaplotype );
            graphs.add(graph);
        }
    }

    private void mergeNodes() {
        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            boolean foundNodesToMerge = true;
            while( foundNodesToMerge ) {
                foundNodesToMerge = false;
                Set<DeBruijnEdge> outEdges = null;
                Set<DeBruijnEdge> inEdges = null;
                final ArrayList<DeBruijnVertex> verticesToRemove = new ArrayList<DeBruijnVertex>();
                DeBruijnVertex incomingVertex = null;
                DeBruijnVertex thisVertex = null;
                for( final DeBruijnEdge e : graph.edgeSet() ) {
                    incomingVertex = graph.getEdgeSource(e);
                    thisVertex = graph.getEdgeTarget(e);
                    if( graph.inDegreeOf(thisVertex) == 1 && graph.outDegreeOf(incomingVertex) == 1) {
                        outEdges = graph.outgoingEdgesOf(thisVertex);
                        inEdges = graph.incomingEdgesOf(incomingVertex);
                        if( edgeSetsHaveSameStatus(e, inEdges, outEdges) ) {
                            foundNodesToMerge = true;
                            verticesToRemove.add(thisVertex);
                            verticesToRemove.add(incomingVertex);
                            break;
                        }
                    }
                }
                if( foundNodesToMerge ) {
                    final byte[] newVertexBases = ArrayUtils.addAll(incomingVertex.getSequence(), thisVertex.getSuffix());
                    final DeBruijnVertex addedVertex = addToGraphIfNew(graph, newVertexBases, thisVertex.kmer);
                    for( final DeBruijnEdge e : outEdges ) {
                        final DeBruijnEdge newEdge = new DeBruijnEdge(e.getIsRef());
                        newEdge.setMultiplicity( e.getMultiplicity() );
                        graph.addEdge(addedVertex, graph.getEdgeTarget(e), newEdge);
                    }
                    for( final DeBruijnEdge e : inEdges ) {
                        final DeBruijnEdge newEdge = new DeBruijnEdge(e.getIsRef());
                        newEdge.setMultiplicity( e.getMultiplicity() );
                        graph.addEdge(graph.getEdgeSource(e), addedVertex, newEdge);
                    }
                    graph.removeAllVertices( verticesToRemove );
                }
            }
        }
    }
    
    private static boolean edgeSetsHaveSameStatus( final DeBruijnEdge edge, final Set<DeBruijnEdge> e1, final Set<DeBruijnEdge> e2 ) {
        final Boolean refStatus = edge.getIsRef();
        final Boolean lowConfStatus = edge.getMultiplicity() > PRUNE_FACTOR;
        if( refStatus && lowConfStatus ) { return true; }
        for( final DeBruijnEdge e : e1 ) {
            if( e.getIsRef() != refStatus ) { return false; }
            if( e.getMultiplicity() > PRUNE_FACTOR != lowConfStatus ) { return false; }
        }
        for( final DeBruijnEdge e : e2 ) {
            if( e.getIsRef() != refStatus ) { return false; }
            if( e.getMultiplicity() > PRUNE_FACTOR != lowConfStatus ) { return false; }
        }
        return true;
    }

    private void pruneGraphs() {
        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            final ArrayList<DeBruijnEdge> edgesToRemove = new ArrayList<DeBruijnEdge>();
            for( final DeBruijnEdge e : graph.edgeSet() ) {
                if( e.getMultiplicity() <= PRUNE_FACTOR && !e.getIsRef() ) {
                    edgesToRemove.add(e);
                }
            }
            graph.removeAllEdges(edgesToRemove);
            final ArrayList<DeBruijnVertex> verticesToRemove = new ArrayList<DeBruijnVertex>();
            for( final DeBruijnVertex v : graph.vertexSet() ) {
                if( graph.inDegreeOf(v) == 0 && graph.outDegreeOf(v) == 0) {
                    verticesToRemove.add(v);
                }
            }
            graph.removeAllVertices(verticesToRemove);
        }
    }

    private static void createGraphFromSequences( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final ArrayList<GATKSAMRecord> reads, final int KMER_LENGTH, final Haplotype refHaplotype ) {
        final byte[] refSequence = refHaplotype.getBases();
        if( refSequence.length > KMER_LENGTH + KMER_OVERLAP ) {
            final int kmersInSequence = refSequence.length - KMER_LENGTH + 1;
            for (int i = 0; i < kmersInSequence - 1; i++) {
                // get the kmers
                final byte[] kmer1 = new byte[KMER_LENGTH];
                System.arraycopy(refSequence, i, kmer1, 0, KMER_LENGTH);
                final byte[] kmer2 = new byte[KMER_LENGTH];
                System.arraycopy(refSequence, i+1, kmer2, 0, KMER_LENGTH);

                addEdgeToGraph( graph, kmer1, kmer2, true );
            }
        }

        for ( final GATKSAMRecord read : reads ) {
            final byte[] sequence = read.getReadBases();
            if( sequence.length > KMER_LENGTH + KMER_OVERLAP ) {
                final int kmersInSequence = sequence.length - KMER_LENGTH + 1;
                for (int i = 0; i < kmersInSequence - 1; i++) {
                    // get the kmers
                    final byte[] kmer1 = new byte[KMER_LENGTH];
                    System.arraycopy(sequence, i, kmer1, 0, KMER_LENGTH);
                    final byte[] kmer2 = new byte[KMER_LENGTH];
                    System.arraycopy(sequence, i+1, kmer2, 0, KMER_LENGTH);

                    addEdgeToGraph( graph, kmer1, kmer2, false );
                }
            }
        }
    }

    private static void addEdgeToGraph( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final byte[] kmer1, final byte[] kmer2, boolean isRef ) {

        final DeBruijnVertex v1 = addToGraphIfNew( graph, kmer1, kmer1.length );
        final DeBruijnVertex v2 = addToGraphIfNew( graph, kmer2, kmer2.length );

        final Set<DeBruijnEdge> edges = graph.outgoingEdgesOf(v1);
        DeBruijnEdge targetEdge = null;
        for ( final DeBruijnEdge edge : edges ) {
            if ( graph.getEdgeTarget(edge).equals(v2) ) {
                targetEdge = edge;
                break;
            }
        }

        if ( targetEdge == null ) {
            graph.addEdge(v1, v2, new DeBruijnEdge( isRef ));
        } else {
            if( isRef ) {
                targetEdge.setIsRef(true);
            } 
            targetEdge.setMultiplicity(targetEdge.getMultiplicity() + 1);
        }
    }

    private static DeBruijnVertex addToGraphIfNew( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final byte[] kmer, final int kmerLength ) {

        // the graph.containsVertex() method is busted, so here's a hack around it
        final DeBruijnVertex newV = new DeBruijnVertex(kmer, kmerLength);
        for ( final DeBruijnVertex v : graph.vertexSet() ) {
            if ( v.equals(newV) )
                return v;
        }

        graph.addVertex(newV);
        return newV;
    }
    
    private void printGraphs( final boolean wasPruned ) {
        int count = 0;
        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            GRAPH_WRITER.println("digraph kmer" + count++ + (wasPruned ? "Pruned":"") +" {");
            for( final DeBruijnEdge edge : graph.edgeSet() ) {
                if( edge.getMultiplicity() > PRUNE_FACTOR ) {
                    GRAPH_WRITER.println("\t" + graph.getEdgeSource(edge).toString() + " -> " + graph.getEdgeTarget(edge).toString() + " [" + (edge.getMultiplicity() <= PRUNE_FACTOR ? "style=dotted,color=grey" : "label=\""+ edge.getMultiplicity() +"\"") + "];");
                }
                if( edge.getIsRef() ) {
                    GRAPH_WRITER.println("\t" + graph.getEdgeSource(edge).toString() + " -> " + graph.getEdgeTarget(edge).toString() + " [color=red];");
                }
                if( !edge.getIsRef() && edge.getMultiplicity() <= PRUNE_FACTOR ) { System.out.println("Graph pruning warning!"); }
            }
            for( final DeBruijnVertex v : graph.vertexSet() ) {
                final String label = ( graph.inDegreeOf(v) == 0 ? v.toString() : v.getSuffixString() );
                GRAPH_WRITER.println("\t" + v.toString() + " [label=\"" + label + "\"]");
            }
            GRAPH_WRITER.println("}");
        }
    }

    private ArrayList<Haplotype> findBestPaths( final Haplotype refHaplotype ) {
        final ArrayList<Haplotype> returnHaplotypes = new ArrayList<Haplotype>();
        returnHaplotypes.add( refHaplotype );

        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            final ArrayList<KBestPaths.Path> bestPaths = KBestPaths.getKBestPaths(graph, NUM_BEST_PATHS_PER_KMER);
            for ( final KBestPaths.Path path : bestPaths ) {
                final Haplotype h = new Haplotype( path.getBases( graph ), path.getScore() );
                if( h.getBases() != null ) {
                    if( h.getBases().length >= refHaplotype.getBases().length ) { // the haplotype needs to be at least as long as the active region
                        if( !returnHaplotypes.contains(h) ) { // no reason to add a new haplotype if the bases are the same as one already present
                            returnHaplotypes.add( h );
                        }
                    }
                }
            }
        }

        if( DEBUG ) { 
            System.out.println("Found " + returnHaplotypes.size() + " candidate haplotypes to evaluate.");
        }

        return returnHaplotypes;
    }
}