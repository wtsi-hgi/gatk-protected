package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.jgrapht.graph.DefaultDirectedGraph;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks, rpoplin
 * Date: Mar 14, 2011
 */

public class SimpleDeBruijnAssembler extends LocalAssemblyEngine {

    private static final int KMER_OVERLAP = 6; // the additional size of a valid chunk of sequence, used to string together k-mers
    private static final int NUM_BEST_PATHS_PER_KMER_GRAPH = 7;
    private static final byte MIN_QUALITY = (byte) 10;

    private final boolean DEBUG;
    private final PrintStream GRAPH_WRITER;
    private final ArrayList<DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>> graphs = new ArrayList<DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>>();

    private int PRUNE_FACTOR = 1;
    
    public SimpleDeBruijnAssembler( final boolean debug, final PrintStream graphWriter ) {
        super();
        DEBUG = debug;
        GRAPH_WRITER = graphWriter;
    }

    public ArrayList<Haplotype> runLocalAssembly( final ArrayList<GATKSAMRecord> reads, final Haplotype refHaplotype, final int PRUNE_FACTOR ) {
        this.PRUNE_FACTOR = PRUNE_FACTOR;

        // create the graphs
        createDeBruijnGraphs( reads, refHaplotype );

        // clean up the graphs by pruning and merging
        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            pruneGraph( graph );
            eliminateNonRefPaths( graph );
            mergeNodes( graph );
        }

        if( GRAPH_WRITER != null ) {
            printGraphs();
        }

        // find the best paths in the graphs
        return findBestPaths( refHaplotype );
    }

    private void createDeBruijnGraphs( final ArrayList<GATKSAMRecord> reads, final Haplotype refHaplotype ) {
        graphs.clear();

        // create the graph
        for( int kmer = 7; kmer <= 85; kmer += 6 ) {
            final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
            if( createGraphFromSequences( graph, reads, kmer, refHaplotype ) ) {
                graphs.add(graph);
            }
        }
    }

    protected static void mergeNodes( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph ) {
        boolean foundNodesToMerge = true;
        while( foundNodesToMerge ) {
            foundNodesToMerge = false;
            final ArrayList<DeBruijnVertex> verticesToRemove = new ArrayList<DeBruijnVertex>();
            for( final DeBruijnEdge e : graph.edgeSet() ) {
                final DeBruijnVertex outgoingVertex = graph.getEdgeTarget(e);
                final DeBruijnVertex incomingVertex = graph.getEdgeSource(e);
                if( !outgoingVertex.equals(incomingVertex) && graph.inDegreeOf(outgoingVertex) == 1 && graph.outDegreeOf(incomingVertex) == 1) {
                    final Set<DeBruijnEdge> outEdges = graph.outgoingEdgesOf(outgoingVertex);
                    final Set<DeBruijnEdge> inEdges = graph.incomingEdgesOf(incomingVertex);
                    foundNodesToMerge = true;
                    verticesToRemove.add(outgoingVertex);
                    verticesToRemove.add(incomingVertex);
                    if( inEdges.size() == 1 && outEdges.size() == 1 ) {
                        inEdges.iterator().next().setMultiplicity( inEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() / 2 ) );
                        outEdges.iterator().next().setMultiplicity( outEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() / 2 ) );
                    } else if( inEdges.size() == 1 ) {
                        inEdges.iterator().next().setMultiplicity( inEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() - 1 ) );
                    } else if( outEdges.size() == 1 ) {
                        outEdges.iterator().next().setMultiplicity( outEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() - 1 ) );
                    }

                    final byte[] newVertexBases = ArrayUtils.addAll(incomingVertex.getSequence(), outgoingVertex.getSuffix());
                    final DeBruijnVertex addedVertex = new DeBruijnVertex( newVertexBases, outgoingVertex.kmer );
                    graph.addVertex(addedVertex);
                    for( final DeBruijnEdge edge : outEdges ) {
                        final DeBruijnEdge newEdge = new DeBruijnEdge(edge.getIsRef());
                        newEdge.setMultiplicity( edge.getMultiplicity() );
                        graph.addEdge(addedVertex, graph.getEdgeTarget(edge), newEdge);
                    }
                    for( final DeBruijnEdge edge : inEdges ) {
                        final DeBruijnEdge newEdge = new DeBruijnEdge(edge.getIsRef());
                        newEdge.setMultiplicity( edge.getMultiplicity() );
                        graph.addEdge(graph.getEdgeSource(edge), addedVertex, newEdge);
                    }
                    break;
                }
            }
            graph.removeAllVertices( verticesToRemove );
        }
    }

    protected void pruneGraph( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph ) {
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

    protected static void eliminateNonRefPaths( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph ) {
        final ArrayList<DeBruijnVertex> verticesToRemove = new ArrayList<DeBruijnVertex>();
        boolean done = false;
        while( !done ) {
            done = true;
            for( final DeBruijnVertex v : graph.vertexSet() ) {
                if( graph.inDegreeOf(v) == 0 || graph.outDegreeOf(v) == 0 ) {
                    boolean isRefNode = false;
                    for( final DeBruijnEdge e : graph.edgesOf(v) ) {
                        if( e.getIsRef() ) {
                            isRefNode = true;
                            break;
                        }
                    }
                    if( !isRefNode ) {
                        done = false;
                        verticesToRemove.add(v);
                    }
                }
            }
            graph.removeAllVertices(verticesToRemove);
            verticesToRemove.clear();
        }
    }

    private static boolean createGraphFromSequences( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final ArrayList<GATKSAMRecord> reads, final int KMER_LENGTH, final Haplotype refHaplotype ) {
        final byte[] refSequence = refHaplotype.getBases();
        if( refSequence.length > KMER_LENGTH + KMER_OVERLAP ) {
            final int kmersInSequence = refSequence.length - KMER_LENGTH + 1;
            for (int i = 0; i < kmersInSequence - 1; i++) {
                // get the kmers
                final byte[] kmer1 = new byte[KMER_LENGTH];
                System.arraycopy(refSequence, i, kmer1, 0, KMER_LENGTH);
                final byte[] kmer2 = new byte[KMER_LENGTH];
                System.arraycopy(refSequence, i+1, kmer2, 0, KMER_LENGTH);
                if( !addKmersToGraph(graph, kmer1, kmer2, true) ) { return false; }
            }
        }

        for( final GATKSAMRecord read : reads ) {
            final byte[] sequence = read.getReadBases();
            final byte[] qualities = read.getBaseQualities();
            if( sequence.length > KMER_LENGTH + KMER_OVERLAP ) {
                final int kmersInSequence = sequence.length - KMER_LENGTH + 1;
                for( int iii = 0; iii < kmersInSequence - 1; iii++ ) {                    
                    // if the qualities of all the bases in the kmers are high enough
                    boolean badKmer = false;
                    for( int jjj = iii; jjj < iii + KMER_LENGTH + 1; jjj++) {
                        if( qualities[jjj] < MIN_QUALITY ) {
                            badKmer = true;
                            break;
                        }
                    }
                    if( !badKmer ) {
                        // get the kmers
                        final byte[] kmer1 = new byte[KMER_LENGTH];
                        System.arraycopy(sequence, iii, kmer1, 0, KMER_LENGTH);
                        final byte[] kmer2 = new byte[KMER_LENGTH];
                        System.arraycopy(sequence, iii+1, kmer2, 0, KMER_LENGTH);

                        addKmersToGraph(graph, kmer1, kmer2, false);
                    }
                }
            }
        }
        return true;
    }

    protected static boolean addKmersToGraph( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final byte[] kmer1, final byte[] kmer2, final boolean isRef ) {

        final int numVertexBefore = graph.vertexSet().size();
        final DeBruijnVertex v1 = new DeBruijnVertex( kmer1, kmer1.length );
        graph.addVertex(v1);
        final DeBruijnVertex v2 = new DeBruijnVertex( kmer2, kmer2.length );
        graph.addVertex(v2);
        if( isRef && graph.vertexSet().size() == numVertexBefore ) { return false; }

        final DeBruijnEdge targetEdge = graph.getEdge(v1, v2);
        if ( targetEdge == null ) {
            graph.addEdge(v1, v2, new DeBruijnEdge( isRef ));
        } else {
            if( isRef ) {
                targetEdge.setIsRef( true );
            }
            targetEdge.setMultiplicity(targetEdge.getMultiplicity() + 1);
        }
        return true;
    }

    private void printGraphs() {
        int count = 0;
        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            GRAPH_WRITER.println("digraph kmer" + count++ +" {");
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

    @Ensures({"result.contains(refHaplotype)"})
    private ArrayList<Haplotype> findBestPaths( final Haplotype refHaplotype ) {
        final ArrayList<Haplotype> returnHaplotypes = new ArrayList<Haplotype>();
        returnHaplotypes.add( refHaplotype );

        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            for ( final KBestPaths.Path path : KBestPaths.getKBestPaths(graph, NUM_BEST_PATHS_PER_KMER_GRAPH) ) {
                final Haplotype h = new Haplotype( path.getBases( graph ), path.getScore() );
                if( !returnHaplotypes.contains(h) ) { // no reason to add a new haplotype if the bases are the same as one already present
                    returnHaplotypes.add( h );
                }
            }
        }

        if( DEBUG ) { 
            if( returnHaplotypes.size() > 1 ) {
                System.out.println("Found " + returnHaplotypes.size() + " candidate haplotypes to evaluate every read against.");
            } else {
                System.out.println("Found only the reference haplotype in the assembly graph.");
            }
        }

        return returnHaplotypes;
    }
}