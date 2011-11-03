package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;
import org.jgrapht.graph.DefaultDirectedGraph;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 14, 2011
 */
public class SimpleDeBruijnAssembler extends LocalAssemblyEngine {

    private static final boolean DEBUG = false;

    // the additional size of a valid chunk of sequence, used to string together k-mers
    private static final int KMER_OVERLAP = 6;

    // the deBruijn graph object
    private final ArrayList<DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>> graphs = new ArrayList<DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>>();

    public SimpleDeBruijnAssembler(PrintStream out, IndexedFastaSequenceFile referenceReader) {
        super(out, referenceReader);
    }

    public ArrayList<Haplotype> runLocalAssembly(final ArrayList<SAMRecord> reads, final Haplotype refHaplotype) {
        // create the graphs
        createDeBruijnGraph( reads );

        // find the best paths in the graphs
        return findBestPaths( refHaplotype );
    }

    private void createDeBruijnGraph(final ArrayList<SAMRecord> reads) {
        graphs.clear();
        // create the graph
        for( int kmer = 7; kmer <= 101; kmer += 8 ) {
            final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
            createGraphFromSequences( graph, reads, kmer );
            graphs.add(graph);
        }
    }

    private static void createGraphFromSequences( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final ArrayList<SAMRecord> reads, final int KMER_LENGTH ) {
        for ( final SAMRecord read : reads ) {
            final byte[] sequence = read.getReadBases();
            if( sequence.length > KMER_LENGTH + KMER_OVERLAP ) {
                final int kmersInSequence = sequence.length - KMER_LENGTH + 1;
                for (int i = 0; i < kmersInSequence - 1; i++) {
                    // get the kmers
                    final byte[] kmer1 = new byte[KMER_LENGTH];
                    System.arraycopy(sequence, i, kmer1, 0, KMER_LENGTH);
                    final byte[] kmer2 = new byte[KMER_LENGTH];
                    System.arraycopy(sequence, i+1, kmer2, 0, KMER_LENGTH);

                    addEdgeToGraph( graph, kmer1, kmer2 );

                    // TODO -- eventually, we'll need to deal with reverse complementing the sequences ???
                }
            }
        }
    }

    private static void addEdgeToGraph( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final byte[] kmer1, final byte[] kmer2 ) {

        final DeBruijnVertex v1 = addToGraphIfNew( graph, kmer1 );
        final DeBruijnVertex v2 = addToGraphIfNew( graph, kmer2 );

        final Set<DeBruijnEdge> edges = graph.outgoingEdgesOf(v1);
        DeBruijnEdge targetEdge = null;
        for ( DeBruijnEdge edge : edges ) {
            if ( graph.getEdgeTarget(edge).equals(v2) ) {
                targetEdge = edge;
                break;
            }
        }

        if ( targetEdge == null )
            graph.addEdge(v1, v2, new DeBruijnEdge());
        else
            targetEdge.setMultiplicity(targetEdge.getMultiplicity() + 1);
    }

    private static DeBruijnVertex addToGraphIfNew( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final byte[] kmer ) {

        // the graph.containsVertex() method is busted, so here's a hack around it
        final DeBruijnVertex newV = new DeBruijnVertex(kmer);
        for ( final DeBruijnVertex v : graph.vertexSet() ) {
            if ( v.equals(newV) )
                return v;
        }

        graph.addVertex(newV);
        return newV;
    }

    /*
    private void printGraph() {

        if( getOutputStream() != null ) {
            for ( DeBruijnVertex source : graph.vertexSet() ) {
                if ( graph.inDegreeOf(source) == 0 )
                    getOutputStream().print("* ");
                getOutputStream().print(source + " -> ");
                for ( DeBruijnEdge edge : graph.outgoingEdgesOf(source) ) {
                    getOutputStream().print(graph.getEdgeTarget(edge) + " (" + edge.getMultiplicity() + "), ");
                }
                getOutputStream().println();
            }
            getOutputStream().println("------------\n");
        }
    }
    */

    private ArrayList<Haplotype> findBestPaths( final Haplotype refHaplotype ) {
        final ArrayList<Haplotype> returnHaplotypes = new ArrayList<Haplotype>();
        returnHaplotypes.add( refHaplotype );

        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            final ArrayList<KBestPaths.Path> bestPaths = KBestPaths.getKBestPaths(graph, 13);

            for ( final KBestPaths.Path path : bestPaths ) {
                final Haplotype h = new Haplotype( path.getBases( graph ), path.getScore() );
                if( h.bases != null ) {
                    if( getOutputStream() != null ) {
                        getOutputStream().println(h.toString());
                    }
                    if( !returnHaplotypes.contains(h) ) {
                        returnHaplotypes.add( h );
                    }
                }
            }
        }

        return returnHaplotypes;
    }

}
