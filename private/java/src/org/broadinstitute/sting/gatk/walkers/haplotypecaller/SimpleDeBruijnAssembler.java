package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;
import org.jgrapht.graph.DefaultDirectedGraph;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

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

        // reset the graph
        //graphs = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);

        // clip the reads to get just the base sequences we want
        final ArrayList<byte[]> sequences = clipReads( reads );

        // create the graph
        createDeBruijnGraph( sequences );

        // find the best paths in the graph
        return findBestPaths( refHaplotype );
    }

    // This method takes the base sequences from the SAM records and pulls
    // out runs of bases that are not soft-clipped and are all at least Q20s.
    // Clipped sequences that are overly clipped are not used.
    private ArrayList<byte[]> clipReads(final ArrayList<SAMRecord> reads) {

        final ArrayList<byte[]> sequences = new ArrayList<byte[]>();
        final HashMap<String,SAMRecord> nameMap = new HashMap<String, SAMRecord>();

        for( final SAMRecord read : reads ) {
            if( nameMap.containsKey(read.getReadName()) ) { // assuming only two reads in a pair
                final SAMRecord firstRead = nameMap.remove(read.getReadName());

                // if the reads overlap, note: reads are provided in sorted order by alignment start
                if(read.getAlignmentStart() < firstRead.getAlignmentEnd() && read.getAlignmentStart() > firstRead.getAlignmentStart() && read.getAlignmentEnd() > firstRead.getAlignmentEnd() ) {
                    final int numBases = read.getAlignmentStart() - firstRead.getAlignmentStart() + read.getReadLength();
                    final byte[] bases = new byte[numBases];
                    int iii = 0;
                    for(iii = 0; iii < read.getAlignmentStart() - firstRead.getAlignmentStart(); iii++) {
                        bases[iii] = firstRead.getReadBases()[iii];
                    }
                    for(final Byte b : read.getReadBases()) {
                        bases[iii++] = b;
                    }
                    if( DEBUG ) {
                        System.out.println("Created longer read by merging! " + bases.length);
                        System.out.println(firstRead.getReadString());
                        for(int jjj = 0; jjj < read.getAlignmentStart() - firstRead.getAlignmentStart(); jjj++) {
                            System.out.print(" ");
                        }
                        System.out.println(read.getReadString());
                        String displayString = "";
                        for(int jjj = 0; jjj < bases.length; jjj++) {
                            displayString += (char) bases[jjj];
                        }
                        System.out.println(displayString);
                    }
                    sequences.add( bases );
                } else {
                    sequences.add( read.getReadBases() );
                    sequences.add( firstRead.getReadBases() );
                }

            } else {
                nameMap.put(read.getReadName(), read);
            }
        }

        // actual clipping moved to base walker
        for( final SAMRecord read : nameMap.values() ) {
            //if( read.getReadBases().length > KMER_LENGTH + KMER_OVERLAP ) {
                sequences.add( read.getReadBases() );
            //}
        }

        return sequences;
    }

    private void createDeBruijnGraph(final ArrayList<byte[]> reads) {

        graphs.clear();
        // create the graph
        for( int kmer = 7; kmer <= 101; kmer += 8 ) {
            final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
            createGraphFromSequences( graph, reads, kmer );
            graphs.add(graph);
        }

        // remove nodes with incoming multiplicity of N
        // if ( MIN_MULTIPLICITY_TO_USE > 0 )
        //     removeNodesWithLowMultiplicity();


        // BUGBUG: the merging / cleaning up of nodes doesn't work correctly, need to eventually fix so that graphs can be visualized
        // cleanup graph by merging nodes
        //concatenateNodes();

        // cleanup the node sequences so that they print well
        //cleanupNodeSequences();

        //if ( DEBUG )
        //    printGraph();
    }

    private static void createGraphFromSequences( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final ArrayList<byte[]> reads, final int KMER_LENGTH ) {

        for ( final byte[] sequence : reads ) {
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
    private void concatenateNodes() {

        while ( true ) {
            boolean graphWasModified = false;

            Set<DeBruijnVertex> vertexSet = graph.vertexSet();
            // convert to array because results of the iteration on a set are undefined when the graph is modified
            ArrayList<DeBruijnVertex> vertices = new ArrayList<DeBruijnVertex>(vertexSet);

            for( final DeBruijnVertex v1 : vertices ) {

                // try to merge v1 -> v2
                if ( graph.outDegreeOf(v1) == 1 ) {
                    DeBruijnEdge edge = graph.outgoingEdgesOf(v1).iterator().next();
                    DeBruijnVertex v2 = graph.getEdgeTarget(edge);

                    if ( graph.inDegreeOf(v2) == 1 ) {
                        mergeVertices(v1, v2);
                        graphWasModified = true;
                        break;
                    }
                }

                // try to merge v2 -> v1
                if ( graph.inDegreeOf(v1) == 1 ) {
                    DeBruijnEdge edge = graph.incomingEdgesOf(v1).iterator().next();
                    DeBruijnVertex v2 = graph.getEdgeSource(edge);

                    if ( graph.outDegreeOf(v2) == 1 ) {
                        mergeVertices(v2, v1);
                        graphWasModified = true;
                        break;
                    }
                }
            }

            if ( !graphWasModified )
                break;
        }
    }

    private void mergeVertices(DeBruijnVertex V1, DeBruijnVertex V2) {
        // (Vx -> V1 -> V2 -> Vy)
        //     should now be
        //   (Vx -> V12 -> Vy)

        // create V12
        int additionalSequenceFromV2 = V2.actualSequence.length - KMER_LENGTH + 1;
        byte[] newKmer = new byte[V1.actualSequence.length + additionalSequenceFromV2];
        System.arraycopy(V1.actualSequence, 0, newKmer, 0, V1.actualSequence.length);
        System.arraycopy(V2.actualSequence, KMER_LENGTH - 1, newKmer, V1.actualSequence.length, additionalSequenceFromV2);
        DeBruijnVertex V12 = new DeBruijnVertex(newKmer);
        graph.addVertex(V12);

        // copy edges coming from Vx to V12
        Set<DeBruijnEdge> Ex = graph.incomingEdgesOf(V1);
        for ( DeBruijnEdge edge : Ex ) {
            DeBruijnVertex Vx = graph.getEdgeSource(edge);
            DeBruijnEdge newEdge = new DeBruijnEdge();
            newEdge.setMultiplicity(edge.getMultiplicity());
            graph.addEdge(Vx, V12, newEdge);
        }

        // copy edges going to Vy from V12
        Set<DeBruijnEdge> Ey = graph.outgoingEdgesOf(V2);
        for ( DeBruijnEdge edge : Ey ) {
            DeBruijnVertex Vy = graph.getEdgeTarget(edge);
            DeBruijnEdge newEdge = new DeBruijnEdge();
            newEdge.setMultiplicity(edge.getMultiplicity());
            graph.addEdge(V12, Vy, newEdge);
        }

        // remove V1 and V2 and their associated edges
        graph.removeVertex(V1);
        graph.removeVertex(V2);
    }

    private void cleanupNodeSequences() {

        // remove the first k-1 bases of the kmers
        for ( DeBruijnVertex v :  graph.vertexSet() ) {
            if ( graph.inDegreeOf(v) > 0 )
                v.removePrefix(KMER_LENGTH - 1, true);
        }

        // move common suffixes from incoming nodes to this one

        while ( true ) {

            boolean graphWasModified = false;
            for ( DeBruijnVertex v :  graph.vertexSet() ) {

                if ( graph.inDegreeOf(v) > 1 )  {
                    Set<DeBruijnVertex> connectedVs = new HashSet<DeBruijnVertex>();
                    for ( DeBruijnEdge edge : graph.incomingEdgesOf(v) )
                        connectedVs.add(graph.getEdgeSource(edge));

                    if ( propagateCommonSuffix(v, connectedVs) ) {
                        removeEmptyNodes();
                        graphWasModified = true;
                        break;
                    }
                }
            }

            if ( !graphWasModified )
                break;
        }
    }
    */

    /*
    private void removeEmptyNodes() {

        // remember that results of an iteration on a set are undefined when the graph is modified
        while ( true ) {

            boolean graphWasModified = false;
            for ( DeBruijnVertex v :  graph.vertexSet() ) {
                if ( v.printableSequence.length == 0 ) {
                    removeNode(v);
                    graphWasModified = true;
                    break;
                }
            }

            if ( !graphWasModified )
                break;
        }
    }

    private void removeNode(DeBruijnVertex v) {
        Set<DeBruijnEdge> incoming = graph.incomingEdgesOf(v);
        Set<DeBruijnEdge> outgoing = graph.outgoingEdgesOf(v);

        // make edges from all incoming nodes to all outgoing nodes
        for ( DeBruijnEdge Ex : incoming ) {
            DeBruijnVertex Vx = graph.getEdgeSource(Ex);
            for ( DeBruijnEdge Ey : outgoing ) {
                DeBruijnVertex Vy = graph.getEdgeTarget(Ey);

                DeBruijnEdge newEdge = new DeBruijnEdge();
                newEdge.setMultiplicity(Ex.getMultiplicity());
                graph.addEdge(Vx, Vy, newEdge);
            }
        }

        // remove v and its associated edges
        graph.removeVertex(v);
    }

    private boolean propagateCommonSuffix(DeBruijnVertex Vx, Set<DeBruijnVertex> incoming) {

        // find the common matching suffix
        byte[] match = null;
        for ( DeBruijnVertex v : incoming ) {
            if ( match == null ) {
                match = v.printableSequence;
            } else {
                int idx = 0;
                while ( idx < match.length && idx < v.printableSequence.length && match[match.length - idx - 1] == v.printableSequence[v.printableSequence.length - idx - 1] )
                    idx++;

                if ( idx < match.length ) {
                    match = new byte[idx];
                    System.arraycopy(v.printableSequence, v.printableSequence.length - idx, match, 0, idx);
                }
            }
        }

        // if there is a common suffix...
        if ( match != null && match.length > 0 ) {

            // remove the suffix from the end of the incoming nodes...
            for ( DeBruijnVertex v : incoming )
                v.removeSuffix(match.length, false);

            // ...and put it at the front of this node
            Vx.addPrefix(match, false);
            return true;
        }

        return false;
    }

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
            final ArrayList<KBestPaths.Path> bestPaths = KBestPaths.getKBestPaths(graph, 10);

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

    private void assignReadsToGraph(ArrayList<byte[]> reads) {

        // TODO -- implement me

    }

    /****
    private void removeNodesWithLowMultiplicity() {

        Set<DeBruijnVertex> vertexSet = graph.vertexSet();
        // convert to array because results of the iteration on a set are undefined when the graph is modified
        ArrayList<DeBruijnVertex> vertices = new ArrayList<DeBruijnVertex>(vertexSet);

        for (int i = 0; i < vertices.size(); i++) {

            DeBruijnVertex v = vertices.get(i);
            if ( graph.inDegreeOf(v) == 1 &&
                    graph.incomingEdgesOf(v).iterator().next().getMultiplicity() < MIN_MULTIPLICITY_TO_USE )
                removeNode(v);
        }
    }
    ****/
}
