package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.jgrapht.graph.DefaultDirectedGraph;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Collections;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 23, 2011
 */
// Class for finding the K best paths (as determined by the sum of multiplicities of the edges) in a graph.
// This is different from most graph traversals because we want to test paths from any source node to any sink node.
public class KBestPaths {

    // static access only
    protected KBestPaths() { }

    protected static class MyInt { public int val = 0;}

    // class to keep track of paths
    protected static class Path {

        // the last vertex seen in the path
        private DeBruijnVertex lastVertex;

        // the list of edges comprising the path
        private List<DeBruijnEdge> edges;

        // the scores for the path
        private int totalScore = 0, lowestEdge = -1;

        public Path(DeBruijnVertex initialVertex) {
            lastVertex = initialVertex;
            edges = new ArrayList<DeBruijnEdge>(0);
        }

        public Path(Path p, DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, DeBruijnEdge edge) {
            lastVertex = graph.getEdgeTarget(edge);
            edges = new ArrayList<DeBruijnEdge>(p.edges);
            edges.add(edge);
            totalScore = p.totalScore + edge.getMultiplicity();
            lowestEdge = ( p.lowestEdge == -1 ) ? edge.getMultiplicity() : Math.min(p.lowestEdge, edge.getMultiplicity());
        }

        public boolean containsEdge(DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, DeBruijnEdge edge) {
            for ( DeBruijnEdge e : edges ) {
                if ( e.equals(graph, edge))
                    return true;
            }

            return false;
        }

        public List<DeBruijnEdge> getEdges() { return edges; }

        public int getScore() { return totalScore; }

        public int getLowestEdge() { return lowestEdge; }

        public DeBruijnVertex getLastVertexInPath() { return lastVertex; }

        // assumes uncleaned edges and vertices, so each edge just adds one more base to the string
        public byte[] getBases(final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph) {
            if(edges.size() == 0) { return lastVertex.printableSequence; }

            int length = 0;
            length += graph.getEdgeSource( edges.get(0) ).printableSequence.length;
            for ( final DeBruijnEdge e : edges ) {
                length += 1;
            }

            byte[] bases = new byte[length];
            int curPos = 0;
            for( final byte b : graph.getEdgeSource( edges.get(0) ).printableSequence ) {
                bases[curPos++] = b;
            }
            for ( final DeBruijnEdge e : edges ) {
                bases[curPos++] = graph.getEdgeTarget( e ).printableSequence[graph.getEdgeTarget( e ).printableSequence.length-1];
            }
            return bases;
        }

        // assumes uncleaned edges and vertices, so each edge just adds one more base to the string
        public byte[] getBases(final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final int length) {

            byte[] bases = new byte[length];
            int curPos = 0;

            if( edges.size() == 0 ) {
                for( final byte b : lastVertex.printableSequence ) {
                    bases[curPos++] = b;
                }

            } else {
                for( final byte b : graph.getEdgeSource( edges.get(0) ).printableSequence ) {
                    bases[curPos++] = b;
                }
                for ( final DeBruijnEdge e : edges ) {
                    bases[curPos++] = graph.getEdgeTarget( e ).printableSequence[graph.getEdgeTarget( e ).printableSequence.length-1];
                }
            }
            if(length - curPos > 40) { return null; } // BUGBUG: not sure why there are several very short paths
            for( int iii = curPos; iii < length; iii++ ) {
                bases[iii] = (byte) 'N';
            }
            return bases;
        }
    }

    protected static class PathComparator implements Comparator<Path> {
        public int compare(final Path path1, final Path path2) {
            return path1.totalScore - path2.totalScore;
        }
    }

    public static List<Path> getKBestPaths(DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, int k) {
        PriorityQueue<Path> bestPaths = new PriorityQueue<Path>(k, new PathComparator());

        // run a DFS for best paths
        for ( final DeBruijnVertex v : graph.vertexSet() ) {
            if ( graph.inDegreeOf(v) == 0 ) {
                findBestPaths(graph, new Path(v), k, bestPaths);
            }
        }

        return new ArrayList<Path>(bestPaths);
    }

    private static void findBestPaths(DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, Path path, int k, PriorityQueue<Path> bestPaths) {
        findBestPaths(graph, path, k, bestPaths, new MyInt());
    }

    private static void findBestPaths(DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, Path path, int k, PriorityQueue<Path> bestPaths, MyInt n) {

        // did we hit the end of a path?
        if ( allOutgoingEdgesHaveBeenVisited(graph, path) ) {
            if ( bestPaths.size() < k ) {
                bestPaths.add(path);
            } else if ( bestPaths.peek().totalScore < path.totalScore ) {
                bestPaths.remove();
                bestPaths.add(path);
            }

        } else if( n.val > 500000) {
            // do nothing, just return
        } else {
            // recursively run DFS
            final ArrayList<DeBruijnEdge> edgeArrayList = new ArrayList<DeBruijnEdge>();
            edgeArrayList.addAll(graph.outgoingEdgesOf(path.lastVertex));
            Collections.sort(edgeArrayList);
            Collections.reverse(edgeArrayList);
            for ( final DeBruijnEdge edge : edgeArrayList ) {
                // make sure the edge is not already in the path
                if ( path.containsEdge(graph, edge) )
                    continue;

                final Path newPath = new Path(path, graph, edge);
                n.val++;
                findBestPaths(graph, newPath, k, bestPaths, n);

            }
        }
    }

    private static boolean allOutgoingEdgesHaveBeenVisited(DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, Path path) {
        for ( final DeBruijnEdge edge : graph.outgoingEdgesOf(path.lastVertex) ) {
            if ( !path.containsEdge(graph, edge) ) {
                return false;
            }
        }
        return true;
    }
}
