package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.jgrapht.graph.DefaultDirectedGraph;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 23, 2011
 */

// simple edge class for connecting nodes in the graph
public class DeBruijnEdge implements Comparable<DeBruijnEdge> {

    private int multiplicity;
    private boolean isRef;

    public DeBruijnEdge() {
        multiplicity = 1;
        isRef = false;
    }

    public DeBruijnEdge( final boolean isRef ) {
        multiplicity = (isRef ? 0 : 1);
        this.isRef = isRef;
    }

    public int getMultiplicity() {
        return multiplicity;
    }

    public void setMultiplicity( final int value ) {
        multiplicity = value;
    }

    public boolean getIsRef() {
        return isRef;
    }

    public void setIsRef( final boolean isRef ) {
        this.isRef = isRef;
    }

    public boolean equals(DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, DeBruijnEdge edge) {
        return (graph.getEdgeSource(this) == graph.getEdgeSource(edge)) && (graph.getEdgeTarget(this) == graph.getEdgeTarget(edge));
    }

    public int compareTo(final DeBruijnEdge that) {
        return this.multiplicity - that.multiplicity;
    }

}
