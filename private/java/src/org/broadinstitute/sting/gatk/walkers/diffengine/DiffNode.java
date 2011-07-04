/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.diffengine;

import org.broadinstitute.sting.utils.Utils;
import scala.Left;
import scala.xml.Null;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/4/11
 * Time: 12:55 PM
 *
 * An interface that must be implemented to allow us to calculate differences
 * between structured objects
 */
public class DiffNode extends DiffElement {
    public final static DiffNode ROOT = new DiffNode("ROOT", null, null, null);

    private List<DiffLeaf> leaves;
    private List<DiffNode> nodes;

    private static List<DiffLeaf> emptyLeaves() { return new ArrayList<DiffLeaf>(); }
    private static List<DiffNode> emptyNodes() { return new ArrayList<DiffNode>(); }

    private DiffNode(String name, DiffNode parent, List<DiffLeaf> leaves, List<DiffNode> nodes) {
        super(name, parent);
        this.leaves = leaves;
        this.nodes = nodes;
    }

    public static DiffNode withLeavesAndNodes(String name, DiffNode parent, List<DiffLeaf> leaves, List<DiffNode> nodes) {
        return new DiffNode(name, parent, leaves, nodes);
    }

    public static DiffNode rooted(String name) {
        return empty(name, ROOT);
    }

    public static DiffNode empty(String name, DiffNode parent) {
        return withLeavesAndNodes(name, parent, emptyLeaves(), emptyNodes());
    }

    public static DiffNode withLeaves(String name, DiffNode parent, List<DiffLeaf> leaves) {
        return withLeavesAndNodes(name, parent, leaves, emptyNodes());
    }

    public static DiffNode withNodes(String name, DiffNode parent, List<DiffNode> nodes) {
        return withLeavesAndNodes(name, parent, emptyLeaves(), nodes);
    }

    public boolean isRoot() { return this == ROOT; }

    public void add(DiffLeaf leaf) {
        leaves.add(leaf);
    }

    public void add(String name, Object value) {
        add(new DiffLeaf(name, this, value));
    }

    public void add(DiffNode node) {
        nodes.add(node);
    }

    public String toString() {
        return toString(0);
    }

    public String toString(int offset) {
        String off = offset > 0 ? Utils.dupString(' ', offset) : "";
        String off2 = Utils.dupString(' ', offset+2);
        StringBuilder b = new StringBuilder();

        b.append(String.format("%s%s%n", off, super.toString()));
        for ( DiffLeaf leaf : leaves )
            b.append(off2).append(leaf.toString()).append('\n');
        for ( DiffNode node : nodes )
            b.append(node.toString(offset + 4));

        return b.toString();
    }

    public List<DiffLeaf> getLeaves() {
        return leaves;
    }

    public List<DiffNode> getNodes() {
        return nodes;
    }

    public DiffLeaf getLeaf(String name) {
        for ( DiffLeaf leaf : getLeaves() )
            if ( leaf.getName().equals(name) )
                return leaf;
        return null;
    }

    public DiffNode getNode(String name) {
        for ( DiffNode node : getNodes() )
            if ( node.getName().equals(name) )
                return node;
        return null;
    }
}
