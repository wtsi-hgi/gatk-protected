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
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
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
    public final static DiffNode ROOT = new DiffNode("ROOT", null, null);
    public boolean isRoot() { return this == ROOT; }

    private Map<String, DiffElement> leaves;

    private static Map<String, DiffElement> emptyLeaves() { return new HashMap<String, DiffElement>(); }

    private DiffNode(String name, DiffNode parent, Map<String, DiffElement> leaves) {
        super(name, parent);
        this.leaves = leaves;
    }

    public static DiffNode rooted(String name) {
        return empty(name, ROOT);
    }

    public static DiffNode empty(String name, DiffNode parent) {
        return withLeaves(name, parent, emptyLeaves());
    }

    public static DiffNode withLeaves(String name, DiffNode parent, List<DiffElement> leaves) {
        DiffNode node = empty(name, parent);
        node.add(leaves);
        return node;
    }

    public static DiffNode withLeaves(String name, DiffNode parent, Map<String, DiffElement> leaves) {
        return new DiffNode(name, parent, leaves);
    }

    public void add(DiffElement elt) {
        leaves.put(elt.getName(), elt);
    }

    public void add(Collection<DiffElement> elts) {
        for ( DiffElement e : elts )
            add(e);
    }

    public void add(String name, Object value) {
        add(new DiffLeaf(name, this, value));
    }


    public String toString() {
        return toString(0);
    }

    public String toString(int offset) {
        String off = offset > 0 ? Utils.dupString(' ', offset) : "";
        String off2 = Utils.dupString(' ', offset+2);
        StringBuilder b = new StringBuilder();

        b.append(String.format("%s%s%n", off, super.toString()));
        for ( DiffLeaf leaf : getLeaves() )
            b.append(off2).append(leaf.toString()).append('\n');
        for ( DiffNode node : getNodes() )
            b.append(node.toString(offset + 4));

        return b.toString();
    }

    public Collection<String> getElementNames() {
        return leaves.keySet();
    }

    public Collection<DiffElement> getElements() {
        return leaves.values();
    }

    private <T> Collection<T> getElements(Class<T> restriction) {
        List<T> l = new ArrayList<T>();
        for ( DiffElement e : getElements() ) {
            if ( restriction.isInstance(e) )
                l.add((T)e);
        }
        return l;
    }

    public Collection<DiffLeaf> getLeaves() {
        return getElements(DiffLeaf.class);
    }

    public Collection<DiffNode> getNodes() {
        return getElements(DiffNode.class);
    }

    public DiffElement getElement(String name) {
        for ( DiffElement elt : getElements() )
            if ( elt.getName().equals(name) )
                return elt;
        return null;
    }

    public <T> T getElement(String name, Class<T> restriction) {
        DiffElement elt = getElement(name);
        if ( elt == null )
            return null;
        else if ( restriction.isInstance(elt) )
            return (T)elt;
        else
            throw new ReviewedStingException("Requested object of type " + restriction + " but found binding to class " + elt.getClass());
    }

    public DiffLeaf getLeaf(String name) {
        return getElement(name, DiffLeaf.class);
    }

    public DiffNode getNode(String name) {
        return getElement(name, DiffNode.class);
    }

    // --------------------------------------------------------------------------------
    //
    // fromString and toOneLineString
    //
    // --------------------------------------------------------------------------------

    public static DiffElement fromString(String tree) {
        return fromString(tree, ROOT);
    }

    /**
     * Doesn't support full tree structure parsing
     * @param tree
     * @param parent
     * @return
     */
    private static DiffElement fromString(String tree, DiffNode parent) {
        // X=(A=A B=B C=(D=D))
        String[] parts = tree.split("=", 2);
        if ( parts.length != 2 )
            throw new ReviewedStingException("Unexpected tree structure: " + tree + " parts=" + parts);
        String name = parts[0];
        String value = parts[1];

        if ( value.length() == 0 )
            throw new ReviewedStingException("Illegal tree structure: " + value + " at " + tree);

        if ( value.charAt(0) == '(' ) {
            if ( ! value.endsWith(")") )
                throw new ReviewedStingException("Illegal tree structure.  Missing ): " + value + " at " + tree);
            String subtree = value.substring(1, value.length()-1);
            DiffNode rec = DiffNode.empty(name, parent);
            String[] subParts = subtree.split(" ");
            for ( String subPart : subParts ) {
                rec.add(fromString(subPart, rec));
            }
            return rec;
        } else {
            return new DiffLeaf(name, parent, value);
        }
    }

    @Override
    public String toOneLineString() {
        StringBuilder b = new StringBuilder();

        b.append(super.toOneLineString());
        b.append('(');
        List<String> parts = new ArrayList<String>();
        for ( DiffElement elt : getElements() )
            parts.add(elt.toOneLineString());
        b.append(Utils.join(" ", parts));
        b.append(')');

        return b.toString();
    }
}
