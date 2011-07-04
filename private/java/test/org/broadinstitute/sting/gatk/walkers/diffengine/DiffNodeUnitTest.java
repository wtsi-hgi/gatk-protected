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

// our package
package org.broadinstitute.sting.gatk.walkers.diffengine;


// the imports for unit testing.


import org.apache.xmlbeans.impl.tool.Diff;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * Basic unit test for DifferableReaders in reduced reads
 */
public class DiffNodeUnitTest extends BaseTest {
    // Data is:
    // MY_ROOT
    //   fields: A=A, B=B
    //   nodes: C, D
    //   C: fields: E=E, nodes: none
    //   D: fields: F=F, G=G, nodes: none
    static DiffNode MY_ROOT = DiffNode.rooted("MY_ROOT");
    static DiffLeaf LEAF_A = new DiffLeaf("A", MY_ROOT, "A");
    static DiffLeaf LEAF_B = new DiffLeaf("B", MY_ROOT, "B");
    static DiffNode NODE_C = DiffNode.empty("C", MY_ROOT);
    static DiffNode NODE_D = DiffNode.empty("D", MY_ROOT);
    static DiffLeaf LEAF_E = new DiffLeaf("E", NODE_C, "E");
    static DiffLeaf LEAF_F = new DiffLeaf("F", NODE_D, "F");
    static DiffLeaf LEAF_G = new DiffLeaf("G", NODE_D, "G");

    static {
        MY_ROOT.add(LEAF_A);
        MY_ROOT.add(LEAF_B);
        MY_ROOT.add(NODE_C);
        MY_ROOT.add(NODE_D);
        NODE_C.add(LEAF_E);
        NODE_D.add(LEAF_F);
        NODE_D.add(LEAF_G);
    }


    // --------------------------------------------------------------------------------
    //
    // Element testing routines
    //
    // --------------------------------------------------------------------------------

    private class ElementTest {
        public DiffElement elt;
        public String name;
        public String fullName;
        public DiffElement parent;

        private ElementTest(DiffElement elt, DiffElement parent, String name, String fullName) {
            this.elt = elt;
            this.name = name;
            this.fullName = fullName;
            this.parent = parent;
        }
    }

    @DataProvider(name = "elementdata")
    public Object[][] createElementData() {
        List<ElementTest> params = new ArrayList<ElementTest>();

        params.add(new ElementTest(MY_ROOT, DiffNode.ROOT, "MY_ROOT", "MY_ROOT"));
        params.add(new ElementTest(NODE_C, MY_ROOT, "C", "MY_ROOT.C"));
        params.add(new ElementTest(NODE_D, MY_ROOT, "D", "MY_ROOT.D"));
        params.add(new ElementTest(LEAF_A, MY_ROOT, "A", "MY_ROOT.A"));
        params.add(new ElementTest(LEAF_B, MY_ROOT, "B", "MY_ROOT.B"));
        params.add(new ElementTest(LEAF_E, NODE_C, "E", "MY_ROOT.C.E"));
        params.add(new ElementTest(LEAF_F, NODE_D, "F", "MY_ROOT.D.F"));
        params.add(new ElementTest(LEAF_G, NODE_D, "G", "MY_ROOT.D.G"));

        List<Object[]> params2 = new ArrayList<Object[]>();
        for ( ElementTest x : params ) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "elementdata")
    public void testElementMethods(ElementTest test) {
        Assert.assertNotNull(test.elt.getName());
        Assert.assertNotNull(test.elt.getParent());
        Assert.assertEquals(test.elt.getName(), test.name);
        Assert.assertEquals(test.elt.getParent(), test.parent);
        Assert.assertEquals(test.elt.fullyQualifiedName(), test.fullName);
    }

    // --------------------------------------------------------------------------------
    //
    // leaf testing routines
    //
    // --------------------------------------------------------------------------------

    private class LeafTest {
        public DiffLeaf leaf;
        public Object value;

        private LeafTest(DiffLeaf leaf, Object value) {
            this.leaf = leaf;
            this.value = value;
        }
    }

    @DataProvider(name = "leafdata")
    public Object[][] createLeafData() {
        List<LeafTest> params = new ArrayList<LeafTest>();

        params.add(new LeafTest(LEAF_A, "A"));
        params.add(new LeafTest(LEAF_B, "B"));
        params.add(new LeafTest(LEAF_E, "E"));
        params.add(new LeafTest(LEAF_F, "F"));
        params.add(new LeafTest(LEAF_G, "G"));

        List<Object[]> params2 = new ArrayList<Object[]>();
        for ( LeafTest x : params ) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "leafdata")
    public void testLeafMethods(LeafTest test) {
        Assert.assertNotNull(test.leaf.getValue());
        Assert.assertEquals(test.leaf.getValue(), test.value);
    }

    // --------------------------------------------------------------------------------
    //
    // Node testing routines
    //
    // --------------------------------------------------------------------------------

    private class NodeTest {
        public DiffNode node;
        public List<String> fields;
        public List<String> nodes;

        private NodeTest(DiffNode node, List<String> fields, List<String> nodes) {
            this.node = node;
            this.fields = fields;
            this.nodes = nodes;
        }
    }

    @DataProvider(name = "nodedata")
    public Object[][] createData1() {
        List<NodeTest> params = new ArrayList<NodeTest>();

        params.add(new NodeTest(MY_ROOT, Arrays.asList("A", "B"), Arrays.asList("C", "D")));
        params.add(new NodeTest(NODE_C, Arrays.asList("E"), Collections.<String>emptyList()));
        params.add(new NodeTest(NODE_D, Arrays.asList("F", "G"), Collections.<String>emptyList()));

        List<Object[]> params2 = new ArrayList<Object[]>();
        for ( NodeTest x : params ) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "nodedata")
    public void testNodeAccessors(NodeTest test) {
        Assert.assertNotNull(test.node.getLeaves());
        Assert.assertNotNull(test.node.getNodes());

        for ( String field : test.fields ) {
            DiffLeaf fieldNode = test.node.getLeaf(field);
            Assert.assertNotNull(fieldNode, "Failed to find field " + fieldNode + " in " + test.node);
            Assert.assertEquals(fieldNode.getName(), field);
        }

        for ( String node : test.nodes ) {
            DiffNode subNode = test.node.getNode(node);
            Assert.assertNotNull(subNode, "Failed to find node " + node + " in " + test.node);
            Assert.assertEquals(subNode.getName(), node);
        }
    }

    @Test(enabled = true, dataProvider = "nodedata")
    public void testCounts(NodeTest test) {
        Assert.assertEquals(test.node.getLeaves().size(), test.fields.size());
        Assert.assertEquals(test.node.getNodes().size(), test.nodes.size());
    }
}