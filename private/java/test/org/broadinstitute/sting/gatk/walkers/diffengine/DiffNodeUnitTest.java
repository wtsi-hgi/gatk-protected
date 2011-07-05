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


import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

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
    static DiffValue Value_A = new DiffValue("A", MY_ROOT, "A");
    static DiffValue Value_B = new DiffValue("B", MY_ROOT, "B");
    static DiffNode NODE_C = DiffNode.empty("C", MY_ROOT);
    static DiffNode NODE_D = DiffNode.empty("D", MY_ROOT);
    static DiffValue Value_E = new DiffValue("E", NODE_C, "E");
    static DiffValue Value_F = new DiffValue("F", NODE_D, "F");
    static DiffValue Value_G = new DiffValue("G", NODE_D, "G");

    static {
        MY_ROOT.add(Value_A);
        MY_ROOT.add(Value_B);
        MY_ROOT.add(NODE_C);
        MY_ROOT.add(NODE_D);
        NODE_C.add(Value_E);
        NODE_D.add(Value_F);
        NODE_D.add(Value_G);
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

        private ElementTest(DiffValue elt, DiffValue parent, String name, String fullName) {
            this(elt.getBinding(), parent.getBinding(), name, fullName);
        }

        private ElementTest(DiffElement elt, DiffElement parent, String name, String fullName) {
            this.elt = elt;
            this.name = name;
            this.fullName = fullName;
            this.parent = parent;
        }

        public String toString() {
            return String.format("ElementTest elt=%s name=%s fullName=%s parent=%s",
                    elt.toOneLineString(), name, fullName, parent.getName());
        }
    }

    @DataProvider(name = "elementdata")
    public Object[][] createElementData() {
        List<ElementTest> params = new ArrayList<ElementTest>();

        params.add(new ElementTest(MY_ROOT.getBinding(), DiffElement.ROOT, "MY_ROOT", "MY_ROOT"));
        params.add(new ElementTest(NODE_C, MY_ROOT, "C", "MY_ROOT.C"));
        params.add(new ElementTest(NODE_D, MY_ROOT, "D", "MY_ROOT.D"));
        params.add(new ElementTest(Value_A, MY_ROOT, "A", "MY_ROOT.A"));
        params.add(new ElementTest(Value_B, MY_ROOT, "B", "MY_ROOT.B"));
        params.add(new ElementTest(Value_E, NODE_C, "E", "MY_ROOT.C.E"));
        params.add(new ElementTest(Value_F, NODE_D, "F", "MY_ROOT.D.F"));
        params.add(new ElementTest(Value_G, NODE_D, "G", "MY_ROOT.D.G"));

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
    // DiffValue testing routines
    //
    // --------------------------------------------------------------------------------

    private class LeafTest {
        public DiffValue diffvalue;
        public Object value;

        private LeafTest(DiffValue diffvalue, Object value) {
            this.diffvalue = diffvalue;
            this.value = value;
        }

        public String toString() {
            return String.format("LeafTest diffvalue=%s value=%s", diffvalue.toOneLineString(), value);
        }
    }

    @DataProvider(name = "leafdata")
    public Object[][] createLeafData() {
        List<LeafTest> params = new ArrayList<LeafTest>();

        params.add(new LeafTest(Value_A, "A"));
        params.add(new LeafTest(Value_B, "B"));
        params.add(new LeafTest(Value_E, "E"));
        params.add(new LeafTest(Value_F, "F"));
        params.add(new LeafTest(Value_G, "G"));

        List<Object[]> params2 = new ArrayList<Object[]>();
        for ( LeafTest x : params ) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "leafdata")
    public void testLeafMethods(LeafTest test) {
        Assert.assertNotNull(test.diffvalue.getValue());
        Assert.assertEquals(test.diffvalue.getValue(), test.value);
    }

    // --------------------------------------------------------------------------------
    //
    // Node testing routines
    //
    // --------------------------------------------------------------------------------

    private class NodeTest {
        public DiffNode node;
        public Set<String> fields;
        public Set<String> subnodes;
        public Set<String> allNames;

        private NodeTest(DiffNode node, List<String> fields, List<String> subnodes) {
            this.node = node;
            this.fields = new HashSet<String>(fields);
            this.subnodes = new HashSet<String>(subnodes);
            this.allNames = new HashSet<String>(fields);
            allNames.addAll(subnodes);
        }

        public String toString() {
            return String.format("NodeTest node=%s fields=%s subnodes=%s",
                    node.toOneLineString(), fields, subnodes);
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
        Assert.assertNotNull(test.node.getElements());

        for ( String name : test.allNames ) {
            DiffElement elt = test.node.getElement(name);
            Assert.assertNotNull(elt, "Failed to find field " + elt + " in " + test.node);
            Assert.assertEquals(elt.getName(), name);
            Assert.assertEquals(elt.getValue().isAtomic(), test.fields.contains(name), "Failed atomic/compound expectation: " + test.node);
        }
    }

    // NOTE: add routines are being implicitly tested by the creation of the data structures

    @Test(enabled = true, dataProvider = "nodedata")
    public void testCounts(NodeTest test) {
        Assert.assertEquals(test.node.getElements().size(), test.allNames.size());
        Assert.assertEquals(test.node.getElementNames(), test.allNames);
    }

    // --------------------------------------------------------------------------------
    //
    // fromString testing routines
    //
    // --------------------------------------------------------------------------------

    private class FromStringTest {
        public String string;
        public DiffElement expected;

        private FromStringTest(String string, DiffElement expected) {
            this.string = string;
            this.expected = expected;
        }

        public String toString() {
            return String.format("FromStringTest string=%s expected=%s", string, expected.toOneLineString());
        }
    }

    @DataProvider(name = "fromstringdata")
    public Object[][] createFromData() {
        List<FromStringTest> params = new ArrayList<FromStringTest>();

        params.add(new FromStringTest("A=A", Value_A.getBinding()));
        params.add(new FromStringTest("B=B", Value_B.getBinding()));
        params.add(new FromStringTest("C=(E=E)", NODE_C.getBinding()));
        params.add(new FromStringTest("D=(F=F G=G)", NODE_D.getBinding()));

        List<Object[]> params2 = new ArrayList<Object[]>();
        for ( FromStringTest x : params ) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "fromstringdata")
    public void parseFromString(FromStringTest test) {
        logger.warn("Testing from string: " + test.string);
        DiffElement elt = DiffNode.fromString(test.string);
        Assert.assertEquals(elt.toOneLineString(), test.expected.toOneLineString());
    }
}