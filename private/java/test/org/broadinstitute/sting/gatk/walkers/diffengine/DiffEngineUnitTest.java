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


import com.sun.xml.internal.ws.api.pipe.Engine;
import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Basic unit test for DifferableReaders in reduced reads
 */
public class DiffEngineUnitTest extends BaseTest {
    DiffEngine engine;

    @BeforeClass(enabled = true)
    public void createDiffEngine() {
        engine = new DiffEngine(10);
    }


    // --------------------------------------------------------------------------------
    //
    // Element testing routines
    //
    // --------------------------------------------------------------------------------

    private class DifferenceTest {
        public DiffElement tree1, tree2;
        public List<String> differences;

        private DifferenceTest(String tree1, String tree2) {
            this(tree1, tree2, Collections.<String>emptyList());
        }

        private DifferenceTest(String tree1, String tree2, String difference) {
            this(tree1, tree2, Arrays.asList(difference));
        }

        private DifferenceTest(String tree1, String tree2, List<String> differences) {
            this.tree1 = DiffNode.fromString(tree1);
            this.tree2 = DiffNode.fromString(tree2);
            this.differences = differences;
        }

        public String toString() {
            return String.format("tree1=%s tree2=%s diff=%s",
                    tree1.toOneLineString(), tree2.toOneLineString(), differences);
        }
    }

    @DataProvider(name = "trees")
    public Object[][] createTrees() {
        List<DifferenceTest> params = new ArrayList<DifferenceTest>();

        params.add(new DifferenceTest("A=X", "A=X"));
        params.add(new DifferenceTest("A=X", "A=Y", "A:X!=Y"));
        params.add(new DifferenceTest("A=X", "B=X", Arrays.asList("A:X!=MISSING", "B:MISSING!=X")));
        params.add(new DifferenceTest("A=(X=1)", "B=(X=1)", Arrays.asList("A:(X=1)!=MISSING", "B:MISSING!=(X=1)")));
        params.add(new DifferenceTest("A=(X=1)", "A=(X=1)"));
        params.add(new DifferenceTest("A=(X=1 Y=2)", "A=(X=1 Y=2)"));
        params.add(new DifferenceTest("A=(X=1 Y=2 B=(Z=3))", "A=(X=1 Y=2 B=(Z=3))"));
        params.add(new DifferenceTest("A=(X=1)", "A=(X=2)", "A.X:1!=2"));
        params.add(new DifferenceTest("A=(X=1 Y=2 B=(Z=3))", "A=(X=1 Y=2 B=(Z=4))", "A.B.Z:3!=4"));
        params.add(new DifferenceTest("A=(X=1)", "A=(X=1 Y=2)", "A.Y:MISSING!=2"));
        params.add(new DifferenceTest("A=(X=1 Y=2 B=(Z=3))", "A=(X=1 Y=2)", "A.B:(Z=3)!=MISSING"));

        List<Object[]> params2 = new ArrayList<Object[]>();
        for ( DifferenceTest x : params ) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "trees")
    public void testDiffs(DifferenceTest test) {
        logger.warn("Test tree1: " + test.tree1.toOneLineString());
        logger.warn("Test tree2: " + test.tree2.toOneLineString());

        List<Difference> diffs = engine.diff(test.tree1, test.tree2);
        logger.warn("Test expected diff : " + test.differences);
        logger.warn("Observed diffs     : " + diffs);
    }

    @Test(enabled = true)
    public void testLongestCommonPostfix() {
        testLongestCommonPostfixHelper("A", "A", 1);
        testLongestCommonPostfixHelper("A", "B", 0);
        testLongestCommonPostfixHelper("A.B", "A.B", 2);
        testLongestCommonPostfixHelper("A.B.C", "A.B.C", 3);
        testLongestCommonPostfixHelper("A.B.C", "X.B.C", 2);
        testLongestCommonPostfixHelper("A.B.C", "X.Y.C", 1);
        testLongestCommonPostfixHelper("A.B.C", "X.Y.Z", 0);
        testLongestCommonPostfixHelper("A.B.C", "A.X.C", 1);
        testLongestCommonPostfixHelper("A.B.C", "A.X.Z", 0);
        testLongestCommonPostfixHelper("A.B.C", "A.B.Z", 0);
    }

    public void testLongestCommonPostfixHelper(String p1, String p2, int expected) {
        String[] parts1 = p1.split("\\.");
        String[] parts2 = p2.split("\\.");
        int obs = DiffEngine.longestCommonPostfix(parts1, parts2);
        Assert.assertEquals(obs, expected, "p1=" + p1 + " p2=" + p2 + " failed");
    }

    @Test(enabled = true, dependsOnMethods = "testLongestCommonPostfix")
    public void testSummarizePath() {
        testSummarizePathHelper("A", "A", "A");
        testSummarizePathHelper("A", "B", "*");
        testSummarizePathHelper("A.B", "A.B", "A.B");
        testSummarizePathHelper("A.B", "X.B", "*.B");
        testSummarizePathHelper("A.B", "X.Y", "*.*");
        testSummarizePathHelper("A.B.C", "A.B.C", "A.B.C");
        testSummarizePathHelper("A.B.C", "X.B.C", "*.B.C");
        testSummarizePathHelper("A.B.C", "X.Y.C", "*.*.C");
        testSummarizePathHelper("A.B.C", "X.Y.Z", "*.*.*");
        testSummarizePathHelper("A.B.C", "A.X.C", "*.*.C");
        testSummarizePathHelper("A.B.C", "A.X.Z", "*.*.*");
        testSummarizePathHelper("A.B.C", "A.B.Z", "*.*.*");
    }

    public void testSummarizePathHelper(String p1, String p2, String expected) {
        String[] parts1 = p1.split("\\.");
        String[] parts2 = p2.split("\\.");
        int obs = DiffEngine.longestCommonPostfix(parts1, parts2);
        String path = DiffEngine.summarizedPath(parts2, obs);
        Assert.assertEquals(path, expected, "p1=" + p1 + " p2=" + p2 + " failed");
    }
}