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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Basic unit test for DifferableReaders in reduced reads
 */
public class DiffEngineUnitTest extends BaseTest {
    // --------------------------------------------------------------------------------
    //
    // Element testing routines
    //
    // --------------------------------------------------------------------------------

    private class DifferenceTest {
        public DiffElement tree1, tree2;
        public List<String> differences;

        private DifferenceTest(String tree1, String tree2, List<String> differences) {
            this.tree1 = DiffNode.fromString(tree1);
            this.tree2 = DiffNode.fromString(tree2);
            this.differences = differences;
        }
    }

    @DataProvider(name = "trees")
    public Object[][] createTrees() {
        List<DifferenceTest> params = new ArrayList<DifferenceTest>();

        params.add(new DifferenceTest("A=X", "A=X",
                Collections.<String>emptyList()));

        params.add(new DifferenceTest("A=X", "A=Y",
                Arrays.asList("A:X!=Y")));

        params.add(new DifferenceTest("A=X", "B=X",
                Arrays.asList("A:X!=MISSING", "B:MISSING!=X")));

        params.add(new DifferenceTest("A=(X=1)", "A=(X=1)",
                Collections.<String>emptyList()));

        params.add(new DifferenceTest("A=(X=1)", "A=(X=2)",
                Arrays.asList("A.X:1!=2")));

        params.add(new DifferenceTest("A=(X=1)", "A=(X=1 Y=2)",
                Arrays.asList("A.Y:MISSING!=2")));

        List<Object[]> params2 = new ArrayList<Object[]>();
        for ( DifferenceTest x : params ) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "trees")
    public void testDiffs(DifferenceTest test) {
        logger.warn("Test tree1: " + test.tree1.toOneLineString());
        logger.warn("Test tree2: " + test.tree2.toOneLineString());
        logger.warn("Test diff : " + test.differences);
    }
}