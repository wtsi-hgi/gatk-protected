/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SelectHeadersIntegrationTest extends WalkerTest {
    private static String testfile = privateTestDir + "NA12878.hg19.example1.vcf";

    public static String baseTestString(String args) {
        return "-T SelectHeaders -R " + hg19Reference + " -L 1 -o %s --no_cmdline_in_header" + args;
    }

    @Test
    public void testSelectHeaderName() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -hn FILTER -hn INFO -hn FORMAT --variant " + testfile),
                1,
                Arrays.asList("bb5cbf5c82b3a8d2630ff9055361903d")
        );

        executeTest("testSelectHeaderName--" + testfile, spec);
    }

    @Test
    public void testSelectHeaderExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -he '(FILTER|INFO|FORMAT)' --variant " + testfile),
                1,
                Arrays.asList("bb5cbf5c82b3a8d2630ff9055361903d")
        );

        executeTest("testSelectHeaderExpression--" + testfile, spec);
    }

    @Test
    public void testExcludeHeaderName() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -xl_hn CombineVariants --variant " + testfile),
                1,
                Arrays.asList("009d0ded389c895a3d9f47513fc0e842")
        );

        executeTest("testExcludeHeaderName--" + testfile, spec);
    }

    @Test
    public void testIncludeReferenceName() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -irn --variant " + testfile),
                1,
                Arrays.asList("1f070a462362c1fc66f88d9438b53533")
        );

        executeTest("testIncludeReferenceName--" + testfile, spec);
    }

    @Test
    public void testIncludeIntervals() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -iln --variant " + testfile),
                1,
                Arrays.asList("864168eb5c3c4ac27fe5b62b2f5cdb4a")
        );

        executeTest("testIncludeIntervals--" + testfile, spec);
    }


    @Test
    public void testComplexSelection() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -he '(FILTER|INFO|FORMAT)' -irn -iln --variant " + testfile),
                1,
                Arrays.asList("91158df3e27ec9dbd3b84fe3e8342a21")
        );

        executeTest("testComplexSelection--" + testfile, spec);
    }

    @Test
    public void testParallelization2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -hn FILTER -hn INFO -hn FORMAT --variant " + testfile + " -nt 2"),
                1,
                Arrays.asList("bb5cbf5c82b3a8d2630ff9055361903d")
        );
        executeTest("testParallelization (2 threads)--" + testfile, spec);
    }

    @Test
    public void testParallelization4() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -hn FILTER -hn INFO -hn FORMAT --variant " + testfile + " -nt 4"),
                1,
                Arrays.asList("bb5cbf5c82b3a8d2630ff9055361903d")
        );

        executeTest("testParallelization (4 threads)--" + testfile, spec);
    }
}
