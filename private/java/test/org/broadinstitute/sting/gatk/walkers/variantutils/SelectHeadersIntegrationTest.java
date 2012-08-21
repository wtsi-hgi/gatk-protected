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
                Arrays.asList("c338b844c3cd1fe8fbe3df4d7b07321e")
        );

        executeTest("testSelectHeaderName--" + testfile, spec);
    }

    @Test
    public void testSelectHeaderExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -he '(FILTER|INFO|FORMAT)' --variant " + testfile),
                1,
                Arrays.asList("c338b844c3cd1fe8fbe3df4d7b07321e")
        );

        executeTest("testSelectHeaderExpression--" + testfile, spec);
    }

    @Test
    public void testExcludeHeaderName() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -xl_hn CombineVariants --variant " + testfile),
                1,
                Arrays.asList("7b8e82c246307bdb2460f3c1e776cfdc")
        );

        executeTest("testExcludeHeaderName--" + testfile, spec);
    }

    @Test
    public void testIncludeIntervals() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -iln --variant " + testfile),
                1,
                Arrays.asList("582e1617b76a99b394eafed30ed02eee")
        );

        executeTest("testIncludeIntervals--" + testfile, spec);
    }


    @Test
    public void testComplexSelection() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -he '(FILTER|INFO|FORMAT)' -iln --variant " + testfile),
                1,
                Arrays.asList("1b2eef1f362d2445428f49476556fcaa")
        );

        executeTest("testComplexSelection--" + testfile, spec);
    }

    @Test
    public void testParallelization2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -hn FILTER -hn INFO -hn FORMAT --variant " + testfile + " -nt 2"),
                1,
                Arrays.asList("c338b844c3cd1fe8fbe3df4d7b07321e")
        );
        executeTest("testParallelization (2 threads)--" + testfile, spec);
    }

    @Test
    public void testParallelization4() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -hn FILTER -hn INFO -hn FORMAT --variant " + testfile + " -nt 4"),
                1,
                Arrays.asList("c338b844c3cd1fe8fbe3df4d7b07321e")
        );

        executeTest("testParallelization (4 threads)--" + testfile, spec);
    }
}
