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
                Arrays.asList("8b921c155d03bb9ba4d671c17729fea8")
        );

        executeTest("testSelectHeaderName--" + testfile, spec);
    }

    @Test
    public void testSelectHeaderExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -he '(FILTER|INFO|FORMAT)' --variant " + testfile),
                1,
                Arrays.asList("8b921c155d03bb9ba4d671c17729fea8")
        );

        executeTest("testSelectHeaderExpression--" + testfile, spec);
    }

    @Test
    public void testExcludeHeaderName() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -xl_hn CombineVariants --variant " + testfile),
                1,
                Arrays.asList("ee19a3859dce7a5cafe813c62f9f2efc")
        );

        executeTest("testExcludeHeaderName--" + testfile, spec);
    }

    @Test
    public void testIncludeReferenceName() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -irn --variant " + testfile),
                1,
                Arrays.asList("c48896abe8898577e9ac52556f8778b9")
        );

        executeTest("testIncludeReferenceName--" + testfile, spec);
    }

    @Test
    public void testIncludeIntervals() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -iln --variant " + testfile),
                1,
                Arrays.asList("624909676e654d33d953b0bac5b74bb6")
        );

        executeTest("testIncludeIntervals--" + testfile, spec);
    }


    @Test
    public void testComplexSelection() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -he '(FILTER|INFO|FORMAT)' -irn -iln --variant " + testfile),
                1,
                Arrays.asList("ec94f04c294cf4779fb892a460f9412e")
        );

        executeTest("testComplexSelection--" + testfile, spec);
    }

    @Test
    public void testParallelization2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -hn FILTER -hn INFO -hn FORMAT --variant " + testfile + " -nt 2"),
                1,
                Arrays.asList("8b921c155d03bb9ba4d671c17729fea8")
        );
        executeTest("testParallelization (2 threads)--" + testfile, spec);
    }

    @Test
    public void testParallelization4() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -hn FILTER -hn INFO -hn FORMAT --variant " + testfile + " -nt 4"),
                1,
                Arrays.asList("8b921c155d03bb9ba4d671c17729fea8")
        );

        executeTest("testParallelization (4 threads)--" + testfile, spec);
    }
}
