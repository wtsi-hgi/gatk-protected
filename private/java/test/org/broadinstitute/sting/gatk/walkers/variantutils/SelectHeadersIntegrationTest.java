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
    public static String baseTestString(String args) {
        return "-T SelectHeaders -R " + hg19Reference + " -L 1 -o %s --no_cmdline_in_header" + args;
    }

    @Test
    public void testSelectHeaderName() {
        String testfile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -hn FILTER -hn INFO -hn FORMAT --variant " + testfile),
                1,
                Arrays.asList("de03c3c170398c5657ec5b9cbb56fc9b")
        );

        executeTest("testSelectHeaderName--" + testfile, spec);
    }

    @Test
    public void testSelectHeaderExpression() {
        String testfile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -he '(FILTER|INFO|FORMAT)' --variant " + testfile),
                1,
                Arrays.asList("de03c3c170398c5657ec5b9cbb56fc9b")
        );

        executeTest("testSelectHeaderExpression--" + testfile, spec);
    }

    @Test
    public void testExcludeHeaderName() {
        String testfile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -xl_hn CombineVariants --variant " + testfile),
                1,
                Arrays.asList("851f8a46f64930e644ffd364f134870f")
        );

        executeTest("testExcludeHeaderName--" + testfile, spec);
    }

    @Test
    public void testIncludeReferenceName() {
        String testfile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -irn --variant " + testfile),
                1,
                Arrays.asList("91de3d2db9cac963073150192d8b289e")
        );

        executeTest("testIncludeReferenceName--" + testfile, spec);
    }

    @Test
    public void testIncludeIntervals() {
        String testfile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -iln --variant " + testfile),
                1,
                Arrays.asList("26c15aac907b6d62ef9edb85eeb89b9f")
        );

        executeTest("testIncludeIntervals--" + testfile, spec);
    }


    @Test
    public void testComplexSelection() {
        String testfile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -he '(FILTER|INFO|FORMAT)' -irn -iln --variant " + testfile),
                1,
                Arrays.asList("6e243821aa868c9a56f938419aa698f2")
        );

        executeTest("testComplexSelection--" + testfile, spec);
    }

    @Test
    public void testParallelization() {
        String testfile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec;

        spec = new WalkerTestSpec(
                baseTestString(" -hn FILTER -hn INFO -hn FORMAT --variant " + testfile + " -nt 2"),
                1,
                Arrays.asList("de03c3c170398c5657ec5b9cbb56fc9b")
        );
        executeTest("testParallelization (2 threads)--" + testfile, spec);

        spec = new WalkerTestSpec(
                baseTestString(" -hn FILTER -hn INFO -hn FORMAT --variant " + testfile + " -nt 4"),
                1,
                Arrays.asList("de03c3c170398c5657ec5b9cbb56fc9b")
        );

        executeTest("testParallelization (4 threads)--" + testfile, spec);
    }
}
