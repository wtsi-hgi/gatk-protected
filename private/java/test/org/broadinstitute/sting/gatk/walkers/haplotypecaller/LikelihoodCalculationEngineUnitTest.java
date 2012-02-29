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

// our package
package org.broadinstitute.sting.gatk.walkers.haplotypecaller;


// the imports for unit testing.


import net.sf.picard.util.QualityUtil;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


public class LikelihoodCalculationEngineUnitTest extends BaseTest {
    final static boolean EXTENSIVE_TESTING = false;
    LikelihoodCalculationEngine engine;

    @BeforeSuite
    public void before() {
        engine = new LikelihoodCalculationEngine(0, 0, false);
    }

    // --------------------------------------------------------------------------------
    //
    // Provider
    //
    // --------------------------------------------------------------------------------

    private class BasicLikelihoodTestProvider extends TestDataProvider {
        final String ref, read;
        final byte[] refBasesWithContext, readBasesWithContext;
        final int baseQual, insQual, delQual, gcp;
        final int expectedQual;
        final static String CONTEXT = "ACGTTGCA";
        final static int DEFAULT_GCP = 10;

        public BasicLikelihoodTestProvider(final String ref, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual) {
            this(ref, read, baseQual, insQual, delQual, expectedQual, DEFAULT_GCP);
        }

        public BasicLikelihoodTestProvider(final String ref, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual, final int gcp) {
            super(BasicLikelihoodTestProvider.class, String.format("ref=%s read=%s b/i/d/c quals = %d/%d/%d/%d e[qual]=%d", ref, read, baseQual, insQual, delQual, gcp, expectedQual));
            this.baseQual = baseQual;
            this.delQual = delQual;
            this.insQual = insQual;
            this.gcp = gcp;
            this.read = read;
            this.ref = ref;
            this.expectedQual = expectedQual;

            refBasesWithContext = asBytes(ref);
            readBasesWithContext = asBytes(read);
        }

        public double expectedLogL() {
            return expectedQual / -10;
        }

        public double tolerance() {
            return 0.1; // TODO FIXME arbitrary
        }

        public double calcLogL() {
            final int X_METRIC_LENGTH = readBasesWithContext.length+1;
            final int Y_METRIC_LENGTH = refBasesWithContext.length+1;

            double[][] matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
            double[][] XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
            double[][] YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];

            double logL = engine.computeReadLikelihoodGivenHaplotypeAffineGaps(
                    refBasesWithContext, readBasesWithContext,
                    qualAsBytes(baseQual), qualAsDouble(insQual), qualAsDouble(delQual),
                    qualAsDouble(gcp), 0, matchMetricArray, XMetricArray, YMetricArray);

            return logL;
        }

        private final byte[] asBytes(final String bases) {
            return (CONTEXT + bases + CONTEXT).getBytes();
        }

        private byte[] qualAsBytes(final int phredQual) {
            final byte phredQuals[] = new byte[readBasesWithContext.length];

            // initialize everything to MASSIVE_QUAL so it cannot be moved by HMM
            Arrays.fill(phredQuals, Byte.MAX_VALUE);

            // update just the bases corresponding to the provided micro read with the quality scores
            for ( int i = 0; i < read.length(); i++)
                phredQuals[i + CONTEXT.length()] = (byte)phredQual;

            return phredQuals;
        }

        private double[] qualAsDouble(final int phredQual) {
            final double phredQuals[] = new double[readBasesWithContext.length];

            // initialize everything to MASSIVE_QUAL so it cannot be moved by HMM
            Arrays.fill(phredQuals, -1000.0);

            // update just the bases corresponding to the provided micro read with the quality scores
            for ( int i = 0; i < read.length(); i++)
                phredQuals[i + CONTEXT.length()] = phredQual / -10.0;

            return phredQuals;
        }
    }

    @DataProvider(name = "BasicLikelihoodTestProvider")
    public Object[][] makeBasicLikelihoodTests() {
        // context on either side is ACGTTGCA REF ACGTTGCA
        {
            final int baseQual = 30;
            final int insQual = 40;
            final int delQual = 50;
            new BasicLikelihoodTestProvider("C", "A", baseQual, insQual, delQual, baseQual);
            new BasicLikelihoodTestProvider("C", "C", baseQual, insQual, delQual, 0);
            new BasicLikelihoodTestProvider("C", "G", baseQual, insQual, delQual, baseQual);
            new BasicLikelihoodTestProvider("C", "T", baseQual, insQual, delQual, baseQual);
            new BasicLikelihoodTestProvider("C", "", baseQual, insQual, delQual, delQual);
            new BasicLikelihoodTestProvider("C", "CC", baseQual, insQual, delQual, insQual);
        }

        // test all combinations
        final List<Integer> baseQuals = EXTENSIVE_TESTING ? Arrays.asList(10, 20, 30, 40, 50) : Arrays.asList(30);
        final List<Integer> indelQuals = EXTENSIVE_TESTING ? Arrays.asList(10, 20, 30, 40, 50) : Arrays.asList(30);
        final List<Integer> gcps = EXTENSIVE_TESTING ? Arrays.asList(10, 20, 30) : Arrays.asList(10);
        final List<Integer> sizes = EXTENSIVE_TESTING ? Arrays.asList(2,3,4,5,6,7,8,9,10) : Arrays.asList(2, 5);

        for ( final int baseQual : baseQuals ) {
            for ( final int indelQual : indelQuals ) {
                for ( final int gcp : gcps ) {

                    // test substitutions
                    for ( final byte refBase : BaseUtils.BASES ) {
                        for ( final byte readBase : BaseUtils.BASES ) {
                            final String ref  = new String(new byte[]{refBase});
                            final String read = new String(new byte[]{readBase});
                            final int expected = refBase == readBase ? 0 : baseQual;
                            new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp);
                        }
                    }

                    // test insertions and deletions
                    for ( final int size : sizes ) {
                        for ( final byte base : BaseUtils.BASES ) {
                            final int expected = indelQual + (size - 2) * gcp;

                            for ( boolean insertionP : Arrays.asList(true, false)) {
                                final String small = Utils.dupString((char)base, 1);
                                final String big = Utils.dupString((char)base, size);

                                final String ref = insertionP ? small : big;
                                final String read = insertionP ? big : small;

                                new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp);
                            }
                        }
                    }
                }
            }
        }

        return BasicLikelihoodTestProvider.getTests(BasicLikelihoodTestProvider.class);
    }

    @Test(dataProvider = "BasicLikelihoodTestProvider", enabled = false)
    public void testBasicLikelihoods(BasicLikelihoodTestProvider cfg) {
        double calculatedLogL = cfg.calcLogL();
        double expectedLogL = cfg.expectedLogL();
        logger.warn(String.format("Test: logL calc=%.2f expected=%.2f for %s", calculatedLogL, expectedLogL, cfg.toString()));
        Assert.assertEquals(calculatedLogL, expectedLogL, cfg.tolerance());
    }
}