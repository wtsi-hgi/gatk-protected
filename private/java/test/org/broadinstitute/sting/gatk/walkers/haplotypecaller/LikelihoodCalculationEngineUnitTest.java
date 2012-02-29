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


import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


public class LikelihoodCalculationEngineUnitTest extends BaseTest {
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
        final int baseQual, insQual, delQual;
        final int gcp = 10;
        final int expectedQual;
        final static String CONTEXT = "ACGTACGT";

        public BasicLikelihoodTestProvider(final String ref, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual) {
            super(BasicLikelihoodTestProvider.class, String.format("ref=%s read=%s b/i/d quals = %d/%d/%d e[qual]=%d", ref, read, baseQual, insQual, delQual, expectedQual));
            this.baseQual = baseQual;
            this.delQual = delQual;
            this.insQual = insQual;
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
            Arrays.fill(phredQuals, (byte)phredQual);
            return phredQuals;
        }

        private double[] qualAsDouble(final int phredQual) {
            final double log10ErrorRates[] = new double[readBasesWithContext.length];
            final double log10ErrorRate = phredQual / -10.0;
            Arrays.fill(log10ErrorRates, log10ErrorRate);
            return log10ErrorRates;
        }
    }

    @DataProvider(name = "BasicLikelihoodTestProvider")
    public Object[][] makeBasicLikelihoodTests() {
        // context on either side is ACGTACGT REF ACGTACGT
        int baseQual = 30;
        int insQual = 40;
        int delQual = 50;
        new BasicLikelihoodTestProvider("C", "A", baseQual, insQual, delQual, baseQual);
        new BasicLikelihoodTestProvider("C", "C", baseQual, insQual, delQual, 0);
        new BasicLikelihoodTestProvider("C", "G", baseQual, insQual, delQual, baseQual);
        new BasicLikelihoodTestProvider("C", "T", baseQual, insQual, delQual, baseQual);
        new BasicLikelihoodTestProvider("C", "", baseQual, insQual, delQual, delQual);
        new BasicLikelihoodTestProvider("C", "CC", baseQual, insQual, delQual, insQual);

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