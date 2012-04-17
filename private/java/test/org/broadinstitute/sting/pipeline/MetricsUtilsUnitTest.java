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

package org.broadinstitute.sting.pipeline;

import net.sf.picard.metrics.MetricBase;
import net.sf.picard.metrics.MultilevelMetrics;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class MetricsUtilsUnitTest {
    @SuppressWarnings("unused")
    class TestMetrics extends MetricBase {
        public int __HIDDEN;
        public double FOUND1;
        public String _FOUND2;
    }

    @SuppressWarnings("unused")
    class TestMultilevelMetrics extends MultilevelMetrics {
        public int __HIDDEN;
        public double FOUND_BASE1;
        public String _FOUND_BASE2;

        @Override
        public String toString() {
            return String.format("{SAMPLE=%s, LIBRARY=%s}", formatName(SAMPLE), formatName(LIBRARY));
        }
    }

    @DataProvider(name = "metricsFields")
    public Object[][] getMetricsFields() {
       return new Object[][] {
               new Object[] { TestMetrics.class, Arrays.asList("FOUND1", "_FOUND2") },
               new Object[] { TestMultilevelMetrics.class, Arrays.asList("FOUND_BASE1", "_FOUND_BASE2", "SAMPLE", "READ_GROUP", "LIBRARY") },
       };
    }

    @Test(dataProvider = "metricsFields")
    public void testGetFieldNames(Class<? extends MetricBase> metricsClass, Collection<String> expectedFields) {
        List<String> actual = MetricsUtils.getFieldNames(metricsClass);
        List<String> expected = new ArrayList<String>(expectedFields);
        Collections.sort(actual);
        Collections.sort(expected);
        Assert.assertEquals(actual, expected, metricsClass.getName());
    }

    @DataProvider(name = "libraryMetrics")
    public Object[][] getLibraryMetrics() {
        TestMultilevelMetrics sampleOnly = new TestMultilevelMetrics();
        sampleOnly.SAMPLE = "sample";

        TestMultilevelMetrics sampleLibrary1 = new TestMultilevelMetrics();
        sampleLibrary1.SAMPLE = "sample";
        sampleLibrary1.LIBRARY = "library1";

        TestMultilevelMetrics sampleLibrary2 = new TestMultilevelMetrics();
        sampleLibrary2.SAMPLE = "sample";
        sampleLibrary2.LIBRARY = "library2";

        TestMultilevelMetrics sampleLibraryNull = new TestMultilevelMetrics();
        sampleLibraryNull.SAMPLE = "sample";
        sampleLibraryNull.LIBRARY = "null"; // Make sure string 'null' does not match value <null>.

        return new Object[][] {
                new Object[] { Arrays.asList(sampleOnly), null, sampleOnly },
                new Object[] { Arrays.asList(sampleOnly), "library1", null },
                new Object[] { Arrays.asList(sampleLibrary1), "library1", sampleLibrary1 },
                new Object[] { Arrays.asList(sampleLibrary1), null, null },
                new Object[] { Arrays.asList(sampleLibraryNull), "null", sampleLibraryNull },
                new Object[] { Arrays.asList(sampleLibraryNull), null, null },
                new Object[] { Arrays.asList(sampleOnly, sampleLibrary1, sampleLibrary2, sampleLibraryNull), null, sampleOnly },
                new Object[] { Arrays.asList(sampleOnly, sampleLibrary1, sampleLibrary2, sampleLibraryNull), "library1", sampleLibrary1 },
                new Object[] { Arrays.asList(sampleOnly, sampleLibrary1, sampleLibrary2, sampleLibraryNull), "library2", sampleLibrary2 },
                new Object[] { Arrays.asList(sampleOnly, sampleLibrary1, sampleLibrary2, sampleLibraryNull), "null", sampleLibraryNull }
        };
    }

    @Test(dataProvider = "libraryMetrics")
    public void testGetLibraryMetric(List<TestMultilevelMetrics> metrics, String libraryName, TestMultilevelMetrics libraryMetric) {
        Assert.assertSame(MetricsUtils.getLibraryMetric(metrics, libraryName), libraryMetric, "getting metric with library set to " + formatName(libraryName));
    }

    @Test(dataProvider = "libraryMetrics")
    public void testFilterMetrics(List<TestMultilevelMetrics> metrics, String libraryName, TestMultilevelMetrics libraryMetric) {
        if (libraryMetric == null) {
            Assert.assertEquals(MetricsUtils.filterMetrics(metrics, "LIBRARY", libraryName), Collections.emptyList(), "filtering metrics with library set to " + formatName(libraryName));
        } else {
            List<TestMultilevelMetrics> actual = MetricsUtils.filterMetrics(metrics, "LIBRARY", libraryName);
            Assert.assertEquals(actual.size(), 1, "filtering metrics count with library set to " + libraryName);
            Assert.assertSame(MetricsUtils.filterMetrics(metrics, "LIBRARY", libraryName).get(0), libraryMetric, "filtering metrics with library set to " + formatName(libraryName));
        }
    }

    @DataProvider(name = "duplicatedMetrics")
    public Object[][] getDuplicatedMetrics() {
        TestMultilevelMetrics sampleOnly = new TestMultilevelMetrics();
        sampleOnly.SAMPLE = "sample";

        TestMultilevelMetrics sampleLibrary = new TestMultilevelMetrics();
        sampleLibrary.SAMPLE = "sample";
        sampleLibrary.LIBRARY = "libraryDuped";

        return new Object[][] {
                new Object[] { Arrays.asList(sampleOnly, sampleOnly), null, sampleOnly, 2 },
                new Object[] { Arrays.asList(sampleLibrary, sampleLibrary), "libraryDuped", sampleLibrary, 2 },
                new Object[] { Arrays.asList(sampleLibrary, sampleOnly, sampleLibrary), "libraryDuped", sampleLibrary, 2 },
                new Object[] { Arrays.asList(sampleOnly, sampleLibrary, sampleOnly, sampleLibrary, sampleOnly), "libraryDuped", sampleLibrary, 2 },
        };
    }

    @Test(dataProvider = "duplicatedMetrics", expectedExceptions = IllegalStateException.class)
    public void testErrorOnDuplicatedMetric(List<TestMultilevelMetrics> metrics, String libraryName, TestMultilevelMetrics libraryMetric, @SuppressWarnings("unused") int libraryCount) {
        Assert.assertSame(MetricsUtils.getLibraryMetric(metrics, libraryName), libraryMetric, "getting metric with library set to " + formatName(libraryName));
    }

    @Test(dataProvider = "duplicatedMetrics")
    public void testFilterDuplicatedMetrics(List<TestMultilevelMetrics> metrics, String libraryName, TestMultilevelMetrics libraryMetric, int libraryCount) {
        List<TestMultilevelMetrics> actuals = MetricsUtils.filterMetrics(metrics, "LIBRARY", libraryName);
        Assert.assertEquals(actuals.size(), libraryCount, "filtering metrics count with library set to " + formatName(libraryName));
        for (TestMultilevelMetrics actual: actuals)
            Assert.assertSame(actual, libraryMetric, "filtering metrics with library set to " + formatName(libraryName));
    }

    private static String formatName(String name) {
        return name == null ? "<null>" : String.format("'%s'", name);
    }
}
