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
import org.broadinstitute.sting.utils.classloader.JVMUtils;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;

/**
 * Utility methods for Picard Metrics
 */
public class MetricsUtils {
    /**
     * Returns the sample metric or null if it doesn't exist.
     * @param metrics List of metrics
     * @param <BEAN> Type of the metric
     * @return Sample metric or null
     * @throws IllegalStateException If more than one sample metric was found in the list
     */
    public static <BEAN extends MetricBase> BEAN getSampleMetric(List<BEAN> metrics) {
        return getLibraryMetric(metrics, null);
    }

    /**
     * Returns the library metric or null if it doesn't exist.
     * @param metrics List of metrics
     * @param library The name of the library or null to return the sample metric.
     * @param <BEAN> Type of the metric
     * @return Library metric or null
     * @throws IllegalStateException If more than one metric was found in the list matching library
     */
    public static <BEAN extends MetricBase> BEAN getLibraryMetric(List<BEAN> metrics, String library) {
        List<BEAN> filteredMetrics = filterMetrics(metrics, "LIBRARY", library);
        if (filteredMetrics.size() > 1)
            throw new IllegalStateException(String.format("%d metrics returned", filteredMetrics.size()));
        if (filteredMetrics.size() == 0)
            return null;
        return filteredMetrics.get(0);
    }

    /**
     * Filters a list of metrics by field name and value.
     * @param metrics List of metrics
     * @param fieldName Field name
     * @param value Field value
     * @param <BEAN> Type of the metric
     * @return List of filtered metrics
     */
    public static <BEAN extends MetricBase> List<BEAN> filterMetrics(List<BEAN> metrics, String fieldName, String value) {
        List<BEAN> filteredMetrics = new ArrayList<BEAN>();
        for (BEAN metric : metrics) {
            Field field = JVMUtils.findField(metric.getClass(), fieldName);
            if (field != null) {
                Object fieldValue = JVMUtils.getFieldValue(field, metric);
                if (value == null) {
                    if (fieldValue == null)
                        filteredMetrics.add(metric);
                } else if (fieldValue != null) {
                    if (value.equals(fieldValue.toString()))
                        filteredMetrics.add(metric);
                }
            }
        }
        return filteredMetrics;
    }

    /**
     * Returns the field names in the metrics class that don't start with a double underscore "__".
     * @param metricClass Class of the metric
     * @param <BEAN> Class of the metric
     * @return List of field names
     */
    public static <BEAN extends MetricBase> List<String> getFieldNames(Class<BEAN> metricClass) {
        ArrayList<String> fieldNames = new ArrayList<String>();
        for (Field field: metricClass.getFields()) {
            String name = field.getName();
            if (!name.startsWith("__"))
                fieldNames.add(name);
        }
        return fieldNames;
    }
}
