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

package org.broadinstitute.sting.pipeline;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.FileNotFoundException;

public class PicardAggregationUtilsUnitTest {
    public static final String PROJECT = "C474";
    public static final String SAMPLE = "NA19651";
    public static final String MISSING_PROJECT = "C0";
    public static final String MISSING_SAMPLE = "0";
    public static final String SLASH_PROJECT = "C279";
    public static final String SLASH_SAMPLE = "MTI-578-3-4/MTI-1015";
    public static final String SPACE_PROJECT = "C507";
    public static final String SPACE_SAMPLE = "FG-CR 6";
    private int latestVersion = -1;

    @Test
    public void testGetLatestVersion() {
        latestVersion = PicardAggregationUtils.getLatestVersion(PROJECT, SAMPLE);
        System.out.println(String.format("Latest version for %s %s is %d", PROJECT, SAMPLE, latestVersion));
        Assert.assertTrue(latestVersion > 0);
    }

    @Test(dependsOnMethods = "testGetLatestVersion")
    public void testGetSampleBam() throws Exception {
        String test = PicardAggregationUtils.getSampleBam(PROJECT, SAMPLE);
        String latest = PicardAggregationUtils.getSampleBam(PROJECT, SAMPLE, latestVersion);
        Assert.assertEquals(test, latest);
    }

    @Test(dependsOnMethods = "testGetLatestVersion")
    public void testGetSampleDir() throws Exception {
        String test = PicardAggregationUtils.getSampleDir(PROJECT, SAMPLE);
        String latest = PicardAggregationUtils.getSampleDir(PROJECT, SAMPLE, latestVersion);
        Assert.assertEquals(test, latest);
    }

    @Test(dependsOnMethods = "testGetLatestVersion")
    public void testIsFinished() {
        Assert.assertTrue(PicardAggregationUtils.isFinished(PROJECT, SAMPLE, latestVersion));
        Assert.assertFalse(PicardAggregationUtils.isFinished(PROJECT, SAMPLE, latestVersion + 1));
    }

    @Test(expectedExceptions = FileNotFoundException.class)
    public void testMissingSampleBam() throws Exception {
        PicardAggregationUtils.getSampleBam(MISSING_PROJECT, MISSING_SAMPLE);
    }

    @Test(expectedExceptions = FileNotFoundException.class)
    public void testMissingSampleDir() throws Exception {
        PicardAggregationUtils.getSampleDir(MISSING_PROJECT, MISSING_SAMPLE);
    }

    @Test
    public void testLatestVersionMissing() {
        Assert.assertEquals(PicardAggregationUtils.getLatestVersion(MISSING_PROJECT, MISSING_SAMPLE), 0);
    }

    @Test
    public void testSafeFileNames() throws FileNotFoundException {
        Assert.assertTrue(PicardAggregationUtils.getLatestVersion(SLASH_PROJECT, SLASH_SAMPLE) > 0);
        Assert.assertTrue(PicardAggregationUtils.getLatestVersion(SPACE_PROJECT, SPACE_SAMPLE) > 0);
        PicardAggregationUtils.getSampleDir(SLASH_PROJECT, SLASH_SAMPLE);
        PicardAggregationUtils.getSampleBam(SLASH_PROJECT, SLASH_SAMPLE);
        PicardAggregationUtils.getSampleDir(SPACE_PROJECT, SPACE_SAMPLE);
        PicardAggregationUtils.getSampleBam(SPACE_PROJECT, SPACE_SAMPLE);
    }
}
