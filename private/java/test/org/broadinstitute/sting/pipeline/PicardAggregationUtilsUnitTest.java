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

import org.apache.commons.io.FileUtils;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

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
        File sampleBam = new File(PicardAggregationUtils.getSampleBam(PROJECT, SAMPLE, latestVersion));
        Assert.assertTrue(sampleBam.exists());

        File latestSampleBam = new File(PicardAggregationUtils.getSampleBam(PROJECT, SAMPLE));
        Assert.assertTrue(latestSampleBam.exists());
        Assert.assertEquals(latestSampleBam, sampleBam);

        File currentSampleBam = new File(PicardAggregationUtils.getCurrentSampleBam(PROJECT, SAMPLE));
        Assert.assertTrue(currentSampleBam.exists());
        Assert.assertEquals(currentSampleBam.getCanonicalFile(), sampleBam.getCanonicalFile());
    }

    @Test(dependsOnMethods = "testGetLatestVersion")
    public void testGetSampleDir() throws Exception {
        File sampleDir = new File(PicardAggregationUtils.getSampleDir(PROJECT, SAMPLE, latestVersion));
        Assert.assertTrue(sampleDir.exists());

        File latestSampleDir = new File(PicardAggregationUtils.getSampleDir(PROJECT, SAMPLE));
        Assert.assertTrue(latestSampleDir.exists());
        Assert.assertEquals(latestSampleDir, sampleDir);

        File currentSampleDir = new File(PicardAggregationUtils.getCurrentSampleDir(PROJECT, SAMPLE));
        Assert.assertTrue(currentSampleDir.exists());
        Assert.assertEquals(currentSampleDir.getCanonicalFile(), sampleDir.getCanonicalFile());
    }

    @Test(dependsOnMethods = "testGetLatestVersion")
    public void testIsFinished() {
        Assert.assertTrue(PicardAggregationUtils.isFinished(PROJECT, SAMPLE, latestVersion));
        Assert.assertFalse(PicardAggregationUtils.isFinished(PROJECT, SAMPLE, latestVersion + 1));
    }

    @Test
    public void testLatestVersionMissing() {
        Assert.assertEquals(PicardAggregationUtils.getLatestVersion(MISSING_PROJECT, MISSING_SAMPLE), 0);
    }

    @Test
    public void testSafeFileNames() throws FileNotFoundException {
        int slashLatest = PicardAggregationUtils.getLatestVersion(SLASH_PROJECT, SLASH_SAMPLE);
        int spaceLatest = PicardAggregationUtils.getLatestVersion(SPACE_PROJECT, SPACE_SAMPLE);
        Assert.assertTrue(slashLatest > 0);
        Assert.assertTrue(spaceLatest > 0);
        PicardAggregationUtils.getSampleDir(SLASH_PROJECT, SLASH_SAMPLE, slashLatest);
        PicardAggregationUtils.getSampleBam(SLASH_PROJECT, SLASH_SAMPLE, slashLatest);
        PicardAggregationUtils.getSampleDir(SPACE_PROJECT, SPACE_SAMPLE, spaceLatest);
        PicardAggregationUtils.getSampleBam(SPACE_PROJECT, SPACE_SAMPLE, spaceLatest);
    }

    @Test
    public void testParseTsv() throws IOException {
        File tsv = writeTsv(PROJECT, SAMPLE);
        List<PicardSample> picardSamples = PicardAggregationUtils.parseSamples(tsv);
        Assert.assertEquals(picardSamples.size(), 1);

        PicardIntervals intervals = PicardAggregationUtils.readAnalysisIntervals(picardSamples);
        Assert.assertEquals(intervals.getReference(), BaseTest.hg19Reference);
        // BaseTest.hg19Intervals and the test are the same intervals, but in a different location. For now just testing the name.
        Assert.assertEquals(new File(intervals.getTargets()).getName(), new File(BaseTest.hg19Intervals).getName());

        List<String> bams = PicardAggregationUtils.getSampleBams(picardSamples);
        for (String bam: bams)
            Assert.assertTrue(new File(bam).exists(), "bam does not exist: " + bam);
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testParseThrowOnBadTsv() throws IOException {
        File tsv = writeTsv(MISSING_PROJECT, MISSING_SAMPLE);
        List<PicardSample> picardSamples = PicardAggregationUtils.parseSamples(tsv);
    }

    @Test
    public void testParseCatchOnBadTsv() throws IOException {
        File tsv = writeTsv(MISSING_PROJECT, MISSING_SAMPLE);
        List<PicardSample> picardSamples = PicardAggregationUtils.parseSamples(tsv, false);
        Assert.assertEquals(picardSamples.size(), 0);
    }

    @Test
    public void testParseTsvWithPicardComments() throws Exception {
        File tsv = writeTsv("C460", "HG01359");
        List<PicardSample> picardSamples = PicardAggregationUtils.parseSamples(tsv);

        PicardIntervals intervals = PicardAggregationUtils.readAnalysisIntervals(picardSamples);
        Assert.assertEquals(intervals.getReference(), BaseTest.hg19Reference);
        // BaseTest.hg19Intervals and the test are the same intervals, but in a different location. For now just testing the name.
        Assert.assertEquals(new File(intervals.getTargets()).getName(), new File(BaseTest.hg19Intervals).getName());

        List<String> bams = PicardAggregationUtils.getSampleBams(picardSamples);
        for (String bam: bams)
            Assert.assertTrue(new File(bam).exists(), "bam does not exist: " + bam);
    }

    private File writeTsv(String project, String sample) throws IOException {
        File tsv = BaseTest.createTempFile("pipeline", ".tsv");
        FileUtils.writeLines(tsv, Collections.singletonList(project + "\t" + sample));
        return tsv;
    }
}
