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

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.FileNotFoundException;

import static org.broadinstitute.sting.pipeline.PicardAggregationUtilsUnitTest.*;

public class PicardAnalysisFilesUnitTest extends BaseTest {
    @Test
    public void testParseLatest() throws Exception {
        PicardAnalysisFiles files = new PicardAnalysisFiles(PROJECT, SAMPLE);
        Assert.assertNotNull(files.getPath());
        files = new PicardAnalysisFiles(PROJECT, SAMPLE, PicardAggregationUtils.getLatestVersion(PROJECT, SAMPLE));
        Assert.assertNotNull(files.getPath());
    }

    @Test
    public void testParseValid() throws Exception {
        PicardAnalysisFiles file = new PicardAnalysisFiles(BaseTest.validationDataLocation + "picard_analysis_file.txt");
        Assert.assertEquals(file.getReferenceSequence(), "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta");
        Assert.assertEquals(file.getTargetIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list");
        Assert.assertEquals(file.getBaitIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.interval_list");
    }

    @Test
    public void testParseValidWithComments() throws Exception {
        PicardAnalysisFiles file = new PicardAnalysisFiles(BaseTest.validationDataLocation + "picard_analysis_file_with_comments.txt");
        Assert.assertEquals(file.getReferenceSequence(), "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta");
        Assert.assertEquals(file.getTargetIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list");
        Assert.assertEquals(file.getBaitIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.interval_list");
    }

    @Test
    public void testParseValidWithoutReference() throws Exception {
        PicardAnalysisFiles file = new PicardAnalysisFiles(BaseTest.validationDataLocation + "picard_analysis_file_without_reference.txt");
        Assert.assertNull(file.getReferenceSequence());
        Assert.assertEquals(file.getTargetIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list");
        Assert.assertEquals(file.getBaitIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.interval_list");
    }

    @Test
    public void testParseValidWithBlankLine() throws Exception {
        PicardAnalysisFiles file = new PicardAnalysisFiles(BaseTest.validationDataLocation + "picard_analysis_file_with_blank_line.txt");
        Assert.assertNull(file.getReferenceSequence());
        Assert.assertEquals(file.getTargetIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list");
        Assert.assertEquals(file.getBaitIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.interval_list");
    }

    @Test(expectedExceptions = FileNotFoundException.class)
    public void testParseBadPath() throws Exception {
        new PicardAnalysisFiles(BaseTest.validationDataLocation + "non_existent_picard_analysis_file.txt");
    }

    @Test(expectedExceptions = FileNotFoundException.class)
    public void testParseMissingLatest() throws Exception {
        new PicardAnalysisFiles(MISSING_PROJECT, MISSING_SAMPLE);
    }

    @Test(expectedExceptions = FileNotFoundException.class)
    public void testParseMissingVersion() throws Exception {
        new PicardAnalysisFiles(PROJECT, SAMPLE, PicardAggregationUtils.getLatestVersion(PROJECT, SAMPLE) + 2);
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testParseMultipleReferences() throws Exception {
        PicardAnalysisFiles file = new PicardAnalysisFiles(BaseTest.validationDataLocation + "picard_analysis_file_with_different_references.txt");
        file.getReferenceSequence();
    }
}
