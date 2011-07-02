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

import java.io.File;

public class ReferenceDataUnitTest {
    @Test
    public void testNames() {
        Assert.assertEquals(ReferenceData.HG18.getName(), "hg18");
        Assert.assertEquals(ReferenceData.HG19.getName(), "hg19");
    }

    @Test
    public void testFilesExist() {
        for (ReferenceData data: ReferenceData.values()) {
            Assert.assertTrue(new File(data.getReference()).exists());
            Assert.assertTrue(new File(data.getRefseq()).exists());
            for (int version: data.getDbsnpVersions()) {
                Assert.assertTrue(new File(data.getDbsnp(version)).exists());
            }
        }
    }

    @Test
    public void testDbsnps() {
        Assert.assertTrue(new File(ReferenceData.HG18.getDbsnp(129)).exists());
        Assert.assertTrue(new File(ReferenceData.HG19.getDbsnp(129)).exists());
        Assert.assertTrue(new File(ReferenceData.HG19.getDbsnp(132)).exists());
        Assert.assertNull(ReferenceData.HG19.getDbsnp(130));
    }

    @Test
    public void testDbsnpTypes() {
        Assert.assertEquals(ReferenceData.HG18.getDbsnpType(129), "DBSNP");
        Assert.assertEquals(ReferenceData.HG19.getDbsnpType(129), "VCF");
        Assert.assertEquals(ReferenceData.HG19.getDbsnpType(132), "VCF");
        Assert.assertNull(ReferenceData.HG19.getDbsnpType(130));
    }

    @Test
    public void testGetByReference() {
        Assert.assertEquals(ReferenceData.getByReference(BaseTest.hg18Reference), ReferenceData.HG18);
        Assert.assertEquals(ReferenceData.getByReference(BaseTest.hg19Reference), ReferenceData.HG19);
        Assert.assertEquals(ReferenceData.getByReference("none"), null);
    }
}
