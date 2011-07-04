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

// our package
package org.broadinstitute.sting.gatk.walkers.diffengine;


// the imports for unit testing.


import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import javax.management.StringValueExp;
import java.io.File;
import java.util.*;

/**
 * Basic unit test for DifferableReaders in reduced reads
 */
public class DiffableReaderUnitTest extends BaseTest {
    DiffEngine engine;
//    private class SingleTest {
//        public String bases;
//        public byte mostCountBase;
//        public int mostCommonCount;
//
//        private SingleTest(String bases, char mostCountBase, int mostCommonCount) {
//            this.mostCommonCount = mostCommonCount;
//            this.mostCountBase = (byte)mostCountBase;
//            this.bases = bases;
//        }
//    }
//
//    @DataProvider(name = "data")
//    public Object[][] createData1() {
//        List<SingleTest> params = new ArrayList<SingleTest>();
//
//        params.add(new SingleTest("A", 'A', 1 ));
//        params.add(new SingleTest("AA", 'A', 2 ));
//        params.add(new SingleTest("AC", 'A', 1 ));
//        params.add(new SingleTest("AAC", 'A', 2 ));
//        params.add(new SingleTest("AAA", 'A', 3 ));
//        params.add(new SingleTest("AAAN", 'A', 3 ));
//        params.add(new SingleTest("AAANNNN", 'A', 3 ));
//        params.add(new SingleTest("AACTG", 'A', 2 ));
//        params.add(new SingleTest("D", 'D', 1 ));
//        params.add(new SingleTest("DDAAD", 'D', 3));
//        params.add(new SingleTest("", (char)BaseCounts.MAX_BASE_WITH_NO_COUNTS, 0 ));
//        params.add(new SingleTest("AAIIIAI", 'I', 4 ));
//
//        List<Object[]> params2 = new ArrayList<Object[]>();
//        for ( SingleTest x : params ) params2.add(new Object[]{x});
//        return params2.toArray(new Object[][]{});
//    }

    File vcfFile = new File(testDir + "diffTestMaster.vcf");
    File bamFile = new File(testDir + "exampleBAM.bam");

    @BeforeClass(enabled = true)
    public void createDiffEngine() {
        engine = new DiffEngine(10);
    }

    @Test(enabled = true)
    public void testPluggableDiffableReaders() {
        logger.warn("testPluggableDiffableReaders");
        Map<String, DiffableReader> readers = engine.getReaders();
        Assert.assertNotNull(readers);
        Assert.assertTrue(readers.size() > 0);
        Assert.assertNotNull(readers.get("VCF"));
        for ( Map.Entry<String, DiffableReader> e : engine.getReaders().entrySet() ) {
            logger.warn("Found diffable reader: " + e.getKey());
            Assert.assertEquals(e.getValue().getName(), e.getKey());
            Assert.assertEquals(e.getValue(), engine.getReader(e.getKey()));
        }
    }

    private static void testLeaf(DiffNode rec, String field, Object expected) {
        DiffLeaf leaf = rec.getLeaf(field);
        Assert.assertNotNull(leaf, "Expected to see leaf named " + field + " in rec " + rec);
        Assert.assertEquals(leaf.getValue(), expected, "Expected to leaf named " + field + " to have value " + expected + " in rec " + rec);
    }

    @Test(enabled = true, dependsOnMethods = "testPluggableDiffableReaders")
    public void testVCF1() {
        logger.warn("testVCF1");
        DiffableReader vcfReader = engine.getReader("VCF");
        Assert.assertTrue(vcfReader.canRead(vcfFile));
        Assert.assertFalse(vcfReader.canRead(bamFile));

        DiffNode diff = vcfReader.readFromFile(vcfFile);
        Assert.assertNotNull(diff);
        //logger.warn(diff);

        Assert.assertEquals(diff.getName(), vcfFile.getName());
        Assert.assertSame(diff.getParent(), DiffNode.ROOT);
        Assert.assertEquals(diff.getLeaves().size(), 0);
        Assert.assertEquals(diff.getNodes().size(), 9);

        // chr1    2646    rs62635284      G       A       0.15    PASS    AC=2;AF=1.00;AN=2       GT:AD:DP:GL:GQ  1/1:53,75:3:-12.40,-0.90,-0.00:9.03
        DiffNode rec1 = diff.getNodes().get(0);
        testLeaf(rec1, "CHROM", "chr1");
        testLeaf(rec1, "POS", 2646);
        testLeaf(rec1, "ID", "rs62635284");
        testLeaf(rec1, "REF", Allele.create("G", true));
        testLeaf(rec1, "ALT", new HashSet<Allele>(Arrays.asList(Allele.create("A"))));
        testLeaf(rec1, "QUAL", 0.15);
        testLeaf(rec1, "FILTER", Collections.<Object>emptySet());
        testLeaf(rec1, "AC", "2");
        testLeaf(rec1, "AF", "1.00");
        testLeaf(rec1, "AN", "2");
    }
}