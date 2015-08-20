/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.readutils;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;

/**
 * Created with IntelliJ IDEA.
 * User: delangel
 * Date: 4/13/13
 * Time: 7:28 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReadAdaptorTrimmerIntegrationTest extends WalkerTest {
    private String getBaseCommand(final String BAM) {
        return  "-T ReadAdaptorTrimmer -R " + b37KGReference +
                " -I " + privateTestDir + BAM +
                " -o %s";
    }

    @Test
    public void testBasicTrimmer() {
        WalkerTestSpec spec = new WalkerTestSpec( getBaseCommand("shortInsertTest.bam"),  1, Arrays.asList("c7d7f69e6b532ec693bfbd821c2e9766"));
        executeTest(String.format("testBasicTrimmer"), spec);
    }

    @Test
    public void testSkippingBadPairs() {
        WalkerTestSpec spec = new WalkerTestSpec( getBaseCommand("shortInsertTest2.bam")+" -removeUnpairedReads",  1, Arrays.asList("f7aa76d1a2535774764e06ba610c21de"));
        executeTest(String.format("testSkippingBadPairs"), spec);
    }

}
