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

package org.broadinstitute.sting.gatk.datasources.reads.performance;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Feb 25, 2011
 * Time: 10:16:53 AM
 * To change this template use File | Settings | File Templates.
 */
class IterateOverEachBase extends ReadProcessor {
    private long As;
    private long Cs;
    private long Gs;
    private long Ts;

    public IterateOverEachBase(final BAMProcessingPerformanceMeter performanceMeter) {
        super(performanceMeter);
    }

    @Override
    public String getTestName() { return "iterate over each base"; }
    public void processRead(final SAMRecord read) {
        for(byte base: read.getReadBases()) {
            switch(base) {
                case 'A': As++; break;
                case 'C': Cs++; break;
                case 'G': Gs++; break;
                case 'T': Ts++; break;
            }
        }
    }
}