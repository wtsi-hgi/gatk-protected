/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.gatk.utils.haplotype;

import java.util.Comparator;

/**
 * Compares two haplotypes first by their lengths and then by lexicographic order of their bases.
 *
 * User: btaylor
 * Date: 8/1/13
 * Time: 11:09 AM
 */
public class HaplotypeSizeAndBaseComparator implements Comparator<Haplotype> {
    @Override
    public int compare( final Haplotype hap1, final Haplotype hap2 ) {
        if (hap1.getBases().length < hap2.getBases().length)
            return -1;
        else if (hap1.getBases().length > hap2.getBases().length)
            return 1;
        else
            return hap1.getBaseString().compareTo(hap2.getBaseString());
    }
}
