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

package org.broadinstitute.sting.utils.bcf2;

/**
* [Short one sentence description of this walker]
* <p/>
* <p>
* [Functionality of this walker]
* </p>
* <p/>
* <h2>Input</h2>
* <p>
* [Input description]
* </p>
* <p/>
* <h2>Output</h2>
* <p>
* [Output description]
* </p>
* <p/>
* <h2>Examples</h2>
* <pre>
*    java
*      -jar GenomeAnalysisTK.jar
*      -T $WalkerName
*  </pre>
*
* @author Your Name
* @since Date created
*/
public enum BCFType {
    RESERVED_0,
    INT8(1, -127, 128), // todo -- confirm range
    INT16(2, Short.MIN_VALUE, Short.MAX_VALUE),
    INT32(4, Integer.MIN_VALUE, Integer.MAX_VALUE),
    INT64(8, Long.MIN_VALUE, Long.MAX_VALUE),
    FLOAT(4),
    DOUBLE(8),
    CHAR(1),
    FLAG(1),
    STRING_LITERAL,
    STRING_REF8(1, -127, 128), // todo -- confirm range
    STRING_REF16(2, Short.MIN_VALUE, Short.MAX_VALUE),
    STRING_REF32(4, Integer.MIN_VALUE, Integer.MAX_VALUE),
    COMPACT_GENOTYPE(2),
    RESERVED_14,
    RESERVED_15;

    private final int sizeInBytes;
    private final long minValue, maxValue;

    BCFType() {
        this(-1, 0, 0);
    }

    BCFType(final int sizeInBytes) {
        this(sizeInBytes, 0, 0);
    }

    BCFType(final int sizeInBytes, final long minValue, final long maxValue) {
        this.sizeInBytes = sizeInBytes;
        this.minValue = minValue;
        this.maxValue = maxValue;
    }

    public int getSizeInBytes() {
        return sizeInBytes;
    }

    public int getID() { return ordinal(); }

    public final boolean withinRange(final long v) {
        return v >= minValue && v <= maxValue;
    }
}
