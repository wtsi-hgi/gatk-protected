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

import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

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
    INT8(1, BCF2Constants.INT8_MISSING_VALUE, -127, 128), // todo -- confirm range
    INT16(2, BCF2Constants.INT16_MISSING_VALUE, Short.MIN_VALUE, Short.MAX_VALUE),
    INT32(4, BCF2Constants.INT32_MISSING_VALUE, Integer.MIN_VALUE, Integer.MAX_VALUE),
    INT64_NOT_USED, // (8, Long.MIN_VALUE, Long.MAX_VALUE),
    FLOAT(4, BCF2Constants.FLOAT_MISSING_VALUE),
    DOUBLE_NOT_USED, // (8),
    CHAR_NOT_USED, // (1),
    FLAG(1),
    STRING_LITERAL,
    STRING_REF8(1, null, -127, 128), // todo -- confirm range
    STRING_REF16(2, null, Short.MIN_VALUE, Short.MAX_VALUE),
    STRING_REF32_NOT_USED, // (4, Integer.MIN_VALUE, Integer.MAX_VALUE),
    COMPACT_GENOTYPE(1),
    RESERVED_14,
    RESERVED_15;

    private final Object missingValue;
    private final int sizeInBytes;
    private final long minValue, maxValue;

    BCFType() {
        this(-1, null, 0, 0);
    }

    BCFType(final int sizeInBytes) {
        this(sizeInBytes, null, 0, 0);
    }

    BCFType(final int sizeInBytes, final Object missingValue) {
        this(sizeInBytes, missingValue, 0, 0);
    }

    BCFType(final int sizeInBytes, final Object missingValue, final long minValue, final long maxValue) {
        this.sizeInBytes = sizeInBytes;
        this.missingValue = missingValue;
        this.minValue = minValue;
        this.maxValue = maxValue;
    }

    public int getSizeInBytes() {
        return sizeInBytes;
    }
    public int getID() { return ordinal(); }
    public final boolean withinRange(final long v) { return v >= minValue && v <= maxValue; }
    public Object getMissingValue() { return missingValue; }

    public Object convertIfNecessary(final Object maybeString) {
        if ( maybeString instanceof String ) {
            final String string = (String)maybeString;
            if ( string.equals(VCFConstants.MISSING_VALUE_v4) )
                return getMissingValue();
            else {
                switch ( this ) {
                    case INT8:
                    case INT16:
                    case INT32:  return Integer.valueOf(string);
                    case FLOAT:  return Float.valueOf(string);
                    case FLAG:   return Boolean.valueOf(string);
                    case STRING_LITERAL:
                    case STRING_REF8:
                    case STRING_REF16: return maybeString;
                    default: throw new ReviewedStingException("Cannot convert value to type" + this);
                }
            }
        } else
            return maybeString;
    }
}
