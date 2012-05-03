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

import java.nio.ByteOrder;
import java.nio.charset.Charset;
import java.util.EnumSet;

public class BCF2Constants {

    public static final String VERSION_LINE_FORMAT = "fileformat=BCF2v%d.%d";
    public static final String VERSION_LINE = String.format(VERSION_LINE_FORMAT, 0, 1);
    public static final String DICTIONARY_LINE_FORMAT = "dictionary=%s";
    public static final String DICTIONARY_LINE_ENTRY_SEPARATOR = ",";

    public static final ByteOrder BCF2_BYTE_ORDER = ByteOrder.LITTLE_ENDIAN;
    public static final Charset BCF2_TEXT_CHARSET = Charset.forName("US-ASCII");  // TODO: enforce this!

    public static final int INT8_TYPE_ID =             0x01;
    public static final int INT16_TYPE_ID =            0x02;
    public static final int INT32_TYPE_ID =            0x03;
    public static final int INT64_TYPE_ID =            0x04;
    public static final int FLOAT_TYPE_ID =            0x05;
    public static final int DOUBLE_TYPE_ID =           0x06;
    public static final int CHAR_TYPE_ID =             0x07;
    public static final int FLAG_TYPE_ID =             0x08;
    public static final int STRING_LITERAL_TYPE_ID =   0x09;
    public static final int STRING_REF8_TYPE_ID =      0x0A;
    public static final int STRING_REF16_TYPE_ID =     0x0B;
    public static final int STRING_REF32_TYPE_ID =     0x0C;
    public static final int COMPACT_GENOTYPE_TYPE_ID = 0x0D;
    public static final int MIXED_VECTOR_TYPE_ID =     0x0E;

    public static final int TYPE_DESCRIPTOR_ATOMICITY_FIELD_MASK =    0x80;
    public static final int TYPE_DESCRIPTOR_NUM_ELEMENTS_FIELD_MASK = 0x40 | 0x20 | 0x10;
    public static final int TYPE_DESCRIPTOR_TYPE_ID_MASK =            0x08 | 0x04 | 0x02 | 0x01;

    public static final int INT8_MISSING_VALUE =   0x80;
    public static final int INT16_MISSING_VALUE =  0x8000;
    public static final int INT32_MISSING_VALUE =  0x80000000;
    public static final long INT64_MISSING_VALUE = 0x8000000000000000L;

    public static final int FLOAT_NAN_VALUE =     0x7FC00000;
    public static final int FLOAT_MISSING_VALUE = 0x7F800001;

    public static final long DOUBLE_NAN_VALUE =     0x7FF8000000000000L;
    public static final long DOUBLE_MISSING_VALUE = 0x7FF0000000000001L;

    public static final int FLAG_ABSENT_VALUE =  0x00;
    public static final int FLAG_PRESENT_VALUE = 0x01;

    public static final byte STRING_TERMINATOR = 0x00;

    public static final int GENOTYPE_HIGH_ORDER_BIT_FIELD_MASK = 0x80;
    public static final int GENOTYPE_PHASED_FIELD_MASK =         0x40;
    public static final int GENOTYPE_FIRST_ALLELE_FIELD_MASK =   0x20 | 0x10 | 0x08;
    public static final int GENOTYPE_SECOND_ALLELE_FIELD_MASK =  0x04 | 0x02 | 0x01;

    public static final byte GENOTYPE_HIGH_ORDER_BIT_VALUE =    0x01;
    public static final byte GENOTYPE_MISSING_VALUE_INDICATOR = 0x07;

    public static final EnumSet<BCF2Type> ATOMIC_INTEGRAL_TYPES = EnumSet.of(BCF2Type.INT8,
                                                                             BCF2Type.INT16,
                                                                             BCF2Type.INT32,
                                                                             BCF2Type.INT64);

    public static final EnumSet<BCF2Type> ATOMIC_STRING_TYPES = EnumSet.of(BCF2Type.STRING_LITERAL,
                                                                           BCF2Type.STRING_REF8,
                                                                           BCF2Type.STRING_REF16,
                                                                           BCF2Type.STRING_REF32);

    public static final EnumSet<BCF2Type> VECTOR_STRING_TYPES = EnumSet.of(BCF2Type.VECTOR_STRING_LITERAL,
                                                                           BCF2Type.VECTOR_STRING_REF8,
                                                                           BCF2Type.VECTOR_STRING_REF16,
                                                                           BCF2Type.VECTOR_STRING_REF32);

    public static final EnumSet<BCF2Type> ALL_STRING_TYPES = EnumSet.of(BCF2Type.STRING_LITERAL,
                                                                        BCF2Type.STRING_REF8,
                                                                        BCF2Type.STRING_REF16,
                                                                        BCF2Type.STRING_REF32,
                                                                        BCF2Type.VECTOR_STRING_LITERAL,
                                                                        BCF2Type.VECTOR_STRING_REF8,
                                                                        BCF2Type.VECTOR_STRING_REF16,
                                                                        BCF2Type.VECTOR_STRING_REF32);

    public static final EnumSet<BCF2Type> ALL_TYPES = EnumSet.allOf(BCF2Type.class);
}
