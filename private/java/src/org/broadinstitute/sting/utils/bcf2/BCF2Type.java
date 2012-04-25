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

import java.nio.ByteBuffer;
import java.util.ArrayList;

public enum BCF2Type {
    INT8 (BCF2Constants.INT8_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            byte value = data.get();
            return new BCF2Value(this, value, value == BCF2Constants.INT8_MISSING_VALUE);
        }
    },

    INT16 (BCF2Constants.INT16_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            short value = data.getShort();
            return new BCF2Value(this, value, value == BCF2Constants.INT16_MISSING_VALUE);
        }
    },

    INT32 (BCF2Constants.INT32_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            int value = data.getInt();
            return new BCF2Value(this, value, value == BCF2Constants.INT32_MISSING_VALUE);
        }
    },

    INT64 (BCF2Constants.INT64_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            long value = data.getLong();
            return new BCF2Value(this, value, value == BCF2Constants.INT64_MISSING_VALUE);
        }
    },

    FLOAT (BCF2Constants.FLOAT_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            int floatBitsAsInt = data.getInt();
            return new BCF2Value(this, Float.intBitsToFloat(floatBitsAsInt), floatBitsAsInt == BCF2Constants.FLOAT_MISSING_VALUE);
        }
    },

    DOUBLE (BCF2Constants.DOUBLE_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            long doubleBitsAsLong = data.getLong();
            return new BCF2Value(this, Double.longBitsToDouble(doubleBitsAsLong), doubleBitsAsLong == BCF2Constants.DOUBLE_MISSING_VALUE);
        }
    },

    FLAG (BCF2Constants.FLAG_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return null; // TODO
        }
    },

    STRING_LITERAL (BCF2Constants.STRING_LITERAL_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            StringBuilder builder = new StringBuilder();
            byte currentByte = data.get();
            while ( currentByte != BCF2Constants.STRING_TERMINATOR ) {
                builder.append((char)currentByte);     // TODO: enforce ASCII charset encoding
                currentByte = data.get();
            }
            String decodedString = builder.toString();

            return new BCF2Value(this, decodedString, decodedString.length() == 0);
        }
    },

    STRING_REF8 (BCF2Constants.STRING_REF8_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return null; // TODO
        }
    },

    STRING_REF16 (BCF2Constants.STRING_REF16_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return null; // TODO
        }
    },

    STRING_REF32 (BCF2Constants.STRING_REF32_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return null; // TODO
        }
    },

    COMPACT_GENOTYPE (BCF2Constants.COMPACT_GENOTYPE_TYPE_ID, true) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return null; // TODO
        }
    },

    VECTOR_INT8 (BCF2Constants.INT8_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, INT8);
        }
    },

    VECTOR_INT16 (BCF2Constants.INT16_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, INT16);
        }
    },

    VECTOR_INT32 (BCF2Constants.INT32_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, INT32);
        }
    },

    VECTOR_INT64 (BCF2Constants.INT64_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, INT64);
        }
    },

    VECTOR_FLOAT (BCF2Constants.FLOAT_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, FLOAT);
        }
    },

    VECTOR_DOUBLE (BCF2Constants.DOUBLE_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, DOUBLE);
        }
    },

    VECTOR_FLAG (BCF2Constants.FLAG_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, FLAG);
        }
    },

    VECTOR_STRING_LITERAL (BCF2Constants.STRING_LITERAL_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, STRING_LITERAL);
        }
    },

    VECTOR_STRING_REF8 (BCF2Constants.STRING_REF8_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, STRING_REF8);
        }
    },

    VECTOR_STRING_REF16 (BCF2Constants.STRING_REF16_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, STRING_REF16);
        }
    },

    VECTOR_STRING_REF32 (BCF2Constants.STRING_REF32_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, STRING_REF32);
        }
    },

    VECTOR_COMPACT_GENOTYPE (BCF2Constants.COMPACT_GENOTYPE_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return decodeVector(data, numItems, COMPACT_GENOTYPE);
        }
    },

    MIXED_VECTOR (BCF2Constants.MIXED_VECTOR_TYPE_ID, false) {
        public BCF2Value decode( ByteBuffer data, int numItems ) {
            return null; // TODO
        }
    };

    protected int typeID;
    protected boolean isAtomic;

    BCF2Type( int typeID, boolean isAtomic ) {
        this.typeID = typeID;
        this.isAtomic = isAtomic;
    }

    public int getTypeID() {
        return typeID;
    }

    public boolean isAtomic() {
        return isAtomic;
    }

    public abstract BCF2Value decode( ByteBuffer data, int numItems );

    protected BCF2Value decodeVector( ByteBuffer data, int numItems, BCF2Type subtype ) {

        long trueSize = numItems == 0 ?
                        (Long)BCF2ParsingTable.parseNextValue(data, BCF2Constants.ATOMIC_INTEGRAL_TYPES).getValue() :
                        numItems;

        if ( trueSize == 0 ) {
            return new BCF2Value(this, new ArrayList<BCF2Value>(), true);
        }

        ArrayList<BCF2Value> vector = new ArrayList<BCF2Value>((int)trueSize);

        for ( long i = 0; i < trueSize; i++ ) {
            vector.add(subtype.decode(data, numItems));
        }
        return new BCF2Value(this, vector, false);
    }
}