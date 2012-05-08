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

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.Map;

/**
 * Simple BCF2 encoder
 *
 * @author Your Name
 * @since Date created
 */
public class BCFEncoder {
    public static final int WRITE_BUFFER_INITIAL_SIZE = 16384;
    private ByteArrayOutputStream encodeStream = new ByteArrayOutputStream(WRITE_BUFFER_INITIAL_SIZE);
    //private DataOutputStream record = new DataOutputStream(recordBuffer);   // TODO: DataOutputStream hardcodes big-endian encoding --
    // needs to be little-endian to be in spec

    private final Map<String, Integer> stringDictionary;

    public BCFEncoder() {
        // todo -- real constructor
        stringDictionary = Collections.emptyMap();
    }

    public int getRecordSizeInBytes() {
        return encodeStream.size();
    }

    public byte[] getRecordBytes() {
        byte[] bytes = encodeStream.toByteArray();
        encodeStream.reset();
        return bytes;
    }

    public final void encodeInt(final int value) throws IOException {
        final BCFType type = determineIntegerType(value);
        encodeSingleton(value, type);
    }

    public final void encodeMissing(final BCFType type) throws IOException {
        encodeVector(Collections.emptyList(), type);
    }

    public final void encodeMissingValues(final int size, final BCFType type) throws IOException {
        for ( int i = 0; i < size; i++ )
            encodeValue(type.getMissingValue(), type);
    }

    // todo -- should be specialized for each object type for efficiency
    public final void encodeSingleton(final Object v, final BCFType type) throws IOException {
        encodeVector(Collections.singleton(v), type);
    }

    public final <T extends Object> void encodeVector(final Collection<T> v, final BCFType type) throws IOException {
        encodeType(v.size(), type);
        encodeValues(v, type);
    }

    public final <T extends Object> void encodeValues(final Collection<T> v, final BCFType type) throws IOException {
        for ( final T v1 : v ) {
            encodeValue(v1, type);
        }
    }

    public final <T extends Object> void encodeValue(final T value, final BCFType type) throws IOException {
        switch (type) {
            case RESERVED_0: throw new ReviewedStingException("Type is 0");
            case INT8:
            case INT16:
            case INT32:  encodePrimitive((Integer)value, type); break;
            //case INT64:
            case FLOAT:  encodePrimitive(Float.floatToIntBits((Float)value), type); break;
            //case DOUBLE: encodePrimitive((long)(double)(Double)value, type); break;
            case FLAG: throw new ReviewedStingException("Flag encoding not supported");
            //case CHAR: // char is just an array of string
            case STRING_LITERAL: encodeLiteralString((String)value); break;
            case STRING_REF8:
            case STRING_REF16: encodePrimitive(stringDictionary.get(value), type); break;
            //case STRING_REF32:
            case COMPACT_GENOTYPE: encodeCompactGenotype((Byte)value); break;
            case RESERVED_14:
            case RESERVED_15: throw new ReviewedStingException("Type reserved " + type);
            default: throw new ReviewedStingException("Impossible state");
        }
    }

    public final BCFType determineStringType(final String value) {
        if ( stringDictionary.containsKey(value) ) {
            final int offset = stringDictionary.get(value);
            return determizeBestTypeBySize(offset, TypeDescriptor.DICTIONARY_TYPES_BY_SIZE);
        } else {
            return BCFType.STRING_LITERAL;
        }
    }

    public final BCFType determineIntegerType(final int value) {
        return determizeBestTypeBySize(value, TypeDescriptor.INTEGER_TYPES_BY_SIZE);
    }

    public final void encodeLiteralString(final String string) throws IOException {
        encodeStream.write(string.getBytes(BCF2Constants.BCF2_TEXT_CHARSET));
        encodeStream.write(0x00);
    }

    private final BCFType determizeBestTypeBySize(final int value, final BCFType[] potentialTypesInSizeOrder) {
        for ( final BCFType potentialType : potentialTypesInSizeOrder ) {
            if ( potentialType.withinRange(value) )
                return potentialType;
        }
        // TODO -- throw error
        return null;
    }

    public final void encodePrimitive(final int value, final BCFType type) throws IOException {
        for ( int i = type.getSizeInBytes() - 1; i >= 0; i-- ) {
            final int shift = i * 8;
            int mask = 0xFF << shift;
            byte byteValue = (byte)((mask & value) >> shift);
            encodeStream.write(byteValue);
        }
    }

    public final void encodeCompactGenotype(final byte value) throws IOException {
        encodeStream.write(value);
    }

    private final void encodeType(final int size, final BCFType type) throws IOException {
        final byte typeByte = TypeDescriptor.encodeTypeDescriptor((int)size, type);
        encodeStream.write(typeByte);
        if ( TypeDescriptor.willOverflow(size) )
            encodeSingleton(size, determineIntegerType(size));
    }

    public final void startGenotypeField(final String key, final int size, final BCFType valueType) throws IOException {
        encodeSingleton(key, BCFType.STRING_LITERAL);
        encodeType(size, valueType);
    }
}
