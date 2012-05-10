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

import org.apache.log4j.Logger;
import org.broad.tribble.FeatureCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

public class BCF2Decoder {
    final protected static Logger logger = Logger.getLogger(FeatureCodec.class);
    private final ArrayList<String> dictionary;

    byte[] recordBytes;
    ByteArrayInputStream recordStream;

    public BCF2Decoder(final ArrayList<String> dictionary) {
        if ( dictionary == null ) throw new ReviewedStingException("Dictionary cannot be null");
        this.dictionary = dictionary;
    }

    /**
     * Convenience constructor.  Builds a decoder and initialize it to start
     * reading the next record from stream
     *
     * @param stream
     */
    public BCF2Decoder(final ArrayList<String> dictionary, final InputStream stream) {
        this(dictionary, readRecordBytes(stream));
    }

    /**
     * Create a new decoder ready to read BCF2 data from the byte[] recordBytes
     *
     * @param recordBytes
     */
    public BCF2Decoder(final ArrayList<String> dictionary, final byte[] recordBytes) {
        this(dictionary);
        setRecordBytes(recordBytes);
    }

    /**
     * Reads the next record from input stream and prepare this decoder to decode values from it
     *
     * @param stream
     * @return
     */
    public void readNextBlock(final InputStream stream) {
        setRecordBytes(readRecordBytes(stream));
    }

    /**
     * Skips the next record from input stream, invalidating current block data
     *
     * @param stream
     * @return
     */
    public void skipNextBlock(final InputStream stream) {
        final int sizeInBytes = readBlockSize(stream);
        try {
            final int bytesRead = (int)stream.skip(sizeInBytes);
            validateReadBytes(bytesRead, sizeInBytes);
        } catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("I/O error while reading BCF2 file", e);
        }
        this.recordBytes = null;
        this.recordStream = null;
    }

    /**
     * Returns the byte[] for the block of data we are currently decoding
     * @return
     */
    public byte[] getRecordBytes() {
        return recordBytes;
    }

    /**
     * The size of the current block in bytes
     *
     * @return
     */
    public int getBlockSize() {
        return recordBytes.length;
    }

    /**
     * Use the recordBytes[] to read BCF2 records from now on
     *
     * @param recordBytes
     */
    public void setRecordBytes(final byte[] recordBytes) {
        this.recordBytes = recordBytes;
        this.recordStream = new ByteArrayInputStream(recordBytes);
    }

    public final Object decodeRequiredTypedValue(final String field) {
        final Object result = decodeTypedValue();
        if ( result == null ) {
            throw new UserException.MalformedBCF2("The value for the required field " + field + " is missing");
        } else {
            return result;
        }
    }

    public final Object decodeTypedValue() {
        final byte typeDescriptor = readTypeDescriptor();
        return decodeTypedValue(typeDescriptor);
    }

    public final byte readTypeDescriptor() {
        return readByte(recordStream);
    }

    public final Object decodeTypedValue(final byte typeDescriptor) {
        final int size = TypeDescriptor.sizeIsOverflow(typeDescriptor) ? (Integer)decodeTypedValue() : TypeDescriptor.decodeSize(typeDescriptor);
        final BCFType type = TypeDescriptor.decodeType(typeDescriptor);

        if ( size == 0 ) { return null; }
        else if ( size == 1 ) return decodeValue(type);
        else {
            ArrayList<Object> ints = new ArrayList<Object>(size);
            for ( int i = 0; i < size; i++ ) {
                ints.add(decodeValue(type));
            }
            return ints;
        }
    }

    public final Object decodeValue(final BCFType type) {
        switch (type) {
            case INT8:
            case INT16:
            case INT32:             return decodeInt(type.getSizeInBytes());
            case FLOAT:             return decodeFloatAsDouble();
            case FLAG:              throw new ReviewedStingException("Flag encoding not supported");
            case STRING_LITERAL:    return decodeLiteralString();
            case STRING_REF8:
            case STRING_REF16:      return getDictionaryString(decodeInt(type.getSizeInBytes()));
            case COMPACT_GENOTYPE:  return decodeInt(type.getSizeInBytes()); // todo -- must decode further in outer loop
            default:                throw new ReviewedStingException("BCF2 codec doesn't know how to decode type " + type );
        }
    }

    public final String decodeLiteralString() {
        // TODO -- fast path for missing value
        StringBuilder builder = new StringBuilder();
        while ( true ) {
            byte b = readByte(recordStream);
            if ( b == 0x00 ) break;
            builder.append((char)b);
        }
        return builder.toString();
    }

    private final String getDictionaryString(int offset) {
        return dictionary.get(offset);
    }

    // ----------------------------------------------------------------------
    //
    // Raw primative data types (ints and floats)
    //
    // ----------------------------------------------------------------------

    public final int decodeInt(int bytesForEachInt) {
        return readInt(bytesForEachInt, recordStream);
    }

    public final boolean rawFloatIsMissing(final int rawFloat) {
        return rawFloat == BCF2Constants.FLOAT_MISSING_VALUE;
    }

    public final int decodeRawFloat() {
        return decodeInt(BCFType.FLOAT.getSizeInBytes());
    }

    public final float rawFloatToFloat(final int rawFloat) {
        return Float.intBitsToFloat(rawFloat);
    }

    public final Double decodeFloatAsDouble() {
        final int i = decodeRawFloat();
        return rawFloatIsMissing(i) ? null : (double)rawFloatToFloat(i);
    }

    // ----------------------------------------------------------------------
    //
    // Utility functions
    //
    // ----------------------------------------------------------------------

    /**
     * Read the size of the next block from inputStream
     *
     * @param inputStream
     * @return
     */
    private final static int readBlockSize(final InputStream inputStream) {
        return readInt(4, inputStream);
    }

    /**
     *
     * @param inputStream
     * @return
     */
    private final static byte[] readRecordBytes(final InputStream inputStream) {
        final int sizeInBytes = readBlockSize(inputStream);
        final byte[] record = new byte[sizeInBytes];
        try {
            final int bytesRead = inputStream.read(record);
            validateReadBytes(bytesRead, sizeInBytes);
        } catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("I/O error while reading BCF2 file", e);
        }

        return record;
    }

    private final static void validateReadBytes(final int actuallyRead, final int expected) {
        if ( actuallyRead < expected ) {
            throw new UserException.MalformedBCF2(String.format("Failed to read next complete record: %s",
                    actuallyRead == -1 ?
                            "premature end of input stream" :
                            String.format("expected %d bytes but read only %d", expected, actuallyRead)));
        }
    }

    private final static byte readByte(final InputStream stream) {
        try {
            return (byte)(stream.read() & 0xFF);
        } catch ( IOException e ) {
            throw new ReviewedStingException("readByte failure", e);
        }
    }

    private final static int readInt(int bytesForEachInt, final InputStream stream) {
        int value = 0;
        for ( int i = bytesForEachInt - 1; i >= 0; i-- ) {
            final int b = readByte(stream) & 0xFF;
            final int shift = i * 8;
            value |= b << shift;
        }
        return value;
    }
}
