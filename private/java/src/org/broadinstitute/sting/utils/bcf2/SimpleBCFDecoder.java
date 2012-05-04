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

import org.broad.tribble.BasicFeature;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

/**
 * Simple BCF decoder
 * @author Mark DePristo
 * @since 5/3/12
 */
public class SimpleBCFDecoder implements FeatureCodec<VariantContext> {
    ArrayList<String> dictionary;
    InputStream stream;

    @Override
    public Feature decodeLoc(final PositionalBufferedStream stream) throws IOException {
        prepareByteStream(stream);
        String chr = (String)decodeTypedValue();
        int start = (Integer)decodeTypedValue();
        int stop = start;
        return new BasicFeature(chr, start, stop);
    }

    private final void prepareByteStream(final InputStream inputStream) throws IOException {
        this.stream = inputStream;
        int sizeInBytes = (Integer)decodeTypedValue();
        byte[] record = new byte[sizeInBytes];
        stream.read(record);
        this.stream = new ByteArrayInputStream(record);
    }

    private final Object decodeTypedValue() {
        byte typeByte = readByte(stream);
        int size = typeByte >> 4;
        int type = typeByte & 0xF;

        if ( overflows(size) ) size = (Integer)decodeTypedValue();

        if ( size == 0 ) return null;
        else if ( size == 1 ) return decodeValue(type);
        else {
            ArrayList<Object> ints = new ArrayList<Object>(size);
            for ( int i = 0; i < size; i++ ) {
                ints.add(decodeValue(type));
            }
            return ints;
        }
    }

    private final Object decodeValue(int type) {
        int nBytes = -1;
        switch (type) {
            case 0: throw new ReviewedStingException("Type is 0");
            case 1: return decodeInt(1);
            case 2: return decodeInt(2);
            case 3: return decodeInt(4);
            case 4: return decodeInt(8);
            case 5:
            case 6: return decodeDouble(type == 5);
            case 7: throw new ReviewedStingException("Char encoding not supported");
            case 8: throw new ReviewedStingException("Flag encoding not supported");
            case 9: return decodeLiteralString();
            case 10:
            case 11:
            case 12: return getDictionaryString(decodeInt(type));
            case 13: return decodeInt(type); // todo -- must decode further in outer loop
            case 14:
            case 15: throw new ReviewedStingException("Type reserved " + type);
            default: throw new ReviewedStingException("Impossible state");
        }
    }

    private final String decodeLiteralString() {
        StringBuilder builder = new StringBuilder();
        while ( true ) {
            byte b = readByte(stream);
            if ( b == 0x00 ) break;
            builder.append((char)b);
        }
        return builder.toString();
    }

    private final String getDictionaryString(int offset) {
        return dictionary.get(offset);
    }

    private final int decodeInt(int bytesForEachInt) {
        int value = 0;
        for ( int i = bytesForEachInt - 1; i >= 0; i-- ) {
            final int b = readByte(stream);
            final int shift = i * 8;
            value |= b << shift;
        }
        return value;
    }

//    for ( int i = nBytes - 1; i >= 0; i-- ) {
//        final int shift = i * 8;
//        long mask = 0xFF << shift;
//        byte byteValue = (byte)((mask & value) >> shift);
//        encodeStream.write(byteValue);

    private final double decodeDouble(boolean asFloat) {
        try {
            if ( asFloat )
                return new DataInputStream(stream).readFloat();
            else
                return new DataInputStream(stream).readDouble();
        } catch (IOException e) {
            throw new ReviewedStingException("Couldn't read double", e);
        }
    }

    private final boolean overflows(int nElementsInType) {
        return nElementsInType == 15;
    }

    @Override
    public VariantContext decode(final PositionalBufferedStream stream) throws IOException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public FeatureCodecHeader readHeader(final PositionalBufferedStream stream) throws IOException {
//        VCFCodec vcfCodec = new VCFCodec();
        // TODO -- really need to know where the stream ends in VCFCodec readHeader
        return FeatureCodecHeader.EMPTY_HEADER;
//        return vcfCodec.readHeader(stream);
    }

    @Override
    public Class<VariantContext> getFeatureType() {
        return VariantContext.class;
    }

    @Override
    public boolean canDecode(final String path) {
        return false;
    }

    private final static byte readByte(final InputStream stream) {
        try {
            return (byte)(stream.read() & 0xFF);
        } catch ( IOException e ) {
            throw new ReviewedStingException("readByte failure", e);
        }
    }
}
