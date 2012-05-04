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

import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Simple BCF2 encoder
 *
 * @author Your Name
 * @since Date created
 */
public class SimpleBCFEncoder {
    public SimpleBCFEncoder() {

    }

    public void encode(final VariantContext vc, final OutputStream out) throws IOException {
        ByteArrayOutputStream encodeStream = new ByteArrayOutputStream();
        encode(vc.getChr(), encodeStream);
        encode(vc.getStart(), encodeStream);
        encode(encodeStream.size(), out);
        out.write(encodeStream.toByteArray());
    }

    public final void encode(final String String, final OutputStream encodeStream) throws IOException {
        encodeType(1, BCF2Constants.STRING_LITERAL_TYPE_ID, encodeStream );
        encodeStream.write(String.getBytes());
        encodeStream.write(0x00);
    }

    public final void encode(final int value, final OutputStream encodeStream) throws IOException {
        encode((long)value, 4, encodeStream);
    }

    public final void encode(final long value, final OutputStream encodeStream) throws IOException {
        encode(value, 8, encodeStream);
    }

    public final void encode(final long value, int nBytes, final OutputStream encodeStream) throws IOException {
        encodeType(1, nBytes == 4 ? BCF2Constants.INT32_TYPE_ID : BCF2Constants.INT64_TYPE_ID, encodeStream);

        for ( int i = nBytes - 1; i >= 0; i-- ) {
            final int shift = i * 8;
            long mask = 0xFF << shift;
            byte byteValue = (byte)((mask & value) >> shift);
            encodeStream.write(byteValue);
        }
    }

    private final void encodeType(int size, int valueType, final OutputStream encodeStream) throws IOException {
        int originalSize = size;
        int encodeSize = Math.min(size, 15);
        byte typeByte = (byte)(encodeSize << 4 | (valueType & 0x0F));
        encodeStream.write(typeByte);
        if ( originalSize >= 15 )
            encode(size, encodeStream);
    }
}
