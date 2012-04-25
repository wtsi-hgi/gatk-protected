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
import java.util.EnumSet;
import java.util.Iterator;

public class BCF2RecordIterator implements Iterator<BCF2Value>, Iterable<BCF2Value> {
    private ByteBuffer record;

    public BCF2RecordIterator( ByteBuffer record ) {
        this.record = record;
    }

    public boolean hasNext() {
        return record.hasRemaining();
    }

    // TODO: have these return strongly-typed values when possible to cut down on casting

    public BCF2Value next() {
        return readNextValue(BCF2Constants.ALL_TYPES);
    }

    public BCF2Value nextAtomicInteger() {
        return readNextValue(BCF2Constants.ATOMIC_INTEGRAL_TYPES);
    }

    public BCF2Value nextAtomicString() {
        return readNextValue(BCF2Constants.ATOMIC_STRING_TYPES);
    }

    public BCF2Value nextVectorString() {
        return readNextValue(BCF2Constants.VECTOR_STRING_TYPES);
    }

    public BCF2Value nextString() {
        return readNextValue(BCF2Constants.ALL_STRING_TYPES);
    }

    public BCF2Value nextFloat() {
        return readNextValue(EnumSet.of(BCF2Type.FLOAT));
    }

    public BCF2Value nextDouble() {
        return readNextValue(EnumSet.of(BCF2Type.DOUBLE));
    }

    public BCF2Value nextCustom( EnumSet<BCF2Type> expectedTypes ) {
        return readNextValue(expectedTypes);
    }

    private BCF2Value readNextValue( EnumSet<BCF2Type> expectedTypes ) {
        return BCF2ParsingTable.parseNextValue(record, expectedTypes);
    }

    public Iterator<BCF2Value> iterator() {
        return this;
    }

    public void remove() {
        throw new UnsupportedOperationException("Cannot remove values using a BCF2RecordIterator");
    }
}
