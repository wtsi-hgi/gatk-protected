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

package org.broadinstitute.sting.utils.collections;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

public class ExperimentalNestedIntegerArray<T> extends NestedIntegerArray<T> {

    public ExperimentalNestedIntegerArray(final int... dimensions) {
        super(dimensions);

        // Pre-allocate the second level of the tree to avoid ever having to lock the top level in order to
        // create second-level branches on demand
        for ( int i = 0; i < data.length; i++ ) {
            data[i] = new Object[dimensions[1]];
        }
    }

    @Override
    public void put(final T value, final int... keys) { // WARNING! value comes before the keys!
        if ( keys.length != numDimensions )
            throw new ReviewedStingException("Exactly " + numDimensions + " keys should be passed to this NestedIntegerArray but " + keys.length + " were provided");

        final int numNestedDimensions = numDimensions - 1;
        Object[] myData = data;
        for ( int i = 0; i < numNestedDimensions; i++ ) {
            if ( keys[i] >= dimensions[i] )
                throw new ReviewedStingException("Key " + keys[i] + " is too large for dimension " + i + " (max is " + (dimensions[i]-1) + ")");

            Object[] temp;

            // If we're at the top level, no need to lock, since we know the second level has already been
            // completely filled in
            if ( i == 0 ) {
                temp = (Object[])myData[keys[i]];  // myData[keys[i]] guaranteed to be non-null in first dimension
            }
            else {
                synchronized(myData) {  // yes, we really do want to synchronize on myData, even though it's local,
                                        // as multiple threads could be traversing the same part of the tree simultaneously
                    temp = (Object[])myData[keys[i]];
                    if ( temp == null ) {
                        temp = new Object[dimensions[i+1]];
                        myData[keys[i]] = temp;
                    }
                }
            }

            myData = temp;
        }

        myData[keys[numNestedDimensions]] = value;
    }
}
