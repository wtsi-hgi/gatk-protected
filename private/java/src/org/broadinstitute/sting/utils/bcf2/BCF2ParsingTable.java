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

import org.broadinstitute.sting.utils.exceptions.UserException;

import java.nio.BufferUnderflowException;
import java.nio.ByteBuffer;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;

public class BCF2ParsingTable {
    private static final BCF2Type[] PARSING_TABLE = new BCF2Type[256];

    static {
        populateParsingTable();
    }

    public static BCF2Value parseNextValue( ByteBuffer record, EnumSet<BCF2Type> expectedTypes ) {
        BCF2Value val;

        try {
            byte typeDescriptor = record.get();
            int numItems = (typeDescriptor & BCF2Constants.TYPE_DESCRIPTOR_NUM_ELEMENTS_FIELD_MASK) >> 4;

            val = PARSING_TABLE[typeDescriptor & 0x00FF].decode(record, numItems);
        }
        catch ( BufferUnderflowException e ) {
            throw new UserException.MalformedBCF2("Premature end of BCF2 record");
        }

        if ( ! expectedTypes.contains(val.getType()) ) {
            throw new UserException.MalformedBCF2(String.format("Expected a value of type(s) %s but encountered a value of type %s",
                                                                expectedTypes, val.getType()));
        }

        return val;
    }

    private static void populateParsingTable() {
        Map<Integer, BCF2Type> atomicTypeIDMap = new HashMap<Integer, BCF2Type>();
        Map<Integer, BCF2Type> nonAtomicTypeIDMap = new HashMap<Integer, BCF2Type>();

        for ( BCF2Type type : BCF2Type.values() ) {
            if ( type.isAtomic() ) {
                atomicTypeIDMap.put(type.getTypeID(), type);
            }
            else {
                nonAtomicTypeIDMap.put(type.getTypeID(), type);
            }
        }

        for ( int typeDescriptor = 0; typeDescriptor < 256; typeDescriptor++ ) {
            int typeID = typeDescriptor & BCF2Constants.TYPE_DESCRIPTOR_TYPE_ID_MASK;
            boolean isAtomic = (typeDescriptor & BCF2Constants.TYPE_DESCRIPTOR_ATOMICITY_FIELD_MASK) == 0;

            BCF2Type correspondingType = isAtomic ? atomicTypeIDMap.get(typeID) : nonAtomicTypeIDMap.get(typeID);

            PARSING_TABLE[typeDescriptor] = correspondingType;
        }
    }
}
