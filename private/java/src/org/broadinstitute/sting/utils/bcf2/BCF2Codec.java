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

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;

public class BCF2Codec implements FeatureCodec<VariantContext> {
    private VCFHeader header = null;
    private final ArrayList<String> contigNames = new ArrayList<String>();
    private final ArrayList<String> dictionary = new ArrayList<String>();
    InputStream stream;

    public Feature decodeLoc( final PositionalBufferedStream inputStream ) {
        return decode(inputStream);  // TODO: a less expensive version of decodeLoc() that doesn't use VariantContext
    }

    public VariantContext decode( final PositionalBufferedStream inputStream ) {
        prepareByteStream(inputStream);
        final VariantContextBuilder builder = new VariantContextBuilder();

        decodeChrom(builder);
        decodePos(builder);
        decodeID(builder);
        decodeAlleles(builder);
        decodeQual(builder);
        decodeFilter(builder);
        decodeInfo(builder);
        //decodeGenotypes(iterator, builder);

        return builder.make();
    }

    public Class<VariantContext> getFeatureType() {
        return VariantContext.class;
    }

    public FeatureCodecHeader readHeader( final PositionalBufferedStream inputStream ) {
        AsciiLineReader headerReader = new AsciiLineReader(inputStream);
        String headerLine;
        List<String> headerLines = new ArrayList<String>();
        boolean foundHeaderEnd = false;

        try {
            while ( ! foundHeaderEnd && (headerLine = headerReader.readLine()) != null) {
                if ( headerLine.startsWith(VCFHeader.METADATA_INDICATOR) ) {
                    headerLines.add(headerLine);
                }
                else if ( headerLine.startsWith(VCFHeader.HEADER_INDICATOR) ) {
                    headerLines.add(headerLine);
                    foundHeaderEnd = true;
                }
                else {
                    throw new UserException.MalformedBCF2("Reached end of header without encountering a field layout line");
                }
            }
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("I/O error while reading BCF2 header");
        }

        if ( ! foundHeaderEnd ) {
            throw new UserException.MalformedBCF2("Reached end of header without encountering a field layout line");
        }

        VCFHeader parsedHeader = createHeader(headerLines, VCFHeaderVersion.VCF4_1);
        this.header = parsedHeader;
        return new FeatureCodecHeader(parsedHeader, inputStream.getPosition());  // TODO -- position doesn't work
    }

    private VCFHeader createHeader( List<String> headerStrings, VCFHeaderVersion version ) {
        // Adapted from the AbstractVCFCodec class

        Set<VCFHeaderLine> metaData = new TreeSet<VCFHeaderLine>();
        Set<String> sampleNames = new LinkedHashSet<String>();

        for ( String str : headerStrings ) {
            if ( ! str.startsWith(VCFHeader.METADATA_INDICATOR) ) {
                String[] strings = str.substring(1).split(VCFConstants.FIELD_SEPARATOR);
                if ( strings.length < VCFHeader.HEADER_FIELDS.values().length )
                    throw new UserException.MalformedBCF2("there are not enough columns present in the header line: " + str);

                int arrayIndex = 0;
                for ( VCFHeader.HEADER_FIELDS field : VCFHeader.HEADER_FIELDS.values() ) {
                    try {
                        if ( field != VCFHeader.HEADER_FIELDS.valueOf(strings[arrayIndex]) )
                            throw new UserException.MalformedBCF2("we were expecting column name '" + field + "' but we saw '" + strings[arrayIndex] + "'");
                    } catch ( IllegalArgumentException e ) {
                        throw new UserException.MalformedBCF2("unknown column name '" + strings[arrayIndex] + "'; it does not match a legal column header name.");
                    }
                    arrayIndex++;
                }

                boolean sawFormatTag = false;
                if ( arrayIndex < strings.length ) {
                    if ( !strings[arrayIndex].equals("FORMAT") )
                        throw new UserException.MalformedBCF2("we were expecting column name 'FORMAT' but we saw '" + strings[arrayIndex] + "'");
                    sawFormatTag = true;
                    arrayIndex++;
                }

                while ( arrayIndex < strings.length )
                    sampleNames.add(strings[arrayIndex++]);

                if ( sawFormatTag && sampleNames.size() == 0 )
                    throw new UserException.MalformedBCF2("The FORMAT field was provided but there is no genotype/sample data");

            } else {
                if ( str.startsWith(VCFConstants.INFO_HEADER_START) ) {
                    final VCFInfoHeaderLine info = new VCFInfoHeaderLine(str.substring(7),version);
                    metaData.add(info);
                } else if ( str.startsWith(VCFConstants.FILTER_HEADER_START) ) {
                    final VCFFilterHeaderLine filter = new VCFFilterHeaderLine(str.substring(9), version);
                    metaData.add(filter);
                } else if ( str.startsWith(VCFConstants.FORMAT_HEADER_START) ) {
                    final VCFFormatHeaderLine format = new VCFFormatHeaderLine(str.substring(9), version);
                    metaData.add(format);
                } else if ( str.startsWith(VCFConstants.CONTIG_HEADER_START) ) {
                    final VCFSimpleHeaderLine contig = new VCFSimpleHeaderLine(str.substring(9), version, VCFConstants.CONTIG_HEADER_START.substring(2), null);
                    contigNames.add(contig.getID());
                    metaData.add(contig);
                } else if ( str.startsWith(VCFConstants.ALT_HEADER_START) ) {
                    final VCFSimpleHeaderLine alt = new VCFSimpleHeaderLine(str.substring(6), version, VCFConstants.ALT_HEADER_START.substring(2), Arrays.asList("ID", "Description"));
                    metaData.add(alt);
                } else {
                    int equals = str.indexOf("=");
                    if ( equals != -1 )
                        metaData.add(new VCFHeaderLine(str.substring(2, equals), str.substring(equals+1)));
                }
            }
        }

        header = new VCFHeader(metaData, sampleNames);
        header.buildVCFReaderMaps(new ArrayList<String>(sampleNames));
        return header;
    }

    private final void prepareByteStream(final InputStream inputStream) {
        this.stream = inputStream;
        final int sizeInBytes = decodeInt(4);
        final byte[] record = new byte[sizeInBytes];
        try {
            final int bytesRead = stream.read(record);
            if ( bytesRead < sizeInBytes ) {
                throw new UserException.MalformedBCF2(String.format("Failed to read next complete record: %s",
                        bytesRead == -1 ?
                                "premature end of input stream" :
                                String.format("expected %d bytes but read only %d", sizeInBytes, bytesRead)));
            }
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("I/O error while reading BCF2 file", e);
        }

        this.stream = new ByteArrayInputStream(record);
    }

    private void decodeChrom( final VariantContextBuilder builder ) {
        int contigOffset = (Integer)decodeRequiredTypedValue("CHROM");
        final String contig = lookupContigName(contigOffset);
        builder.chr(contig);
    }

    private void decodePos( final VariantContextBuilder builder ) {
        final int pos = (Integer)decodeRequiredTypedValue("POS");
        builder.start((long)pos);
        builder.stop((long)pos);
    }

    private void decodeID( final VariantContextBuilder builder ) {
        final String id = (String)decodeTypedValue();

        if ( id == null ) {
            builder.noID();
        }
        else {
            builder.id(id);
        }
    }

    private void decodeAlleles( final VariantContextBuilder builder ) {
        final String ref = (String)decodeTypedValue();
        final Object altObject = decodeTypedValue(); // alt can be either atomic or vectors

        // TODO -- probably need inline decoder for efficiency here (no sense in going bytes -> string -> vector -> bytes
        final List<Allele> alleles = new ArrayList<Allele>(2);
        alleles.add(Allele.create(ref, true));

        for ( final String alt : asStrings(altObject)) {
            alleles.add(Allele.create(alt));
        }

        // TODO: call into VCF parser to validate alleles like VCFCodec does
        builder.alleles(alleles);
    }

    private void decodeQual( final VariantContextBuilder builder ) {
        final Object qual = decodeTypedValue();

        if ( qual != null ) {
            builder.log10PError((Double)qual); // todo -- needs to know actual type unfortunately
        }
    }

    private void decodeFilter( final VariantContextBuilder builder ) {
        final Object filters = decodeTypedValue();

        if ( filters == null ) {
            builder.unfiltered();
        }
        else {
            builder.filters(new LinkedHashSet<String>(asStrings(filters)));
        }
    }

    private void decodeInfo( final VariantContextBuilder builder ) {
        final int numInfoFields = (Integer)decodeRequiredTypedValue("Number of info fields");
        final Map<String, Object> infoFieldEntries = new HashMap<String, Object>(numInfoFields);

        for ( int i = 0; i < numInfoFields; i++ ) {
            final String key = (String)decodeTypedValue();
            final Object value = decodeTypedValue();
            infoFieldEntries.put(key, value);
            // TODO: maybe? type-check against the header
        }

        builder.attributes(infoFieldEntries);
    }

    private void decodeGenotypes( final VariantContextBuilder builder ) {
        // TODO
    }

    private final Object decodeRequiredTypedValue(final String field) {
        final Object result = decodeTypedValue();
        if ( result == null ) {
            throw new UserException.MalformedBCF2("The value for the required field " + field + " is missing");
        } else {
            return result;
        }
    }

    private final Object decodeTypedValue() {
        byte typeDescriptor = readByte(stream);
        return decodeTypedValue(typeDescriptor);
    }

    private final Object decodeTypedValue(final byte typeDescriptor) {
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

    private final Object decodeValue(BCFType type) {
        switch (type) {
            case RESERVED_0: throw new ReviewedStingException("Type is 0");
            case INT8:   return decodeInt(type.getSizeInBytes());
            case INT16:  return decodeInt(type.getSizeInBytes());
            case INT32:  return decodeInt(type.getSizeInBytes());
            case INT64:  return (long)decodeInt(type.getSizeInBytes());
            case FLOAT:  return decodeDouble(type.getSizeInBytes());
            case DOUBLE: return decodeDouble(type.getSizeInBytes());
            case CHAR: throw new ReviewedStingException("Char encoding not supported");
            case FLAG: throw new ReviewedStingException("Flag encoding not supported");
            case STRING_LITERAL: return decodeLiteralString();
            case STRING_REF8:
            case STRING_REF16:
            case STRING_REF32: return getDictionaryString((int)decodeInt(type.getSizeInBytes()));
            //case COMPACT_GENOTYPE: return decodeInt(type); // todo -- must decode further in outer loop
            case RESERVED_14:
            case RESERVED_15: throw new ReviewedStingException("Type reserved " + type);
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

    private String lookupContigName( long contigOffset ) {
        if ( contigOffset < contigNames.size() ) {
            return contigNames.get((int)contigOffset);
        }
        else {
            throw new UserException.MalformedBCF2(String.format("No contig at index %d present in the sequence dictionary from the BCF2 header (%s)",
                    (int)contigOffset, contigNames));
        }
    }

    private final int decodeInt(int bytesForEachInt) {
        int value = 0;
        for ( int i = bytesForEachInt - 1; i >= 0; i-- ) {
            final int b = readByte(stream) & 0xFF;
            final int shift = i * 8;
            value |= b << shift;
        }
        return value;
    }

    private final double decodeDouble(final int sizeInBytes) {
        // TODO -- handle missing
        if ( sizeInBytes == 4 )
            return (double)(float)(decodeInt(sizeInBytes) & 0xFFFFFFFF);
        else
            // TODO -- fixme
            return (double)(decodeInt(sizeInBytes));
    }


    @Override
    public boolean canDecode( final String path ) {
        try {
            AsciiLineReader reader = new AsciiLineReader(new PositionalBufferedStream(new FileInputStream(path)));
            String firstLine = reader.readLine();
            if ( BCF2Constants.VERSION_LINE.equals(firstLine) ) {
                return true;
            }
        }
        catch ( IOException e ) {
            return false;
        }

        return false;
    }

    private final Collection<String> asStrings(final Object o) {
        return asCollection(String.class, o);
    }

    private final <T> Collection<T> asCollection(final Class<T> c, final Object o) {
        if ( o == null )
            return Collections.emptyList();
        else if ( o instanceof List ) {
            return (List<T>)o;
        } else {
            return (Set<T>)Collections.singleton(o);
        }
    }

    private final static byte readByte(final InputStream stream) {
        try {
            return (byte)(stream.read() & 0xFF);
        } catch ( IOException e ) {
            throw new ReviewedStingException("readByte failure", e);
        }
    }
}
