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
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.util.*;

public class BCF2Codec implements FeatureCodec<VariantContext> {
    private VCFHeader header = null;
    private List<String> contigNames = new ArrayList<String>();

    public Feature decodeLoc( final PositionalBufferedStream inputStream ) {
        return decode(inputStream);  // TODO: a less expensive version of decodeLoc() that doesn't use VariantContext
    }

    public VariantContext decode( final PositionalBufferedStream inputStream ) {
        ByteBuffer record = loadNextRecord(inputStream);
        BCF2RecordIterator iterator = new BCF2RecordIterator(record);
        VariantContextBuilder builder = new VariantContextBuilder();

        decodeChrom(iterator, builder);
        decodePos(iterator, builder);
        decodeID(iterator, builder);
        decodeAlleles(iterator, builder);
        decodeQual(iterator, builder);
        decodeFilter(iterator, builder);
        decodeInfo(iterator, builder);
        //decodeGenotypes(iterator, builder);

        calculateStop(builder);

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
        return new FeatureCodecHeader(parsedHeader, inputStream.getPosition());
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

    private ByteBuffer loadNextRecord( InputStream inputStream ) {
        int nextRecordSize = getRecordSize(inputStream);
        return loadIntoByteBuffer(inputStream, nextRecordSize);
    }

    private ByteBuffer loadIntoByteBuffer( InputStream in, int bytesRequested ) {
        ByteBuffer buffer;
        byte[] bytesToLoad = new byte[bytesRequested];

        try {
            int bytesRead = in.read(bytesToLoad);
            if ( bytesRead < bytesRequested ) {
                throw new UserException.MalformedBCF2(String.format("Failed to read next complete record: %s",
                                                      bytesRead == -1 ?
                                                      "premature end of input stream" :
                                                      String.format("expected %d bytes but read only %d", bytesRequested, bytesRead)));
            }
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("I/O error while reading BCF2 file", e);
        }

        buffer = ByteBuffer.wrap(bytesToLoad);
        buffer.order(BCF2Constants.BCF2_BYTE_ORDER);
        return buffer;
    }

    private void decodeChrom( BCF2RecordIterator iterator, VariantContextBuilder builder ) {
        BCF2Value contigOffset = iterator.nextAtomicInteger();

        if ( contigOffset.isMissing() ) {
            throw new UserException.MalformedBCF2("Record missing the required CHROM field");
        }

        String contig = lookupContigName((Long)contigOffset.getValue());
        builder.chr(contig);
    }

    private void decodePos( BCF2RecordIterator iterator, VariantContextBuilder builder ) {
        BCF2Value pos = iterator.nextAtomicInteger();

        if ( pos.isMissing() ) {
            throw new UserException.MalformedBCF2("Record missing the required POS field");
        }

        builder.start((Long) pos.getValue());
    }

    private void decodeID( BCF2RecordIterator iterator, VariantContextBuilder builder ) {
        BCF2Value id = iterator.nextAtomicString();

        if ( id.isMissing() ) {
            builder.noID();
        }
        else {
            builder.id((String)id.getValue());
        }
    }

    private void decodeAlleles( BCF2RecordIterator iterator, VariantContextBuilder builder ) {
        BCF2Value ref = iterator.nextString();    // these can both be either atomic or vectors
        BCF2Value alt = iterator.nextString();

        List<Allele> alleles = new ArrayList<Allele>();

        if ( ref.getType().isAtomic() ) {
            alleles.add(Allele.create((String)ref.getValue(), true));
        }
        else {
            // TODO: Why does the spec allow multiple ref alleles?
            List<BCF2Value> rawRefAlleles = unpackVector(ref);
            for ( BCF2Value refAllele : rawRefAlleles ) {
                alleles.add(Allele.create((String)refAllele.getValue(), true));
            }
        }

        if ( alt.getType().isAtomic() ) {
            alleles.add(Allele.create((String)alt.getValue(), false));
        }
        else {
            List<BCF2Value> rawAltAlleles = unpackVector(alt);
            for ( BCF2Value altAllele : rawAltAlleles ) {
                alleles.add(Allele.create((String)altAllele.getValue(), false));
            }
        }

        // TODO: validate alleles like VCFCodec does

        builder.alleles(alleles);
    }

    private void decodeQual( BCF2RecordIterator iterator, VariantContextBuilder builder ) {
        BCF2Value qual = iterator.nextFloat();

        if ( ! qual.isMissing() ) {
            builder.log10PError((Float)qual.getValue());
        }
    }

    private void decodeFilter( BCF2RecordIterator iterator, VariantContextBuilder builder ) {
        BCF2Value filters = iterator.nextString();

        if ( filters.isMissing() ) {
            builder.unfiltered();
        }
        else {
            List<BCF2Value> filterList = unpackVector(filters);
            Set<String> filterSet = new LinkedHashSet<String>(filterList.size());

            for ( BCF2Value filter : filterList ) {
                filterSet.add((String)filter.getValue());
            }

            builder.filters(filterSet);
        }
    }

    private void decodeInfo( BCF2RecordIterator iterator, VariantContextBuilder builder ) {
        BCF2Value numInfoFields = iterator.nextAtomicInteger();
        long numInfoFieldsValue = (Long)numInfoFields.getValue();

        if ( numInfoFields.isMissing() ) {
            throw new UserException.MalformedBCF2("The value for the number of INFO fields cannot be missing");
        }

        Map<String, Object> infoFieldEntries = new HashMap<String, Object>((int)Math.min(numInfoFieldsValue, Integer.MAX_VALUE));

        for ( long i = 0; i < numInfoFieldsValue; i++ ) {
            BCF2Value key = iterator.nextAtomicString();
            BCF2Value value = iterator.next();  // TODO: use nextCustom() to type-check against the header

            infoFieldEntries.put((String)key.getValue(), value.getValue());  // TODO: potentially need to unpack vectors here
        }

        builder.attributes(infoFieldEntries);
    }

    private void decodeGenotypes( BCF2RecordIterator iterator, VariantContextBuilder builder ) {
        // TODO
    }

    private void calculateStop( VariantContextBuilder builder ) {
        // TODO
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

    private List<BCF2Value> unpackVector( BCF2Value vector ) {
        return (List<BCF2Value>)vector.getValue();
    }

    private int getRecordSize( InputStream in ) {
        // TODO: Terrible hack -- need to be able to parse the byte stream BEFORE the record is loaded into a buffer.
        //       For now, assume the record size is a 32-bit int with no type descriptor.

        ByteBuffer recordSizeBuffer = loadIntoByteBuffer(in, 4);
        recordSizeBuffer.order(BCF2Constants.BCF2_BYTE_ORDER);
        return recordSizeBuffer.getInt();
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
}
