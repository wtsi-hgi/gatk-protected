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

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.*;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BCF2Writer {

    public static final int WRITE_BUFFER_INITIAL_SIZE = 16384;

    private DataOutputStream bcf2File;      // Note: do not flush until completely done writing, to avoid
                                            // issues with eventual BGZF support

    private ByteArrayOutputStream recordBuffer = new ByteArrayOutputStream(WRITE_BUFFER_INITIAL_SIZE);
    private DataOutputStream record = new DataOutputStream(recordBuffer);   // TODO: DataOutputStream hardcodes big-endian encoding --
                                                                            // needs to be little-endian to be in spec

    private VCFHeader header;
    private SAMSequenceDictionary sequenceDictionary;
    private Map<String, Integer> stringDictionary = new HashMap<String, Integer>();

    public BCF2Writer( File destination, VCFHeader header, SAMSequenceDictionary sequenceDictionary ) {
        this(openBCF2File(destination), header, sequenceDictionary);
    }

    public BCF2Writer( DataOutputStream out, VCFHeader header, SAMSequenceDictionary sequenceDictionary ) {
        this.bcf2File = out;
        this.header = header;
        this.sequenceDictionary = sequenceDictionary;
    }

    public void writeHeader() {
        try {
            writeHeaderLine(VCFHeader.METADATA_INDICATOR, BCF2Constants.VERSION_LINE);

            for ( VCFHeaderLine line : header.getMetaData() ) {
                if ( VCFHeaderVersion.isFormatString(line.getKey()) ||     // Skip the version line that the VCFHeader class inserts
                     line.getKey().equals(VCFHeader.CONTIG_KEY)) {         // And any existing contig definitions (we'll add our own)
                    continue;
                }

                writeHeaderLine(VCFHeader.METADATA_INDICATOR, line.toString());
            }

            for ( SAMSequenceRecord contigEntry : sequenceDictionary.getSequences() ) {
                VCFSimpleHeaderLine contigHeaderLine = new VCFSimpleHeaderLine(VCFHeader.CONTIG_KEY, contigEntry.getSequenceName(), "");
                writeHeaderLine(VCFHeader.METADATA_INDICATOR, contigHeaderLine.toString());
            }

            writeHeaderLine(VCFHeader.HEADER_INDICATOR, constructHeaderFieldLayoutLine(header));
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile("Error writing BCF2 header", e);
        }
    }

    public void add( VariantContext vc ) {
        try {
            buildChrom(vc);
            buildPos(vc);
            buildID(vc);
            buildAlleles(vc);
            buildQual(vc);
            buildFilter(vc);
            buildInfo(vc);
            //buildGenotypes(vc);

            writeRecordToOutputFile();
            recordBuffer.reset();     // reuse buffer for next record
        }
        catch ( IOException e ) {
            throw new UserException("Error writing record to BCF2 file: " + vc.toString());
        }
    }

    public void close() {
        try {
            bcf2File.flush();
            bcf2File.close();
        }
        catch ( IOException e ) {
            throw new UserException("Failed to close BCF2 file");
        }
    }

    private static DataOutputStream openBCF2File( File bcf2File ) {
        try {
            return new DataOutputStream(new FileOutputStream(bcf2File));
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile(bcf2File, "Failed to open BCF2 file for writing", e);
        }
    }

    private void writeHeaderLine( String prefix, String line ) throws IOException {
        bcf2File.write(prefix.getBytes(BCF2Constants.BCF2_TEXT_CHARSET));
        bcf2File.write(line.getBytes(BCF2Constants.BCF2_TEXT_CHARSET));
        bcf2File.write("\n".getBytes(BCF2Constants.BCF2_TEXT_CHARSET));
    }

    private String constructHeaderFieldLayoutLine( VCFHeader header ) {
        StringBuilder fieldLayoutLine = new StringBuilder();

        for ( VCFHeader.HEADER_FIELDS field : header.getHeaderFields() ) {
            fieldLayoutLine.append(field.toString());
            fieldLayoutLine.append(VCFConstants.FIELD_SEPARATOR);
        }

        if ( header.hasGenotypingData() ) {
            fieldLayoutLine.append("FORMAT");

            for ( String sample : header.getGenotypeSamples() ) {
                fieldLayoutLine.append(VCFConstants.FIELD_SEPARATOR);
                fieldLayoutLine.append(sample);
            }
        }

        return fieldLayoutLine.toString();
    }

    private void writeRecordToOutputFile() throws IOException {
        bcf2File.writeInt(recordBuffer.size());
        bcf2File.write(recordBuffer.toByteArray());
    }

    private void buildChrom( VariantContext vc ) throws IOException {
        Integer contigIndex = sequenceDictionary.getSequenceIndex(vc.getChr());

        if ( contigIndex == -1 ) {
            throw new UserException(String.format("Contig %s not found in sequence dictionary from reference",
                                                   vc.getChr()));
        }

        writeOptimallySizedInteger(contigIndex);
    }

    private void buildPos( VariantContext vc ) throws IOException {
        writeOptimallySizedInteger(vc.getStart());
    }

    private void buildID( VariantContext vc ) throws IOException {
        writeLiteralString(vc.getID());
    }

    private void buildAlleles( VariantContext vc ) throws IOException {
        writeString(vc.getReference().getDisplayString());

        List<Allele> altAlleles = vc.getAlternateAlleles();

        if ( altAlleles.size() == 0 ) {
            writeString("");
        }
        else if ( altAlleles.size() == 1 ) {
            writeString(altAlleles.get(0).getDisplayString());
        }
        else {
            List<String> altAlleleStrings = new ArrayList<String>();
            for ( Allele altAllele : vc.getAlternateAlleles() ) {
                altAlleleStrings.add(altAllele.getDisplayString());
            }
            writeStringVector(altAlleleStrings);
        }
    }

    private void buildQual( VariantContext vc ) throws IOException {
        if ( ! vc.hasLog10PError() ) {
            writeMissingFloat();
        }
        else {
            writeFloat((float)vc.getPhredScaledQual());
        }
    }

    private void buildFilter( VariantContext vc ) throws IOException {
        List<String> filters = new ArrayList<String>(vc.getFilters());

        if ( filters.size() == 0 ) {
            writeString("");
        }
        else if ( filters.size() == 1 ) {
            writeString(filters.get(0));
        }
        else {
            writeStringVector(filters);
        }
    }

    private void buildInfo( VariantContext vc ) throws IOException {
        int numInfoFields = vc.getAttributes().size();
        writeOptimallySizedInteger(numInfoFields);

        for ( Map.Entry<String, Object> infoFieldEntry : vc.getAttributes().entrySet() ) {
            writeString(infoFieldEntry.getKey());
            writeString(infoFieldEntry.getValue().toString()); // TODO: get the actual types from the header
        }
    }

    private void buildGenotypes( VariantContext vc ) throws IOException {

    }

    // TODO: refactor these methods into the BCF2Type class or similar

    private void writeTypeDescriptor( BCF2Type type, long size ) throws IOException {
        int typeDescriptor = (type.isAtomic() ? 0 : 0x80) |
                             (type.isAtomic() || size > 7 ? 0 : ((byte)size & 0xFF) << 4) |
                             type.getTypeID();

        record.writeByte(typeDescriptor);
        if ( size > 7 ) {
            writeOptimallySizedInteger(size);
        }
    }

    private void writeOptimallySizedInteger( long value ) throws IOException {
        if ( value <= Byte.MAX_VALUE ) {
            writeTypeDescriptor(BCF2Type.INT8, 0);
            record.writeByte((int)value);
        }
        else if ( value <= Short.MAX_VALUE ) {
            writeTypeDescriptor(BCF2Type.INT16, 0);
            record.writeShort((int)value);
        }
        else if ( value <= Integer.MAX_VALUE ) {
            writeTypeDescriptor(BCF2Type.INT32, 0);
            record.writeInt((int) value);
        }
        else {
            writeTypeDescriptor(BCF2Type.INT64, 0);
            record.writeLong(value);
        }
    }

    private void writeString( String value ) throws IOException {
        if ( stringDictionary.containsKey(value) ) {
            writeOptimallySizedDictionaryString(stringDictionary.get(value));
        }
        else {
            writeLiteralString(value);
        }
    }

    private void writeLiteralString( String value ) throws IOException {
        writeTypeDescriptor(BCF2Type.STRING_LITERAL, 0);
        record.write(value.getBytes(BCF2Constants.BCF2_TEXT_CHARSET));
        record.writeByte(0);
    }

    private void writeOptimallySizedDictionaryString( int dictionaryOffset ) throws IOException {
        if ( dictionaryOffset <= Byte.MAX_VALUE ) {
            writeTypeDescriptor(BCF2Type.STRING_REF8, 0);
            record.writeByte(dictionaryOffset);
        }
        else if ( dictionaryOffset <= Short.MAX_VALUE ) {
            writeTypeDescriptor(BCF2Type.STRING_REF16, 0);
            record.writeShort(dictionaryOffset);
        }
        else {
            writeTypeDescriptor(BCF2Type.STRING_REF32, 0);
            record.writeInt(dictionaryOffset);
        }
    }

    private void writeDouble( double value ) throws IOException {
        writeTypeDescriptor(BCF2Type.DOUBLE, 0);
        record.writeDouble(value);
    }

    private void writeMissingDouble() throws IOException {
        writeTypeDescriptor(BCF2Type.DOUBLE, 0);
        ByteBuffer buf = ByteBuffer.allocate(8);
        buf.order(BCF2Constants.BCF2_BYTE_ORDER);
        record.write(buf.putLong(BCF2Constants.DOUBLE_MISSING_VALUE).array());
    }

    private void writeFloat( float value ) throws IOException {
        writeTypeDescriptor(BCF2Type.FLOAT, 0);
        record.writeFloat(value);
    }

    private void writeMissingFloat() throws IOException {
        writeTypeDescriptor(BCF2Type.FLOAT, 0);
        ByteBuffer buf = ByteBuffer.allocate(4);
        buf.order(BCF2Constants.BCF2_BYTE_ORDER);
        record.write(buf.putInt(BCF2Constants.FLOAT_MISSING_VALUE).array());
    }

    // TODO: need a generalized vector writing method
    private void writeStringVector( List<String> values ) throws IOException {
        writeTypeDescriptor(BCF2Type.VECTOR_STRING_LITERAL, values.size());

        for ( String value : values ) {
            record.write(value.getBytes(BCF2Constants.BCF2_TEXT_CHARSET));
            record.writeByte(0);
        }
    }
}
