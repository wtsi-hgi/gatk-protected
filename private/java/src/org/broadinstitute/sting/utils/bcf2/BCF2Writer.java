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
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BCF2Writer {
    private DataOutputStream bcf2File;      // Note: do not flush until completely done writing, to avoid
                                            // issues with eventual BGZF support
    private VCFHeader header;
    private SAMSequenceDictionary sequenceDictionary;
    private Map<String, Integer> stringDictionary = new HashMap<String, Integer>();
    private final BCFEncoder encoder;

    public BCF2Writer( File destination, VCFHeader header, SAMSequenceDictionary sequenceDictionary ) {
        this(openBCF2File(destination), header, sequenceDictionary);
    }

    public BCF2Writer( DataOutputStream out, VCFHeader header, SAMSequenceDictionary sequenceDictionary ) {
        this.bcf2File = out;
        this.header = header;
        this.sequenceDictionary = sequenceDictionary;
        encoder = new BCFEncoder();
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
        bcf2File.writeInt(encoder.getRecordSizeInBytes());
        bcf2File.write(encoder.getRecordBytes());
    }

    private void buildChrom( VariantContext vc ) throws IOException {
        Integer contigIndex = sequenceDictionary.getSequenceIndex(vc.getChr());

        if ( contigIndex == -1 ) {
            throw new UserException(String.format("Contig %s not found in sequence dictionary from reference",
                                                   vc.getChr()));
        }

        encoder.encodeInt(contigIndex);
    }

    private void buildPos( VariantContext vc ) throws IOException {
        encoder.encodeInt(vc.getStart());
    }

    private void buildID( VariantContext vc ) throws IOException {
        encoder.encodeSingleton(vc.getID(), BCFType.STRING_LITERAL);
    }

    private void buildAlleles( VariantContext vc ) throws IOException {
        encoder.encodeSingleton(vc.getReference().getDisplayString(), BCFType.STRING_LITERAL);

        List<Allele> altAlleles = vc.getAlternateAlleles();

        if ( altAlleles.size() == 0 ) {
            encoder.encodeMissing(BCFType.STRING_LITERAL);
        } else {
            List<String> strings = new ArrayList<String>(altAlleles.size());
            for ( final Allele alt : altAlleles )
                strings.add(alt.getDisplayString());
            encoder.encodeVector(strings, BCFType.STRING_LITERAL);
        }
    }

    private void buildQual( VariantContext vc ) throws IOException {
        if ( ! vc.hasLog10PError() ) {
            encoder.encodeMissing(BCFType.DOUBLE);
        }
        else {
            encoder.encodeSingleton(vc.getPhredScaledQual(), BCFType.DOUBLE);
        }
    }

    private void buildFilter( VariantContext vc ) throws IOException {
        if ( vc.isFiltered() ) {
            encoder.encodeVector(vc.getFilters(), BCFType.STRING_LITERAL); // TODO -- need to determine best type
        } else {
            encoder.encodeMissing(BCFType.STRING_LITERAL);
        }
    }

    private void buildInfo( VariantContext vc ) throws IOException {
        int numInfoFields = vc.getAttributes().size();
        encoder.encodeInt(numInfoFields);

        for ( Map.Entry<String, Object> infoFieldEntry : vc.getAttributes().entrySet() ) {
            encoder.encodeSingleton(infoFieldEntry.getKey(), BCFType.STRING_LITERAL);
            // TODO -- determine type from header instead of using string
            encoder.encodeSingleton(infoFieldEntry.getValue().toString(), BCFType.STRING_LITERAL);
        }
    }

    private void buildGenotypes( VariantContext vc ) throws IOException {

    }
}
