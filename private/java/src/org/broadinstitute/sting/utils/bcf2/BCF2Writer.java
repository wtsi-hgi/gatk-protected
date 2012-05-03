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

import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderVersion;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.*;

public class BCF2Writer {

    private DataOutputStream bcf2OutputStream;  // TODO: needs to be little-endian
                                                // Note: do not flush until completely done writing, to avoid
                                                // issues with eventual BGZF support
    private boolean headerWritten;

    public BCF2Writer( File destination ) {
        this(openBCF2File(destination));
    }

    public BCF2Writer( DataOutputStream out ) {
        bcf2OutputStream = out;
        headerWritten = false;
    }

    public void writeHeader( VCFHeader header ) {
        try {
            writeHeaderLine(VCFHeader.METADATA_INDICATOR, BCF2Constants.VERSION_LINE);

            for ( VCFHeaderLine line : header.getMetaData() ) {
                if ( VCFHeaderVersion.isFormatString(line.getKey()) ) // Skip the version line that the VCFHeader class inserts
                    continue;

                writeHeaderLine(VCFHeader.METADATA_INDICATOR, line.toString());
            }

            writeHeaderLine(VCFHeader.HEADER_INDICATOR, constructHeaderFieldLayoutLine(header));
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile("Error writing BCF2 header", e);
        }

        headerWritten = true;
    }

    public void add( VariantContext vc ) {
        if ( ! headerWritten ) {
            throw new ReviewedStingException("Cannot write BCF2 records before the BCF2 header has been written");
        }
        try {
            writeChrom(vc);
            writePos(vc);
            writeID(vc);
            writeAlleles(vc);
            writeQual(vc);
            writeFilter(vc);
            writeInfo(vc);
            writeGenotypes(vc);
        }
        catch ( IOException e ) {
            throw new UserException("Error writing record to BCF2 file: " + vc.toString());
        }
    }

    public void close() {
        try {
            bcf2OutputStream.flush();
            bcf2OutputStream.close();
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
        bcf2OutputStream.writeChars(prefix);
        bcf2OutputStream.writeChars(line);
        bcf2OutputStream.writeChars("\n");
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

    private void writeChrom( VariantContext vc ) throws IOException {

    }

    private void writePos( VariantContext vc ) throws IOException {

    }

    private void writeID( VariantContext vc ) throws IOException {

    }

    private void writeAlleles( VariantContext vc ) throws IOException {

    }

    private void writeQual( VariantContext vc ) throws IOException {

    }

    private void writeFilter( VariantContext vc ) throws IOException {

    }

    private void writeInfo( VariantContext vc ) throws IOException {

    }

    private void writeGenotypes( VariantContext vc ) throws IOException {

    }

    // TODO: refactor these methods into the BCF2Type class or similar

    private void writeTypeDescriptor( BCF2Type type, long size ) throws IOException {
        int typeDescriptor = (type.isAtomic() ? 0 : 0x80) |
                             (type.isAtomic() || size > 7 ? 0 : ((byte)size & 0xFF) << 4) |
                             type.getTypeID();
        bcf2OutputStream.write(typeDescriptor);
    }
}
