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
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;

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

    public void add( final VariantContext initialVC ) {
        final VariantContext vc = initialVC.fullyDecode(header);

        try {
            buildChrom(vc);
            buildPos(vc);
            buildID(vc);
            buildAlleles(vc);
            buildQual(vc);
            buildFilter(vc);
            buildInfo(vc);
            buildSamplesData(vc);

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
            encoder.encodeMissing(BCFType.FLOAT);
        } else {
            encoder.encodeSingleton((float)vc.getPhredScaledQual(), BCFType.FLOAT);
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

    private void buildSamplesData(final VariantContext vc) throws IOException {
        // write size
        List<String> genotypeFields = StandardVCFWriter.calcVCFGenotypeKeys(vc);
        encoder.encodeInt(genotypeFields.size());
        for ( final String field : genotypeFields ) {
            if ( field.equals(VCFConstants.GENOTYPE_KEY) ) {
                addGenotypes(vc);
            } else if ( field.equals(VCFConstants.GENOTYPE_QUALITY_KEY) ) {
                addGQ(vc);
            } else if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY) ) {
                addGenotypeFilters(vc);
            } else {
                addGenericGenotypeField(vc, field);
            }
        }
    }

    private final int getNGenotypeFieldValues(final String field, final VariantContext vc) {
        final VCFCompoundHeaderLine metaData = VariantContext.getMetaDataForField(header, field);
        int nFields = metaData.getCount(vc.getAlternateAlleles().size());
        if ( nFields == -1 ) { // unbounded, need to look at values
            return computeMaxSizeOfGenotypeFieldFromValues(field, vc);
        } else {
            return nFields;
        }
    }

    private final int computeMaxSizeOfGenotypeFieldFromValues(final String field, final VariantContext vc) {
        int size = 1;
        final GenotypesContext gc = vc.getGenotypes();

        for ( final Genotype g : gc ) {
            final Object o = g.getAttribute(field);
            if ( o == null ) continue;
            if ( o instanceof List ) {
                // only do compute if first value is of type list
                final List values = (List)g.getAttribute(field);
                if ( values != null )
                    size = Math.max(size, values.size());
            } else {
                return 1;
            }
        }

        return size;
    }

    private final void addGenericGenotypeField(final VariantContext vc, final String field) throws IOException {
        final int numInFormatField = getNGenotypeFieldValues(field, vc);
        final BCFType type = getBCF2TypeFromHeader(field);

        encoder.startGenotypeField(field, numInFormatField, type);
        for ( final Genotype g : vc.getGenotypes() ) {
            if ( ! g.hasAttribute(field) ) {
                encoder.encodeMissingValues(numInFormatField, type);
            } else {
                final Object val = g.getAttribute(field);
                final Collection<Object> vals = numInFormatField == 1 ? Collections.singleton(val) : (Collection)val;
                encoder.encodeValues(convertToType(vals, type), type);
            }
        }
    }

    private final Collection<Object> convertToType(final Collection<Object> values, final BCFType type) {
        return values;
    }
//        if ( values == null ) return Collections.emptyList();
//        final List<Object> converted = new ArrayList<Object>(values.size());
//        for ( final Object val : values ) {
//            converted.add(type.convertIfNecessary(val));
//        }
//        return converted;
//    }

    private final BCFType getBCF2TypeFromHeader(final String field) {
        // TODO -- should take VC to determine best encoding
        final VCFCompoundHeaderLine metaData = VariantContext.getMetaDataForField(header, field);
        switch ( metaData.getType() ) {
            case Character: return BCFType.STRING_LITERAL;
            case Flag: return BCFType.FLAG;
            case String: return BCFType.STRING_LITERAL;
            case Integer: return BCFType.INT32;
            case Float: return BCFType.FLOAT;
            default: throw new ReviewedStingException("Unexpected type for field" + field);
        }
    }

    private final void addGenotypeFilters(final VariantContext vc) throws IOException {
        encoder.startGenotypeField(VCFConstants.GENOTYPE_FILTER_KEY, 1, BCFType.STRING_LITERAL);
        for ( final Genotype g : vc.getGenotypes() ) {
            if ( g.filtersWereApplied() && g.isFiltered() ) {
                encoder.encodeValue(ParsingUtils.join(";", ParsingUtils.sortList(g.getFilters())), BCFType.STRING_LITERAL);
            } else {
                encoder.encodeMissingValues(1, BCFType.STRING_LITERAL); // todo fixme
            }
        }
    }

    private final void addGQ(final VariantContext vc) throws IOException {
        encoder.startGenotypeField(VCFConstants.GENOTYPE_QUALITY_KEY, 1, BCFType.INT8);
        for ( final Genotype g : vc.getGenotypes() ) {
            if ( g.hasLog10PError() ) {
                final int GQ = (int)Math.round(Math.min(g.getPhredScaledQual(), VCFConstants.MAX_GENOTYPE_QUAL));
                if ( GQ > VCFConstants.MAX_GENOTYPE_QUAL ) throw new ReviewedStingException("Unexpectedly large GQ " + GQ + " at " + vc);
                encoder.encodeValue(GQ, BCFType.INT8);
            } else {
                encoder.encodeMissingValues(1, BCFType.INT8);
            }
        }
    }

    private final void addGenotypes(final VariantContext vc) throws IOException {
        final Map<Allele, String> alleleMap = StandardVCFWriter.buildAlleleMap(vc);

        final int requiredPloidy = 2;
        encoder.startGenotypeField(VCFConstants.GENOTYPE_KEY, requiredPloidy, BCFType.COMPACT_GENOTYPE);
        for ( final Genotype g : vc.getGenotypes() ) {
            if ( g.getPloidy() != requiredPloidy ) throw new ReviewedStingException("Cannot currently handle non-diploid calls!");
            final List<Byte> encoding = new ArrayList<Byte>(requiredPloidy);
            for ( final Allele a : g.getAlleles() ) {
                final int offset = a.isNoCall() ? 0 : (Byte.valueOf(alleleMap.get(a)) + 2);
                encoding.add((byte)(offset << 1 | (g.isPhased() ? 0x01 : 0x00)));
            }
            encoder.encodeValues(encoding, BCFType.COMPACT_GENOTYPE);
        }
    }
}
