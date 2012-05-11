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
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.*;

public class BCF2Writer extends IndexingVCFWriter {
    private final static boolean doNotWriteGenotypes = false;
    private OutputStream outputStream;      // Note: do not flush until completely done writing, to avoid issues with eventual BGZF support
    private VCFHeader header;
    private Map<String, Integer> contigDictionary = new HashMap<String, Integer>();
    private Map<String, Integer> stringDictionary = new LinkedHashMap<String, Integer>();
    private BCFEncoder encoder; // initialized after the header arrives

    public BCF2Writer(final String name, final File location, final OutputStream output, final SAMSequenceDictionary refDict, final boolean enableOnTheFlyIndexing) {
        super(name, location, output, refDict, enableOnTheFlyIndexing);
        this.outputStream = getOutputStream();
    }

    // --------------------------------------------------------------------------------
    //
    // Interface functions
    //
    // --------------------------------------------------------------------------------

    @Override
    public void writeHeader(final VCFHeader header) {
        this.header = header;

        // create the config offsets map
        for ( final VCFContigHeaderLine contig : header.getContigLines())
            contigDictionary.put(contig.getID(), contig.getContigIndex());

        // set up the strings dictionary
        int offset = 0;
        stringDictionary.put(VCFConstants.PASSES_FILTERS_v4, offset++); // special case the special PASS field
        for ( VCFHeaderLine line : header.getMetaData() ) {
            if ( line instanceof VCFIDHeaderLine ) {
                VCFIDHeaderLine idLine = (VCFIDHeaderLine)line;
                stringDictionary.put(idLine.getID(), offset++);
            }
        }

        // add the dictionary ##dictionary=x,y,z line to the header
        final String dictionaryLineValue = Utils.join(BCF2Constants.DICTIONARY_LINE_ENTRY_SEPARATOR, stringDictionary.keySet());
        header.addMetaDataLine(new VCFHeaderLine(BCF2Constants.DICTIONARY_LINE_TAG, dictionaryLineValue));

        // write out the header
        StandardVCFWriter.writeHeader(header, new OutputStreamWriter(outputStream), doNotWriteGenotypes, BCF2Constants.VERSION_LINE, "BCF2 stream");

        // with the string dictionary in hand we can create the encoder
        encoder = new BCFEncoder(stringDictionary);
    }

    @Override
    public void add( final VariantContext initialVC ) {
        final VariantContext vc = initialVC.fullyDecode(header);
        super.add(vc); // allow on the fly indexing

        try {
            final byte[] infoBlock = buildSitesData(vc);
            final byte[] genotypesBlock = buildSamplesData(vc);

            // write the two blocks to disk
            writeBlock(infoBlock, genotypesBlock);
        }
        catch ( IOException e ) {
            throw new UserException("Error writing record to BCF2 file: " + vc.toString(), e);
        }
    }

    @Override
    public void close() {
        try {
            outputStream.flush();
            outputStream.close();
        }
        catch ( IOException e ) {
            throw new UserException("Failed to close BCF2 file");
        }
        super.close();
    }

    // --------------------------------------------------------------------------------
    //
    // implicit block
    //
    // The first four records of BCF are inline untype encoded data of:
    //
    // 4 byte integer chrom offset
    // 4 byte integer start
    // 4 byte integer ref length
    // 4 byte float qual
    //
    // --------------------------------------------------------------------------------
    private byte[] buildSitesData( VariantContext vc ) throws IOException {
        final int contigIndex = contigDictionary.get(vc.getChr());
        if ( contigIndex == -1 )
            throw new UserException(String.format("Contig %s not found in sequence dictionary from reference", vc.getChr()));

        // note use of encodeValue to not insert the typing byte
        encoder.encodeValue(contigIndex, BCFType.INT32);

        // pos
        encoder.encodeValue(vc.getStart(), BCFType.INT32);

        // ref length
        encoder.encodeValue(vc.getEnd() - vc.getStart() + 1, BCFType.INT32);

        // qual
        if ( vc.hasLog10PError() )
            encoder.encodeFloat((float)vc.getPhredScaledQual(), BCFType.FLOAT);
        else
            encoder.encodeMissingValue(BCFType.FLOAT);

        // info fields
        final int nAlleles = vc.getNAlleles();
        final int nInfo = vc.getAttributes().size();
        final int nGenotypeFormatFields = StandardVCFWriter.calcVCFGenotypeKeys(vc).size();
        final int nSamples = vc.getNSamples();

        encoder.encodePrimitive((nAlleles << 16) | (nInfo & 0x00FF), BCFType.INT32);
        encoder.encodePrimitive((nGenotypeFormatFields << 24) | (nSamples & 0x0FFF), BCFType.INT32);

        buildID(vc);
        buildAlleles(vc);
        buildFilter(vc);
        buildInfo(vc);

        return encoder.getRecordBytes();
    }

    private void buildID( VariantContext vc ) throws IOException {
        encoder.encodeString(vc.getID());
    }

    private void buildAlleles( VariantContext vc ) throws IOException {
        for ( final Allele allele : vc.getAlleles() ) {
            final String s = vc.getAlleleWithRefPadding(allele);
            encoder.encodeString(s);
        }
    }

    private void buildFilter( VariantContext vc ) throws IOException {
        if ( vc.isFiltered() ) {
            encoder.encodeStringsByRef(vc.getFilters());
        } else {
            encoder.encodeMissing(BCFType.INT32);
        }
    }

    private void buildInfo( VariantContext vc ) throws IOException {
        for ( Map.Entry<String, Object> infoFieldEntry : vc.getAttributes().entrySet() ) {
            encoder.encodeStringByRef(infoFieldEntry.getKey());
            // TODO -- use real type
            encoder.encodeString(infoFieldEntry.getValue().toString());
        }
    }

    private byte[] buildSamplesData(final VariantContext vc) throws IOException {
        // write size
        List<String> genotypeFields = StandardVCFWriter.calcVCFGenotypeKeys(vc);
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

        return encoder.getRecordBytes();
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

    private final BCFType getBCF2TypeFromHeader(final String field) {
        // TODO -- should take VC to determine best encoding
        final VCFCompoundHeaderLine metaData = VariantContext.getMetaDataForField(header, field);
        switch ( metaData.getType() ) {
            case Character: return BCFType.CHAR;
            //case Flag: return BCFType.FLAG;  // TODO -- HOW TO ENCODE FLAG?
            case String: return BCFType.CHAR;
            case Integer: return BCFType.INT32;
            case Float: return BCFType.FLOAT;
            default: throw new ReviewedStingException("Unexpected type for field" + field);
        }
    }

    private final void addGenotypeFilters(final VariantContext vc) throws IOException {
        encoder.startGenotypeField(VCFConstants.GENOTYPE_FILTER_KEY, 1, BCFType.CHAR);
        for ( final Genotype g : vc.getGenotypes() ) {
            if ( g.filtersWereApplied() && g.isFiltered() ) {
                encoder.encodeString(ParsingUtils.join(";", ParsingUtils.sortList(g.getFilters())));
            } else {
                encoder.encodeMissingValues(1, BCFType.CHAR); // todo fixme
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
        encoder.startGenotypeField(VCFConstants.GENOTYPE_KEY, requiredPloidy, BCFType.INT8);
        for ( final Genotype g : vc.getGenotypes() ) {
            if ( g.getPloidy() != requiredPloidy ) throw new ReviewedStingException("Cannot currently handle non-diploid calls!");
            final List<Integer> encoding = new ArrayList<Integer>(requiredPloidy);
            for ( final Allele a : g.getAlleles() ) {
                final int offset = a.isNoCall() ? 0 : (Byte.valueOf(alleleMap.get(a)) + 2);
                encoding.add(offset << 1 | (g.isPhased() ? 0x01 : 0x00));
            }
            encoder.encodeValues(encoding, BCFType.INT8);
        }
    }

    /**
     * Write the data in the encoder to the outputstream as a length encoded
     * block of data.  After this call the encoder stream will be ready to
     * start a new data block
     *
     * @throws IOException
     */
    private void writeBlock(final byte[] infoBlock, final byte[] genotypesBlock) throws IOException {
        BCFEncoder.encodePrimitive(infoBlock.length, BCFType.INT32, outputStream);
        BCFEncoder.encodePrimitive(genotypesBlock.length, BCFType.INT32, outputStream);
        outputStream.write(infoBlock);
        outputStream.write(genotypesBlock);
    }
}
