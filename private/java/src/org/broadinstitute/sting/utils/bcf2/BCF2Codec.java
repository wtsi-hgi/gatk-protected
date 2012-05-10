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

import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

public class BCF2Codec implements FeatureCodec<VariantContext> {
    final protected static Logger logger = Logger.getLogger(FeatureCodec.class);
    private VCFHeader header = null;
    private final ArrayList<String> contigNames = new ArrayList<String>();
    private final ArrayList<String> dictionary = new ArrayList<String>();
    BCF2Decoder decoder;

    // ----------------------------------------------------------------------
    //
    // Feature codec interface functions
    //
    // ----------------------------------------------------------------------

    public Feature decodeLoc( final PositionalBufferedStream inputStream ) {
        return decode(inputStream);
        // TODO: a less expensive version of decodeLoc() that doesn't use VariantContext
        // TODO: very easy -- just decodeImplicitBlock, and then skip to end of end of sites block
        // TODO: and then skip genotypes block
    }

    public VariantContext decode( final PositionalBufferedStream inputStream ) {
        decoder.readNextBlock(inputStream);
        final VariantContextBuilder builder = new VariantContextBuilder();

        decodeImplicitBlock(builder);
        decodeID(builder);
        final ArrayList<Allele> alleles = decodeAlleles(builder);
        decodeFilter(builder);
        decodeInfo(builder);

        // TODO: call readNextBlock here to get genotypes data
        decodeGenotypes(alleles, builder);

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

        // read the header
        this.header = AbstractVCFCodec.parseHeader(headerLines, VCFHeaderVersion.VCF4_1);

        // create the config offsets
        for ( final VCFContigHeaderLine contig : header.getContigLines())
            contigNames.add(contig.getID());

        // create the string dictionary
        parseDictionary(header);

        // prepare the decoder
        decoder = new BCF2Decoder(dictionary);

        // position right before next line (would be right before first real record byte at end of header)
        return new FeatureCodecHeader(header, inputStream.getPosition());
    }

    @Override
    public boolean canDecode( final String path ) {
        try {
            FileInputStream fis = new FileInputStream(path);
            AsciiLineReader reader = new AsciiLineReader(new PositionalBufferedStream(fis));
            String firstLine = reader.readLine();
            if ( firstLine != null && firstLine.equals(BCF2Constants.VERSION_LINE) ) {
                return true;
            }
        } catch ( FileNotFoundException e ) {
            return false;
        } catch ( IOException e ) {
            return false;
        }

        return false;
    }

    private final void parseDictionary(final VCFHeader header) {
        for ( final VCFHeaderLine line : header.getMetaData() ) {
            if ( line.getKey().equals(BCF2Constants.DICTIONARY_LINE_TAG) ) {
                for ( final String string : line.getValue().split(BCF2Constants.DICTIONARY_LINE_ENTRY_SEPARATOR) )
                    dictionary.add(string);
                break;
            }
        }

        // if we got here we never found a dictionary, or there are no elements in the dictionary
        if ( dictionary.size() == 0 )
            throw new UserException.MalformedBCF2("Dictionary header element was absent or empty");
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

    private final void decodeImplicitBlock(final VariantContextBuilder builder) {
        final int contigOffset = decoder.decodeInt(BCFType.INT32.getSizeInBytes());
        final String contig = lookupContigName(contigOffset);
        builder.chr(contig);

        final int pos = decoder.decodeInt(BCFType.INT32.getSizeInBytes());
        final int refLength = decoder.decodeInt(BCFType.INT32.getSizeInBytes());
        builder.start((long)pos);
        builder.stop((long)(pos + refLength - 1)); // minus one because of our open intervals

        final int qualAsInt = decoder.decodeRawFloat();
        if ( ! decoder.rawFloatIsMissing(qualAsInt) ) {
            builder.log10PError(decoder.rawFloatToFloat(qualAsInt) / -10.0);
        }
    }

    private void decodeID( final VariantContextBuilder builder ) {
        final String id = (String)decoder.decodeTypedValue();

        if ( id == null ) {
            builder.noID();
        }
        else {
            builder.id(id);
        }
    }

    public static ArrayList<Allele> clipAllelesIfNecessary(int position, String ref, ArrayList<Allele> unclippedAlleles) {
        if ( AbstractVCFCodec.isSingleNucleotideEvent(unclippedAlleles) ) {
            ArrayList<Allele> clippedAlleles = new ArrayList<Allele>(unclippedAlleles.size());
            AbstractVCFCodec.clipAlleles(position, ref, unclippedAlleles, clippedAlleles, -1);
            return clippedAlleles;
        } else
            return unclippedAlleles;
    }

    private ArrayList<Allele> decodeAlleles( final VariantContextBuilder builder ) {
        final String ref = (String)decoder.decodeTypedValue();
        builder.referenceBaseForIndel(ref.getBytes()[0]);
        final Object altObject = decoder.decodeTypedValue(); // alt can be either atomic or vectors

        // TODO -- probably need inline decoder for efficiency here (no sense in going bytes -> string -> vector -> bytes
        ArrayList<Allele> alleles = new ArrayList<Allele>(2);
        alleles.add(Allele.create(ref, true));

        for ( final String alt : asStrings(altObject)) {
            alleles.add(Allele.create(alt));
        }

        alleles = clipAllelesIfNecessary((int)builder.getStart(), ref, alleles);
        builder.alleles(alleles);
        return alleles;
    }

    private void decodeFilter( final VariantContextBuilder builder ) {
        final Object filters = decoder.decodeTypedValue();

        if ( filters == null ) {
            builder.unfiltered();
        }
        else {
            builder.filters(new LinkedHashSet<String>(asStrings(filters)));
        }
    }

    private void decodeInfo( final VariantContextBuilder builder ) {
        final int numInfoFields = (Integer)decoder.decodeRequiredTypedValue("Number of info fields");
        final Map<String, Object> infoFieldEntries = new HashMap<String, Object>(numInfoFields);

        for ( int i = 0; i < numInfoFields; i++ ) {
            final String key = (String)decoder.decodeTypedValue();
            final Object value = decoder.decodeTypedValue();
            infoFieldEntries.put(key, value);
        }

        builder.attributes(infoFieldEntries);
    }

    private void decodeGenotypes( ArrayList<Allele> siteAlleles, final VariantContextBuilder builder ) {
        final List<String> samples = new ArrayList<String>(header.getGenotypeSamples());
        final int nSamples = samples.size();
        final int nFields = (Integer)decoder.decodeTypedValue();
        final Map<String, List<Object>> fieldValues = decodeGenotypeFieldValues(nFields, nSamples);

        final List<Genotype> genotypes = new ArrayList<Genotype>(nSamples);
        for ( int i = 0; i < nSamples; i++ ) {
            final String sampleName = samples.get(i);
            List<Allele> alleles = null;
            boolean isPhased = false;
            double log10PError = VariantContext.NO_LOG10_PERROR;
            Set<String> filters = null;
            Map<String, Object> attributes = null;
            double[] log10Likelihoods = null;

            for ( final Map.Entry<String, List<Object>> entry : fieldValues.entrySet() ) {
                final String field = entry.getKey();
                final List<Object> values = entry.getValue();
                if ( field.equals(VCFConstants.GENOTYPE_KEY) ) {
                    alleles = decodeGenotypeAlleles(siteAlleles, (List<Integer>)values.get(i));
                } else if ( field.equals(VCFConstants.GENOTYPE_QUALITY_KEY) ) {
                    final Integer value = (Integer)values.get(i);
                    if ( value != BCFType.INT8.getMissingJavaValue() )
                        log10PError = value / -10.0;
                } else if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY) ) {
                    filters = new HashSet<String>((List<String>)values.get(i));
                } else { // add to attributes
                    if ( attributes == null ) attributes = new HashMap<String, Object>(nFields);
                    attributes.put(field, values.get(i));
                }
            }

            if ( alleles == null ) throw new ReviewedStingException("BUG: no alleles found");

            final Genotype g = new Genotype(sampleName, alleles, log10PError, filters, attributes, isPhased, log10Likelihoods);
            genotypes.add(g);
        }

        builder.genotypes(genotypes);
    }

    private final List<Allele> decodeGenotypeAlleles(final ArrayList<Allele> siteAlleles, final List<Integer> encoded) {
        final List<Allele> gt = new ArrayList<Allele>(encoded.size());
        for ( final int encode : encoded ) {
            final int offset = encode >> 1;
            if ( offset == 0 )
                gt.add(Allele.NO_CALL);
                else if ( offset > 1 ) // skip absent
                gt.add( siteAlleles.get(offset - 2));
        }
        return gt;
    }

    private final Map<String, List<Object>> decodeGenotypeFieldValues(final int nFields, final int nSamples) {
        Map<String, List<Object>> map = new LinkedHashMap<String, List<Object>>(nFields);

        for ( int i = 0; i < nFields; i++ ) {
            final String field = (String)decoder.decodeTypedValue();
            final byte typeDescriptor = decoder.readTypeDescriptor();
            final List<Object> values = new ArrayList<Object>(nSamples);
            for ( int j = 0; j < nSamples; j++ )
                values.add(decoder.decodeTypedValue(typeDescriptor));
            map.put(field, values);
        }

        return map;
    }

    private String lookupContigName( final int contigOffset ) {
        if ( contigOffset < contigNames.size() ) {
            return contigNames.get(contigOffset);
        }
        else {
            throw new UserException.MalformedBCF2(String.format("No contig at index %d present in the sequence dictionary from the BCF2 header (%s)", contigOffset, contigNames));
        }
    }


    // ----------------------------------------------------------------------
    //
    // Utility functions
    //
    // ----------------------------------------------------------------------

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
}
