/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.utils.codecs;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A codec for the VAR file types produced by the Complete Genomics Institute
 *
 * <p>
 * The format specifies variation for a single sample per file.
 * The format includes several standard fields, but only the first eight appear to always be present.
 * Two lines per locus for variants, one for each haplotype.
 * </p>
 *
 * <p>
 * See also: @see <a href="http://vcftools.sourceforge.net/specs.html">VCF specification</a><br>
 * </p>
 *
 * <h2>File format example</h2>
 * <pre>
 *  1138    2       all     chr1    50593   50653   ref     =       =
 *  1139    2       1       chr1    50653   50654   snp     A       G       108             dbsnp.100:rs2691283
 *  1139    2       2       chr1    50653   50654   snp     A       G       108             dbsnp.100:rs2691283
 *  1140    2       all     chr1    50654   50691   ref     =       =
 *  1141    2       1       chr1    50691   50692   ref     C       C       103
 *  1141    2       2       chr1    50691   50692   snp     C       T       66              dbsnp.100:rs2691284
 *  1142    2       all     chr1    50692   51050   ref     =       =
 *  1143    2       all     chr1    51050   51085   no-call =       ?
 *  1144    2       all     chr1    51085   51152   ref     =       =
 *  1145    2       1       chr1    51152   51152   ins             G       47
 *  1145    2       2       chr1    51152   51152   ins             G       118
 *  1146    2       all     chr1    51152   51209   ref     =       =
 * </pre>
 *
 * @author Eric Banks
 * @since 2011
 */
public class CGVarCodec implements FeatureCodec {

    private static final String REF_TYPE = "ref";
    private static final String SNP_TYPE = "snp";
    private static final String DELETION_TYPE = "del";
    private static final String INSERTION_TYPE = "ins";
    private static final String SUBSTITUTION_TYPE = "sub";

    // the minimum number of features in the CG file line
    private static final int minimumFeatureCount = 8;

    /**
     * decode the location only
     * @param line the input line to decode
     * @return a HapMapFeature
     */
    public Feature decodeLoc(String line) {
        return decode(line);
    }

    /**
     * decode the CG record
     * @param line the input line to decode
     * @return a VariantContext
     */
    public Feature decode(String line) {
        String[] array = line.split("\\s+");

        // make sure the split was successful - that we got an appropriate number of fields
        if ( array.length < minimumFeatureCount )
            return null;

        String type = array[6];

        long start = Long.valueOf(array[4]);
        long end;
        Allele ref, alt = null;

        //System.out.println(line);

        if ( type.equals(SNP_TYPE) ) {
            ref = Allele.create(array[7], true);
            alt = Allele.create(array[8], false);
            end = start;
        } else if ( type.equals(INSERTION_TYPE) ) {
            ref = Allele.create(Allele.NULL_ALLELE_STRING, true);
            alt = Allele.create(array[7], false);
            end = start;
        } else if ( type.equals(DELETION_TYPE) ) {
            ref = Allele.create(array[7], true);
            alt = Allele.create(Allele.NULL_ALLELE_STRING, false);
            end = start + ref.length();
        //} else if ( type.equals(REF_TYPE) ) {
        //    ref = Allele.create("N", true); // ref bases aren't accurate
        //    start++;
        //    end = start;
        //} else if ( type.equals(SUBSTITUTION_TYPE) ) {
        //    ref = Allele.create(array[7], true);
        //    alt = Allele.create(array[8], false);
        //    end = start + Math.max(ref.length(), alt.length());
        } else {
            return null; // we don't handle other types
        }

        HashSet<Allele> alleles = new HashSet<Allele>();
        alleles.add(ref);
        if ( alt != null )
            alleles.add(alt);

        HashMap<String, Object> attrs = new HashMap<String, Object>();
        String id = array[array.length - 1];
        if ( id.indexOf("dbsnp") != -1 ) {
            attrs.put(VariantContext.ID_KEY, parseID(id));
        }

        // create a new feature given the array
        return new VariantContext("CGI", array[3], start, end, alleles, VariantContext.NO_NEG_LOG_10PERROR, null, attrs);
    }

    public Class<VariantContext> getFeatureType() {
        return VariantContext.class;
    }

    // There's no spec and no character to distinguish header lines...
    private final static int NUM_HEADER_LINES = 12;
    public Object readHeader(LineReader reader) {
        return null;

        //String headerLine = null;
        //try {
        //    for (int i = 0; i < NUM_HEADER_LINES; i++)
        //        headerLine = reader.readLine();
        //} catch (IOException e) {
        //    throw new IllegalArgumentException("Unable to read a line from the line reader");
        //}
        //return headerLine;
    }

    private static final Pattern DBSNP_PATTERN = Pattern.compile("^dbsnp\\.\\d+:(.*)");
    private String parseID(String raw) {
        StringBuilder sb = null;

        String[] ids = raw.split(";");
        for ( String id : ids ) {
            Matcher matcher = DBSNP_PATTERN.matcher(id);
            if ( matcher.matches() ) {
                String rsID = matcher.group(1);
                if ( sb == null ) {
                    sb = new StringBuilder(rsID);
                } else {
                    sb.append(";");
                    sb.append(rsID);
                }
            }
        }

        return sb == null ? null : sb.toString();
    }
}
