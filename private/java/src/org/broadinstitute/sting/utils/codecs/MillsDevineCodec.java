
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
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * a codec for the Mills/Devine indel tables from GR 2011 paper
 */
public class MillsDevineCodec implements FeatureCodec {

    private static final String DELETION_TYPE = "DEL";
    private static final String INSERTION_TYPE = "INS";

    // the minimum number of features in the CG file line
    private static final int minimumFeatureCount = 10;

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
        String[] array = line.split("\\t");

	//        System.out.println(array.length);
        // make sure the split was successful - that we got an appropriate number of fields
        //if ( array.length < minimumFeatureCount )
	//  return null;
	//System.out.println(array[2]);
        String INDEL_ID     = array[0];
        String CHR          = array[1];
        String START        = array[2];
        String STOP         = array[3];
        String IS_REVERSE   = array[4];
        String LEN          = array[5];
        String DOUBLE_HIT   = array[6];
        String DOUBLE_CENTER= array[7];
        String MEMBER       = array[8];
        String REF_TYPE     = array[9];
	       String CHIMP_TYPE   = array[10];
        String CELERA_TYPE  = array[11];
        String CHIMP_CHR    = array[12];
        String CHIMP_START  = array[13];
        String CHIMP_STOP   = array[14];
        String CELERA_CHR   = array[15];
        String CELERA_START = array[16];
        String CELERA_STOP  = array[17];
	       String NUM_TRACES   = array[18];

	           String GENE_NAME    = array[19];
        String GENE_TYPE    = array[20];
        String GENE_COLOR   = array[21];
        String IS_CODING    = array[22];
        String IS_DISEASE   = array[23];
        String IS_VALIDATED = array[24];
        String SEQ          = array[25];
        String CLASS        = array[26];
        String CLASS_TYPE   = array[27];
	        String HAS_TSD      = array[28];
        String TSD_LEN      = array[29];
        String TSD_LEFT_SEQ = array[30];
        String TSD_RIGHT_SEQ= array[31];
        String DBSNP_SUBMITTED     = array[32];
        String CONTIG       = array[33];
        String BIN          = array[34];


        long start = Long.valueOf(START)-1;

        //long end = Long.valueOf(array[3]);
        long end;

        Allele ref, alt = null;

        //System.out.println(line);
        if ( REF_TYPE.equals(INSERTION_TYPE) ) {
            ref = Allele.create(Allele.NULL_ALLELE_STRING, true);
            alt = Allele.create(SEQ.toUpperCase(), false);
            end = start;
        } else if ( REF_TYPE.equals(DELETION_TYPE) ) {
            ref = Allele.create(SEQ.toUpperCase(), true);
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

        attrs.put(VariantContext.ID_KEY, INDEL_ID);
        attrs.put("CLASS", CLASS);
        attrs.put("NUM_TRACES",NUM_TRACES);
        attrs.put("DOUBLE_HIT",DOUBLE_HIT);
        attrs.put("DOUBLE_CENTER",DOUBLE_CENTER);
        attrs.put("CLASS_TYPE",CLASS_TYPE);
        attrs.put("CHIMP_TYPE",CHIMP_TYPE);
        attrs.put("CELERA_TYPE",CELERA_TYPE);
        attrs.put("LEN",LEN);
        attrs.put("IS_CODING",IS_CODING);
        attrs.put("IS_DISEASE",IS_DISEASE);
        attrs.put("IS_VALIDATED",IS_VALIDATED);

        // create a new feature given the array
//        VariantContext vcCall = new VariantContext("UG_call", loc.getContig(), loc.getStart(), endLoc,
   //             myAlleles, genotypes, phredScaledConfidence/10.0, passesCallThreshold(phredScaledConfidence) ? null : filter, attributes, refContext.getBase());
        VariantContext vc =  new VariantContext("Mills", CHR, start, end, alleles,null, VariantContext.NO_NEG_LOG_10PERROR, null, attrs,"N".getBytes()[0]);
	    //System.out.println(vc.toString());
	/*        if(array[1].equals("3") ) {
	    System.out.format("%s %s %s\n",CHR,START,REF_TYPE);
            System.out.println(vc.toString());
	}
	*/  return vc;
    }

    public Class getFeatureType() {
        return VariantContext.class;
    }

    // There's no spec and no character to distinguish header lines...
    private final static int NUM_HEADER_LINES = 1;
    public Object readHeader(LineReader reader) {
        String headerLine = null;
        try {
            for (int i = 0; i < NUM_HEADER_LINES; i++)
                headerLine = reader.readLine();
        } catch (IOException e) {
            throw new IllegalArgumentException("Unable to read a line from the line reader");
        }
        return null;
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
