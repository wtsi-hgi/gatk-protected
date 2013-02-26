/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.utils.codecs;

import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * TODO GUILLERMO DEL ANGEL
 *
 * <p>
 * a codec for the Mills/Devine indel tables from GR 2011 paper
 * </p>
 *
 * <p>
 * See also: @see <a href="http://vcftools.sourceforge.net/specs.html">VCF specification</a><br>
 * </p>

 * </p>
 *
 * <h2>File format example</h2>
 * <pre>
 *     line 1
 *     line 2
 *     line 3
 * </pre>
 *
 * @author Guillermo del Angel
 * @since 2010
 */
public class MillsDevineCodec extends AsciiFeatureCodec<VariantContext> {

    private static final String DELETION_TYPE = "DEL";
    private static final String INSERTION_TYPE = "INS";

    // the minimum number of features in the CG file line
    private static final int minimumFeatureCount = 10;


    public MillsDevineCodec() {
        super(VariantContext.class);
    }

    /**
     * decode the CG record
     * @param line the input line to decode
     * @return a VariantContext
     */
    public VariantContext decode(String line) {
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
        // We need to give it an arbitrary reference padding base
        if ( REF_TYPE.equals(INSERTION_TYPE) ) {
            ref = Allele.create("N", true);
            alt = Allele.create("N" + SEQ.toUpperCase(), false);
            end = start;
        } else if ( REF_TYPE.equals(DELETION_TYPE) ) {
            ref = Allele.create("N" + SEQ.toUpperCase(), true);
            alt = Allele.create("N", false);
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
        VariantContext vc =  new VariantContextBuilder("Mills", CHR, start, end, alleles).id(INDEL_ID).attributes(attrs).make();
	    //System.out.println(vc.toString());
	/*        if(array[1].equals("3") ) {
	    System.out.format("%s %s %s\n",CHR,START,REF_TYPE);
            System.out.println(vc.toString());
	}
	*/  return vc;
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
