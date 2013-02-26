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

package org.broadinstitute.sting.gatk.features.maf;

import org.apache.log4j.Logger;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.lang.reflect.Field;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Jan 24, 2011
 */

public class MafCodec extends AsciiFeatureCodec<MafFeature> {
     private final static Logger log = Logger.getLogger(MafCodec.class);

     private int expectedTokenCount = -1;


     private Column BUILD_COL = new Column(new String[]{"NCBI_Build","build"},true);
     private Column CHR_COL = new Column(new String[] {"Chromosome","chr"},true);
     private Column START_COL = new Column(new String[] {"Start_position","start"},true);
     private Column END_COL = new Column(new String[]{"End_position","end"},true);
     private Column REF_ALLELE_COL = new Column(new String[] {"Reference_Allele","ref_allele"},true);
     private Column TUMOR_ALLELE1_COL = new Column(new String[] {"Tumor_Seq_Allele1","tum_allele1"},true);
     private Column TUMOR_ALLELE2_COL = new Column(new String[] {"Tumor_Seq_Allele2","tum_allele2"},true);
     private Column TUMOR_SAMPLE_COL = new Column(new String[] {"Tumor_Sample_Barcode","tumor_barcode"},true);
     private Column NORMAL_SAMPLE_COL = new Column(new String[]{"Matched_Norm_Sample_Barcode","normal_barcode"},true);
    // optional fields (absent from maf lite):
     private Column VARTYPE_COL = new Column(new String[]{"Variant_Type","classification"},false);
     private Column STRAND_COL = new Column(new String[]{"Strand","strand"},false);
     private Column HUGO_GENE_COL = new Column(new String[]{"Hugo_Symbol","gene"},false);
     private Column VARCLASS_COL = new Column(new String[]{"Variant_Classification","type"},false);


     public enum MAF_TYPE {
        UNKNOWN,LITE, ANNOTATED
     }

     private static String INS ="INS";
     private static String DEL ="DEL";
     private static String SNP ="SNP";
     private static String MNP ="MNP";

     private MAF_TYPE mafType=MAF_TYPE.UNKNOWN;

     private List<Column> allColumns = null; /// filled dynamically by constructor through introspection. Slow but less typing.

    private boolean tooManyColsWarned = false;
    private boolean tooFewColsWarned = false;

     public MafCodec() {
         super(MafFeature.class);
         allColumns = new ArrayList<Column>(30);
         Field[] fields = this.getClass().getDeclaredFields();
         try {
            for ( Field f : fields ) {
                 if ( f.get(this) instanceof Column ) {
                     allColumns.add((Column)f.get(this));
                 }
            }
         } catch (IllegalAccessException e) {
             throw new StingException("Error in MAFCodec when trying to introspect itself, this is probably a BUG",e);
         }
     }


    /**
     * Decode a line to obtain just its FeatureLoc for indexing -- contig, start, and stop.
     * This method will NOT fill in the additional information available in the maf file
     * @param line the input line to decode
     * @return  Return the FeatureLoc encoded by the line, or null if the line does not represent a feature (e.g. is
     * a comment)
     */
    public Feature decodeLoc(String line) {
           return reallyDecode(line,false);
    }


    /**
     * Fully decode a line, will try extracting as much additional/annotation information from the maf file as it can.
     * @param line the input line to decode
     * @return  Return the FeatureLoc encoded by the line, or null if the line does not represent a feature (e.g. is
     * a comment)
     */
    public MafFeature decode(String line) {
        return reallyDecode(line,true);
    }

    /** Decodes a maf line. If <code>extra</code> is false, will decode only location and return;
     * if <code>extra</code> is true, then extracts everything it can (samples, annotations, etc)
     * @param line
     * @param extra
     * @return
     */
    public MafFeature reallyDecode(String line, boolean extra) {

        // ignore commented-out lines
        if (line.startsWith("#")) return null;

        // split the line
        String[] tokens = line.split("\\t",-1);

        if ( expectedTokenCount == -1 ) { // do this only when we receive the first line and do not know the number of columns yet
             // we have not seen a single line yet, let's initialize the number of fields from the first line:
             expectedTokenCount = tokens.length;
             log.info("MAF: line has "+expectedTokenCount+" fields (columns)");
             if ( expectedTokenCount == 9 ) {
                mafType = MAF_TYPE.LITE;
                log.info("MAF file appears to be MAF Lite");
             } else {
                 if ( expectedTokenCount >= 63 ) {
                     mafType = MAF_TYPE.ANNOTATED;
                     log.info("MAF file appears to be MAF-Annotated");
                 } else {
                     log.info("MAF file has "+expectedTokenCount +" columns in first line, unknown file type");
                 }
             }
             if ( line.contains("Chromosome") && line.contains("Start") && line.contains("Build") ||
                   line.contains("build") && line.contains("start") && line.contains("ref_allele")  ) {
                // a naive way to detect the line with column names

                 setColumnsFromHeader(tokens);
                 log.info("MAF file contains header, all required columns found!");
                 return null;
            } else {
                 switch( mafType ) {
                     case UNKNOWN: throw new UserException.MalformedFile("Can not guess type of the MAF file from number of columns and there is no header");
                     case LITE: setMafLiteCols(); break;
                     case ANNOTATED: setMafAnnotatedCols(); break;
                 }
                 log.info("MAF file has no header; assuming standard column order for the MAF type "+mafType);
            }
        }


        if (tokens.length < expectedTokenCount) {
            if ( ! tooFewColsWarned ) {
                log.error("MAF line contains too few columns ("+tokens.length+"); this error is reported only once.");
                tooFewColsWarned = true;
            }
        }
        if (tokens.length > expectedTokenCount) {
            if ( ! tooManyColsWarned ) {
                log.warn("MAF line contains more columns than expected ("+tokens.length+"); extra columns discarded. This error is shown only once.");
                tooManyColsWarned = true;
            }
        }

        if ( tokens[CHR_COL.getIndex()].equals("Chromosome") || tokens[CHR_COL.getIndex()].equals("chr")) return null; // if someone uses this codec manually and feeds it the header line multiple times...
        // create a new feature from the line:

        int start = 0;
        try {
            start = Integer.parseInt(START_COL.getValue(tokens));
        } catch (NumberFormatException e) {
            throw new UserException.MalformedFile("Missing or non-numeric start position in line:\n"+line,e);
        }
        int stop = 0 ;
        try {
            stop = Integer.parseInt(END_COL.getValue(tokens));
        } catch (NumberFormatException e) {
            throw new UserException.MalformedFile("Missing or non-numeric stop position in line:\n"+line,e);
        }

        String eventType="UNKNOWN";

        String ref = REF_ALLELE_COL.getValue(tokens);
        String alt1 = TUMOR_ALLELE1_COL.getValue(tokens);
        String alt2 = TUMOR_ALLELE2_COL.getValue(tokens);

        if ( ref.equals("-") ) {
            // insertion
            eventType = INS;
            stop-- ; // maf lists stop as first base after insertion, convert internally to vcf style
            // perform some format validation:

            if ( alt1.equals("-") && alt2.equals("-") )
                throw new UserException.MalformedFile("Inconsistency in MAF: both alt alleles reported as ref ('-') for an insertion");

            if ( ! alt1.equals("-") && ! alt2.equals("-") && ! alt1.equals(alt2) )
                throw new UserException.MalformedFile("Inconsistency in MAF: two different (non-ref) alt alleles reported for an insertion");

            if ( stop != start )
                throw new UserException.MalformedFile("Inconsistency in MAF: end position for an insertion is not start+1");

        } else {
            if ( alt1.equals("-") || alt2.equals("-") ) {
                // deletion
                eventType = DEL;
                start--; // maf lists start as the first deleted base; convert internally to vcf style
                // perform some format validation:

                if ( ! alt1.equals("-") && ! alt1.equals(ref) )
                    throw new UserException.MalformedFile("Inconsistency in MAF: non-deleted alt allele is not ref for a deletion");

                if ( ! alt2.equals("-") && ! alt2.equals(ref) )
                    throw new UserException.MalformedFile("Inconsistency in MAF: non-deleted alt allele is not ref for a deletion");

                if ( (stop - start) != ref.length() )
                    throw new UserException.MalformedFile("Inconsistency in MAF: deletion length is not end-start+1");

            } else {
                // no '-' alleles --> it's a snp/mnp
                if ( ref.length() == 1 ) {
                    // it's a snp
                    eventType = SNP;
                    if ( stop != start )
                        throw new UserException.MalformedFile("Inconsistency in MAF: start/end positions not equal for a SNP");
                } else {
                    // it's an mnp
                    eventType = MNP;
                    if ( (stop - start + 1) != ref.length() )
                        throw new UserException.MalformedFile("Inconsistency in MAF: MNP length is not end-start+1");
                }

                if ( alt1.length() != ref.length() || alt2.length() != ref.length() )
                    throw new UserException.MalformedFile("Inconsistency in MAF: lengths of ref and alt alleles for a SNP/MNP differ");
                if ( ! alt1.equals(ref) && ! alt2.equals(ref) && ! alt1.equals(alt2) )
                    throw new UserException.MalformedFile("Inconsistency in MAF: two different non-ref alt alleles reported for a SNP/MNP");
            }
        }
        // if we got vartype column, make sure it makes sense:
        if ( VARTYPE_COL.isSet(tokens) && ! tokens[VARTYPE_COL.getIndex()].equals(eventType) )  {
            // special case: we annotate everything as MNP while MAF can have DNP/TNP, these are fine:
            if ( eventType == MNP && (
                    tokens[VARTYPE_COL.getIndex()].equals("DNP") && ref.length() == 2 ||
                    tokens[VARTYPE_COL.getIndex()].equals("TNP") && ref.length() == 3)
                    ) {}                                                              // these are fine
            else {
                throw new UserException.MalformedFile("Inconsistency in MAF: variant looks like a "+eventType +" but annotated as "+
                    tokens[VARTYPE_COL.getIndex()]);
            }
        }
        MafFeature feature = new MafFeature(CHR_COL.getValue(tokens),start,stop);

        if ( ! extra ) return feature; // ignore additional fields unless we were explicitly asked to read those!
        
        feature.setVariantType(eventType);
        feature.setRefAllele(ref);
        feature.setObservedTumor(alt1,alt2);
        feature.setTumorSample(TUMOR_SAMPLE_COL.getValue(tokens));
        feature.setNormalSample(NORMAL_SAMPLE_COL.getValue(tokens));

        if ( HUGO_GENE_COL.isSet(tokens) ) feature.setHugoGeneSymbol(tokens[HUGO_GENE_COL.getIndex()]);
        if ( VARCLASS_COL.isSet(tokens) ) feature.setVariantClassification(tokens[VARCLASS_COL.getIndex()]);

        return feature;
    }

    /** Set expected column indices for MafLite
     *
     */
    private void setMafLiteCols() {
        BUILD_COL.setIndex(0);
        CHR_COL.setIndex(1);
        START_COL.setIndex(2);
        END_COL.setIndex(3);
        REF_ALLELE_COL.setIndex(4);
        TUMOR_ALLELE1_COL.setIndex(5);
        TUMOR_ALLELE2_COL.setIndex(6);
        TUMOR_SAMPLE_COL.setIndex(7);
        NORMAL_SAMPLE_COL.setIndex(8);
    }

    private void setMafAnnotatedCols() {
        BUILD_COL.setIndex(3);
        CHR_COL.setIndex(4);
        START_COL.setIndex(5);
        END_COL.setIndex(6);
        REF_ALLELE_COL.setIndex(10);
        TUMOR_ALLELE1_COL.setIndex(11);
        TUMOR_ALLELE2_COL.setIndex(12);
        TUMOR_SAMPLE_COL.setIndex(15);
        NORMAL_SAMPLE_COL.setIndex(16);
        VARTYPE_COL.setIndex(9);
        STRAND_COL.setIndex(7);
        VARCLASS_COL.setIndex(8);
        HUGO_GENE_COL.setIndex(0);
    }

    private void setColumnsFromHeader(String[] tokens) {
        Map<String,Integer> colNames = new HashMap<String,Integer>();
        for ( int i = 0 ; i < tokens.length ; i++ ) colNames.put(tokens[i],i);

        for ( Column c : allColumns ) c.setFromMap(colNames);
    }


}


class Column {
    int index ;
    List<String> names;
    boolean required;

    Column(String name, boolean required) {
        this.names = new ArrayList<String>();
        this.names.add(name);
        this.required = required;
        this.index = -1;
    }

    Column(String [] names, boolean required) {
        this.names = new ArrayList<String>();
        for ( int i = 0 ; i < names.length ; i++ ) this.names.add(names[i]);
        this.required = required;
        this.index = -1;
    }

    public String getName() { return names.get(0); }
    public Collection<String> getNames() { return names; }
    public void setName(String name) {
        for ( String n : names ) {
            if ( n.equals(name) ) return;
        }
        this.names.add( name );
    }

    public int getIndex() { return index; }
    public void setIndex(int index) { this.index = index; }
    public String getValue(String[] fields) {
        if ( index < fields.length ) return fields[index];

        if ( required ) throw new UserException.MalformedFile("In MAF file: required column "+getName()+" has index "+index+
                    ", but only "+fields.length+ " fields are present in maf line");
        return null;
    }

    /** Sets this column's index from the provided name->index map (i.e. searches for itself in the map).
     * If column not found, <code>throw_exception</code> is true <i>AND</i> this column is required, then an exception will
     * be thrown right away; otherwise returns quietely even if map does not contain this column.
     * @param m
     * @param throw_exception
     */
    public void setFromMap(Map<String,Integer> m, boolean throw_exception) {
//        Integer i = null;
        for ( String n : names ) {
            if ( m.containsKey(n) ) {
//                if ( i != null )
//                    throw new UserException.MalformedFile("MAF file contains multiple columns with name or alternative names registered for single data field "+getID());
                // go with the first column name found; we assume here that column names have priorities:
                // for instance, if the file has both 'Chromosome' and 'chr' columns, we will just take
                // Chromosome and run with that
                this.index = m.get(n);
                return;
            }
        }
        if ( this.required && throw_exception ) throw new UserException.MalformedFile("Required column "+getName()+" is missing from the maf file");
        this.index = -1;
    }

/**  Sets this column's index from the provided name->index map (i.e. searches for itself in the map).
 * If this column is required but not found in the map, then an exception will
 * be thrown.
 * @param m
 */
    public void setFromMap(Map<String,Integer> m) {
        setFromMap(m,true);
    }

    public boolean isSet() { return index > -1; }

    public boolean isSet(String[] fields) { return index > -1 && index < fields.length; }

}

