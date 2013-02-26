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

package org.broadinstitute.sting.gatk.walkers.annotator.genomicannotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableFeature;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.*;
import java.util.Map.Entry;

/**
 * This plugin for {@link VariantAnnotatorEngine} serves as the core
 * of the {@link GenomicAnnotator}. It finds all records in the -B input files
 * that match the given variant's position and, optionally, the variant's reference and alternate alleles.
 *
 * For details, see:  http://www.broadinstitute.org/gsa/wiki/index.php/GenomicAnnotator
 */
public class GenomicAnnotation extends InfoFieldAnnotation {

    public static final String CHR_COLUMN = "chr";
    public static final String START_COLUMN = "start";
    public static final String END_COLUMN = "end";
    public static final String HAPLOTYPE_REFERENCE_COLUMN = "haplotypeReference";
    public static final String HAPLOTYPE_ALTERNATE_COLUMN = "haplotypeAlternate";

    public static final String NUM_MATCHES_SPECIAL_INFO_FIELD = "numMatchingRecords";

    /** Characters that aren't allowed within VCF info field key-value pairs */
    public static final char[] ILLEGAL_INFO_FIELD_VALUES            = {  ' ', '=', ';' };
    /** Replacement for each character in ILLEGAL_INFO_FIELD_VALUES */
    public static final char[] ILLEGAL_INFO_FIELD_VALUE_SUBSTITUTES = {  '_', '-',  '!'  };


    private void modifyAnnotationsForIndels(VariantContext vc, String featureName, Map<String, String> annotationsForRecord) {
        String inCodingRegionKey = featureName + ".inCodingRegion";
        String referenceCodonKey = featureName + ".referenceCodon";
        String variantCodonKey = featureName + ".variantCodon";
        String codingCoordStrKey = featureName + ".codingCoordStr";
        String proteinCoordStrKey = featureName + ".proteinCoordStr";
        String haplotypeReferenceKey = featureName + "." + HAPLOTYPE_REFERENCE_COLUMN;
        String haplotypeAlternateKey = featureName + "." + HAPLOTYPE_ALTERNATE_COLUMN;
        String functionalClassKey = featureName + ".functionalClass";
        String startKey = featureName + "." + START_COLUMN;
        String endKey = featureName + "." + END_COLUMN;
        String referenceAAKey = featureName + ".referenceAA";
        String variantAAKey = featureName + ".variantAA";
        String changesAAKey = featureName + ".changesAA";

        annotationsForRecord.put(variantCodonKey, "unknown");
        annotationsForRecord.put(codingCoordStrKey, "unknown");
        annotationsForRecord.put(proteinCoordStrKey, "unknown");
        annotationsForRecord.put(referenceAAKey, "unknown");
        annotationsForRecord.put(variantAAKey, "unknown");

        String refAllele = vc.getReference().getDisplayString();
        if (refAllele.length() == 0) { refAllele = "-"; }

        String altAllele = vc.getAlternateAllele(0).toString();
        if (altAllele.length() == 0) { altAllele = "-"; }

        annotationsForRecord.put(haplotypeReferenceKey, refAllele);
        annotationsForRecord.put(haplotypeAlternateKey, altAllele);
        annotationsForRecord.put(startKey, String.format("%d", vc.getStart()));
        annotationsForRecord.put(endKey, String.format("%d", vc.getEnd()));

        boolean isCodingRegion = annotationsForRecord.containsKey(inCodingRegionKey) && annotationsForRecord.get(inCodingRegionKey).equalsIgnoreCase("true") ? true : false;
        boolean isFrameshift = (vc.getIndelLengths().get(0) % 3 == 0) ? false : true;

        String functionalClass;
        if (isCodingRegion) {
            functionalClass = isFrameshift ? "frameshift" : "inframe";
            annotationsForRecord.put(changesAAKey, "true");
        } else {
            functionalClass = "noncoding";
        }

        annotationsForRecord.put(functionalClassKey, functionalClass);
    }

    /**
     * For each -B input file, for each record which overlaps the current locus, generates a
     * set of annotations of the form:
     *
     * bindingName.columnName1=columnValue, bindingName.columnName2=columnValue2, etc.
     *
     * For example: dbSNP.avHet=0.7, dbSNP.ref_allele=A, etc.
     *
     * @return The following is an explanation of this method's return value:
     *
     * The annotations from a matching in a particular file are stored in a Map<String, String>
     * where the key is bindingName.columnName and the value is the columnValue.
     * Since a single input file can have multiple records that overlap the current
     * locus (eg. dbSNP can have multiple entries for the same genomic position), a different
     * Map<String, String> is created for each matching record in a particular file.
     * The set of matching records for each file is then represented as a List<Map<String, String>>
     *
     * The return value of this method is a Map<String, Object> of the form:
     *     rodName1 -> List<Map<String, String>>
     *     rodName2 -> List<Map<String, String>>
     *     rodName3 -> List<Map<String, String>>
     *     ...
     * Where the rodNames are the -B binding names for each file that were specified on the command line (eg. -B bindingName,AnnotatorInputTable,/path/to/file).
     *
     * NOTE: The lists (List<Map<String, String>>) are guaranteed to have size > 0
     * because a  rodName -> List<Map<String, String>>  entry will only
     * be created in Map<String, Object> if the List has at least one element.
     */
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
            final ReferenceContext ref,
            final Map<String, AlignmentContext> stratifiedContexts,
            final VariantContext vc) {

        //iterate over each record that overlaps the current locus, and, if it passes certain filters,
        //add its values to the list of annotations for this locus.
        final Map<String, Object> annotations = new HashMap<String, Object>();
        for(final GATKFeature gatkFeature : tracker.getAllValuesAsGATKFeatures())
        {
            final String name = gatkFeature.getName();
            if( name.equals("variant") || name.equals("interval") ) {
                continue;
            }

            if( ! (gatkFeature.getUnderlyingObject() instanceof AnnotatorInputTableFeature) ) {
                continue; //GenericAnnotation only works with TabularRODs because it needs to be able to select individual columns.
            }

            final Map<String, String> annotationsForRecord = convertRecordToAnnotations( gatkFeature.getName(), ((AnnotatorInputTableFeature) gatkFeature.getUnderlyingObject()).getColumnValues());

            //If this record contains the HAPLOTYPE_REFERENCE_COLUMN and/or HAPLOTYPE_ALTERNATE_COLUMN, check whether the
            //alleles specified match the the variant's reference allele and alternate allele.
            //If they don't match, this record will be skipped, and its values will not be used for annotations.
            //
            //If one of these columns doesn't exist in the current rod, or if its value is * (star), then this is treated as an automatic match.
            //Otherwise, the HAPLOTYPE_REFERENCE_COLUMN is only considered to be matching the variant's reference if the string values of the two
            //are exactly equal (case-insensitive).

            //The HAPLOTYPE_REFERENCE_COLUMN matches the variant's reference allele based on a case-insensitive string comparison.
            //The HAPLOTYPE_ALTERNATE_COLUMN can optionally list more than allele separated by one of these chars: ,\/:|
            // only check this value for SNPs
            String hapAltValue = vc.isSNP() ? annotationsForRecord.get( generateInfoFieldKey(name, HAPLOTYPE_ALTERNATE_COLUMN) ) : null;
            if ( hapAltValue != null && !hapAltValue.equals("*") ) {
                Set<Allele> alternateAlleles = vc.getAlternateAlleles();
                //if(alternateAlleles.isEmpty()) {
                    //handle a site that has been called monomorphic reference
                    //alternateAlleles.add(vc.getReference());
                    //continue;            //TODO If this site is monomorphic in the VC, and the current record specifies a particular alternate allele, skip this record. Right?
                //} else
                if(alternateAlleles.size() > 1) {
                    throw new UserException.MalformedFile("File associated with " + vc.getSource() + " contains record [" + vc + "] contains " + alternateAlleles.size() + " alternate alleles. GenomicAnnotion currently only supports annotating 1 alternate allele.");
                }

                Allele vcAlt;
                if(alternateAlleles.isEmpty()) {
                    vcAlt = vc.getReference();
                } else {
                    vcAlt = alternateAlleles.iterator().next();
                }

                boolean matchFound = false;
                for(String hapAlt : hapAltValue.split("[,\\\\/:|]")) {

                    if(!hapAlt.isEmpty() && vcAlt.basesMatch(hapAlt)) {
                        matchFound = true;
                        break;
                    }
                }
                if(!matchFound) {
                    continue; //skip record - none of its alternate alleles match the variant's alternate allele
                }
            }

            // only check this value for SNPs
            String hapRefValue = vc.isSNP() ? annotationsForRecord.get( generateInfoFieldKey(name, HAPLOTYPE_REFERENCE_COLUMN) ) : null;
            if(hapRefValue != null)
            {
                hapRefValue = hapRefValue.trim();
                if(!hapRefValue.equals("*"))
                {
                    //match against hapolotypeReference.
                    Allele vcRef = vc.getReference();
                    if(!vcRef.basesMatch(hapRefValue)) {
                        continue; //skip record
                    }
                }
            }

            if (vc.isIndel()) {
                modifyAnnotationsForIndels(vc, name, annotationsForRecord);
            }

            //filters passed, so add this record.
            List<Map<String, String>> listOfMatchingRecords = (List<Map<String, String>>) annotations.get( name );
            if(listOfMatchingRecords == null) {
                listOfMatchingRecords = new LinkedList<Map<String,String>>();
                listOfMatchingRecords.add( annotationsForRecord );
                annotations.put(name, listOfMatchingRecords);
            } else {
                listOfMatchingRecords.add( annotationsForRecord );
            }
        }

        return annotations;
    }




    /**
     * Converts the given record to a set of key-value pairs of the form:
     *   bindingName.columnName1=column1Value, bindingName.columnName2=column2Value
     *   (eg. dbSNP.avHet=0.7, dbSNP.ref_allele=A)
     *
     * @param record AnnotatorInputTableFeature corresponding to one record in one -B input file.
     * @param bindingName The binding name of the given AnnotatorInputTableFeature.
     * @return The map of columnName -> columnValue pairs.
     */
    public static Map<String, String> convertRecordToAnnotations( String bindingName, Map<String, String> record) {
        final Map<String, String> result = new HashMap<String, String>();

        for(final Entry<String, String> entry : record.entrySet()) {
            final String value = entry.getValue();
            if(!value.trim().isEmpty()) {
                result.put( generateInfoFieldKey(bindingName, entry.getKey()), scrubInfoFieldValue(entry.getValue()));
            }
        }

        return result;
    }

    /**
     * Combines the 2 values into a full key.
     * @param rodBindingName -B name
     * @param columnName     column name
     * @return info field key
     */
    public static String generateInfoFieldKey(String rodBindingName, String columnName ) {
        return rodBindingName + '.' + columnName;
    }



    /**
     * Replaces any characters that are not allowed in the info field of a VCF file.
     *
     * @param value info field value
     * @return the value with any illegal characters replaced by legal ones.
     */
    private static String scrubInfoFieldValue(String value) {
        for(int i = 0; i < GenomicAnnotation.ILLEGAL_INFO_FIELD_VALUES.length; i++) {
            value = value.replace(GenomicAnnotation.ILLEGAL_INFO_FIELD_VALUES[i], GenomicAnnotation.ILLEGAL_INFO_FIELD_VALUE_SUBSTITUTES[i]);
        }

        return value;
    }



    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine("GenericAnnotation", 1, VCFHeaderLineType.Integer, "For each variant in the 'variants' ROD, finds all entries in the other -B files that overlap the variant's position."));
    }

    public List<String> getKeyNames() {
        return Arrays.asList("GenericAnnotation");
    }

}
