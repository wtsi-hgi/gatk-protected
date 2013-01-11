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

package org.broadinstitute.variant.variantcontext.v13;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * a base class for compound header lines, which include info lines and format lines (so far)
 */
abstract class VCFCompoundHeaderLine extends VCFHeaderLine implements VCFNamedHeaderLine {
    public enum SupportedHeaderLineType {
        INFO(true), FORMAT(false);

        public final boolean allowFlagValues;
        SupportedHeaderLineType(boolean flagValues) {
            allowFlagValues = flagValues;
        }
    }

    // the field types
    private String name;
    private int count = -1;
    private VCFHeaderLineCount countType;
    private String description;
    private VCFHeaderLineType type;

    // access methods
    public String getName() { return name; }
    public String getDescription() { return description; }
    public VCFHeaderLineType getType() { return type; }
    public VCFHeaderLineCount getCountType() { return countType; }
    public int getCount() {
        if ( countType != VCFHeaderLineCount.INTEGER )
            throw new ReviewedStingException("Asking for header line count when type is not an integer");
        return count;
    }

    // utility method
    public int getCount(int numAltAlleles) {
        int myCount;
        switch ( countType ) {
            case INTEGER: myCount = count; break;
            case UNBOUNDED: myCount = -1; break;
            case A: myCount = numAltAlleles; break;
            case G: myCount = ((numAltAlleles + 1) * (numAltAlleles + 2) / 2); break;
            default: throw new ReviewedStingException("Unknown count type: " + countType);
        }
        return myCount;
    }

    public void setNumberToUnbounded() {
        countType = VCFHeaderLineCount.UNBOUNDED;
        count = -1;
    }

    // our type of line, i.e. format, info, etc
    private final SupportedHeaderLineType lineType;

    /**
     * create a VCF format header line
     *
     * @param name         the name for this header line
     * @param count        the count for this header line
     * @param type         the type for this header line
     * @param description  the description for this header line
     * @param lineType     the header line type
     */
    protected VCFCompoundHeaderLine(String name, int count, VCFHeaderLineType type, String description, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "");
        this.name = name;
        this.countType = VCFHeaderLineCount.INTEGER;
        this.count = count;
        this.type = type;
        this.description = description;
        this.lineType = lineType;
        validate();
    }

    /**
     * create a VCF format header line
     *
     * @param name         the name for this header line
     * @param count        the count type for this header line
     * @param type         the type for this header line
     * @param description  the description for this header line
     * @param lineType     the header line type
     */
    protected VCFCompoundHeaderLine(String name, VCFHeaderLineCount count, VCFHeaderLineType type, String description, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "");
        this.name = name;
        this.countType = count;
        this.type = type;
        this.description = description;
        this.lineType = lineType;
        validate();
    }

    /**
     * create a VCF format header line
     *
     * @param line   the header line
     * @param version      the VCF header version
     * @param lineType     the header line type
     *
     */
    protected VCFCompoundHeaderLine(String line, VCFHeaderVersion version, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "");
        Map<String,String> mapping = VCFHeaderLineTranslator.parseLine(version,line, Arrays.asList("ID","Number","Type","Description"));
        name = mapping.get("ID");
        count = -1;
        final String numberStr = mapping.get("Number");
        if ( numberStr.equals(VCFConstants.PER_ALLELE_COUNT) ) {
            countType = VCFHeaderLineCount.A;
        } else if ( numberStr.equals(VCFConstants.PER_GENOTYPE_COUNT) ) {
            countType = VCFHeaderLineCount.G;
        } else if ( ((version == VCFHeaderVersion.VCF4_0 || version == VCFHeaderVersion.VCF4_1) &&
                     numberStr.equals(VCFConstants.UNBOUNDED_ENCODING_v4)) ||
                    ((version == VCFHeaderVersion.VCF3_2 || version == VCFHeaderVersion.VCF3_3) &&
                     numberStr.equals(VCFConstants.UNBOUNDED_ENCODING_v3)) ) {
            countType = VCFHeaderLineCount.UNBOUNDED;
        } else {
            countType = VCFHeaderLineCount.INTEGER;
            count = Integer.valueOf(numberStr);

        }
        type = VCFHeaderLineType.valueOf(mapping.get("Type"));
        if (type == VCFHeaderLineType.Flag && !allowFlagValues())
            throw new IllegalArgumentException("Flag is an unsupported type for this kind of field");

        description = mapping.get("Description");
        if ( description == null && ALLOW_UNBOUND_DESCRIPTIONS ) // handle the case where there's no description provided
            description = UNBOUND_DESCRIPTION;
        
        this.lineType = lineType;

        validate();
    }

    private void validate() {
        if ( name == null || type == null || description == null || lineType == null )
            throw new IllegalArgumentException(String.format("Invalid VCFCompoundHeaderLine: key=%s name=%s type=%s desc=%s lineType=%s", 
                    super.getKey(), name, type, description, lineType ));
    }

    /**
     * make a string representation of this header line
     * @return a string representation
     */
    protected String toStringEncoding() {
        Map<String,Object> map = new LinkedHashMap<String,Object>();
        map.put("ID", name);
        Object number;
        switch ( countType ) {
            case A: number = VCFConstants.PER_ALLELE_COUNT; break;
            case G: number = VCFConstants.PER_GENOTYPE_COUNT; break;
            case UNBOUNDED: number = VCFConstants.UNBOUNDED_ENCODING_v4; break;
            case INTEGER:
            default: number = count;
        }
        map.put("Number", number);
        map.put("Type", type);
        map.put("Description", description);
        return lineType.toString() + "=" + toStringEncoding(map);
    }

    /**
     * returns true if we're equal to another compounder header line
     * @param o a compound header line
     * @return true if equal
     */
    public boolean equals(Object o) {
        if ( !(o instanceof VCFCompoundHeaderLine) )
            return false;
        VCFCompoundHeaderLine other = (VCFCompoundHeaderLine)o;
        return equalsExcludingDescription(other) &&
                description.equals(other.description);
    }

    public boolean equalsExcludingDescription(VCFCompoundHeaderLine other) {
        return count == other.count &&
                countType == other.countType &&
                type == other.type &&
                lineType == other.lineType &&
                name.equals(other.name);
    }

    public boolean sameLineTypeAndName(VCFCompoundHeaderLine other) {
        return  lineType == other.lineType &&
                name.equals(other.name);
    }

    /**
     * do we allow flag (boolean) values? (i.e. booleans where you don't have specify the value, AQ means AQ=true)
     * @return true if we do, false otherwise
     */
    abstract boolean allowFlagValues();

}
