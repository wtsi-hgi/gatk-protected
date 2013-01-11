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


import java.util.*;


/**
 * Common utility routines for VariantContext and Genotype
 *
 * @author depristo
 */
final class InferredGeneticContext {
    public static final double NO_NEG_LOG_10PERROR = -1.0;

    private static Set<String> NO_FILTERS = Collections.unmodifiableSet(new HashSet<String>());
    private static Map<String, Object> NO_ATTRIBUTES = Collections.unmodifiableMap(new HashMap<String, Object>());

    private double negLog10PError = NO_NEG_LOG_10PERROR;
    private String name = null;
    private Set<String> filters = NO_FILTERS;
    private Map<String, Object> attributes = NO_ATTRIBUTES;

//    public InferredGeneticContext(String name) {
//        this.name = name;
//    }
//
//    public InferredGeneticContext(String name, double negLog10PError) {
//        this(name);
//        setNegLog10PError(negLog10PError);
//    }

    public InferredGeneticContext(String name, double negLog10PError, Set<String> filters, Map<String, ?> attributes) {
        this.name = name;
        setNegLog10PError(negLog10PError);
        if ( filters != null )
            setFilters(filters);
        if ( attributes != null )
            setAttributes(attributes);
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * Sets the name
     *
     * @param name    the name associated with this information
     */
    public void setName(String name) {
        if ( name == null ) throw new IllegalArgumentException("Name cannot be null " + this);
        this.name = name;
    }


    // ---------------------------------------------------------------------------------------------------------
    //
    // Filter
    //
    // ---------------------------------------------------------------------------------------------------------

    public Set<String> getFilters() {
        return Collections.unmodifiableSet(filters);
    }

    public boolean isFiltered() {
        return filters.size() > 0;
    }

    public boolean isNotFiltered() {
        return ! isFiltered();
    }

    public void addFilter(String filter) {
        if ( filters == NO_FILTERS ) // immutable -> mutable
            filters = new HashSet<String>(filters);

        if ( filter == null ) throw new IllegalArgumentException("BUG: Attempting to add null filter " + this);
        if ( getFilters().contains(filter) ) throw new IllegalArgumentException("BUG: Attempting to add duplicate filter " + filter + " at " + this);
        filters.add(filter);
    }

    public void addFilters(Collection<String> filters) {
        if ( filters == null ) throw new IllegalArgumentException("BUG: Attempting to add null filters at" + this);
        for ( String f : filters )
            addFilter(f);
    }

    public void clearFilters() {
        filters = new HashSet<String>();
    }

    public void setFilters(Collection<String> filters) {
        clearFilters();
        addFilters(filters);
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with log error rates
    //
    // ---------------------------------------------------------------------------------------------------------

    public boolean hasNegLog10PError() {
        return getNegLog10PError() != NO_NEG_LOG_10PERROR;
    }

    /**
     * @return the -1 * log10-based error estimate
     */
    public double getNegLog10PError() { return negLog10PError; }
    public double getPhredScaledQual() { return getNegLog10PError() * 10; }

    public void setNegLog10PError(double negLog10PError) {
        if ( negLog10PError < 0 && negLog10PError != NO_NEG_LOG_10PERROR ) throw new IllegalArgumentException("BUG: negLog10PError cannot be < than 0 : " + negLog10PError);
        if ( Double.isInfinite(negLog10PError) ) throw new IllegalArgumentException("BUG: negLog10PError should not be Infinity");
        if ( Double.isNaN(negLog10PError) ) throw new IllegalArgumentException("BUG: negLog10PError should not be NaN");

        this.negLog10PError = negLog10PError;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with attributes
    //
    // ---------------------------------------------------------------------------------------------------------
    public void clearAttributes() {
        attributes = new HashMap<String, Object>();
    }

    /**
     * @return the attribute map
     */
    public Map<String, Object> getAttributes() {
        return Collections.unmodifiableMap(attributes);
    }

    // todo -- define common attributes as enum

    public void setAttributes(Map<String, ?> map) {
        clearAttributes();
        putAttributes(map);
    }

    public void putAttribute(String key, Object value) {
        putAttribute(key, value, false);
    }

    public void putAttribute(String key, Object value, boolean allowOverwrites) {
        if ( ! allowOverwrites && hasAttribute(key) )
            throw new IllegalStateException("Attempting to overwrite key->value binding: key = " + key + " this = " + this);

        if ( attributes == NO_ATTRIBUTES ) // immutable -> mutable
            attributes = new HashMap<String, Object>();
        
        attributes.put(key, value);
    }

    public void removeAttribute(String key) {
        if ( attributes == NO_ATTRIBUTES ) // immutable -> mutable
            attributes = new HashMap<String, Object>();
        attributes.remove(key);
    }

    public void putAttributes(Map<String, ?> map) {
        if ( map != null ) {
            // for efficiency, we can skip the validation if the map is empty
            if ( attributes.size() == 0 ) {
                if ( attributes == NO_ATTRIBUTES ) // immutable -> mutable
                    attributes = new HashMap<String, Object>();
                attributes.putAll(map);
            } else {
                for ( Map.Entry<String, ?> elt : map.entrySet() ) {
                    putAttribute(elt.getKey(), elt.getValue(), false);
                }
            }
        }
    }

    public boolean hasAttribute(String key) {
        return attributes.containsKey(key);
    }

    public int getNumAttributes() {
        return attributes.size();
    }

    /**
     * @param key    the attribute key
     *
     * @return the attribute value for the given key (or null if not set)
     */
    public Object getAttribute(String key) {
        return attributes.get(key);
    }

    public Object getAttribute(String key, Object defaultValue) {
        if ( hasAttribute(key) )
            return attributes.get(key);
        else
            return defaultValue;
    }

    public String getAttributeAsString(String key, String defaultValue) {
        Object x = getAttribute(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof String ) return (String)x;
        return String.valueOf(x); // throws an exception if this isn't a string
    }

    public int getAttributeAsInt(String key, int defaultValue) {
        Object x = getAttribute(key);
        if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) return defaultValue;
        if ( x instanceof Integer ) return (Integer)x;
        return Integer.valueOf((String)x); // throws an exception if this isn't a string
    }

    public double getAttributeAsDouble(String key, double defaultValue) {
        Object x = getAttribute(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof Double ) return (Double)x;
        return Double.valueOf((String)x); // throws an exception if this isn't a string
    }

    public boolean getAttributeAsBoolean(String key, boolean defaultValue) {
        Object x = getAttribute(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof Boolean ) return (Boolean)x;
        return Boolean.valueOf((String)x); // throws an exception if this isn't a string
    }

//    public String getAttributeAsString(String key)      { return (String.valueOf(getAttribute(key))); } // **NOTE**: will turn a null Object into the String "null"
//    public int getAttributeAsInt(String key)            { Object x = getAttribute(key); return x instanceof Integer ? (Integer)x : Integer.valueOf((String)x); }
//    public double getAttributeAsDouble(String key)      { Object x = getAttribute(key); return x instanceof Double ? (Double)x : Double.valueOf((String)x); }
//    public boolean getAttributeAsBoolean(String key)      { Object x = getAttribute(key); return x instanceof Boolean ? (Boolean)x : Boolean.valueOf((String)x); }
//    public Integer getAttributeAsIntegerNoException(String key)  { try {return getAttributeAsInt(key);} catch (Exception e) {return null;} }
//    public Double getAttributeAsDoubleNoException(String key)    { try {return getAttributeAsDouble(key);} catch (Exception e) {return null;} }
//    public String getAttributeAsStringNoException(String key)    { if (getAttribute(key) == null) return null; return getAttributeAsString(key); }
//    public Boolean getAttributeAsBooleanNoException(String key)  { try {return getAttributeAsBoolean(key);} catch (Exception e) {return null;} }
}