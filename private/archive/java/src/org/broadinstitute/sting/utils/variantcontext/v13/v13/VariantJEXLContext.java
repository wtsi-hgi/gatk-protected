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

import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.MapContext;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author aaron
 * @author depristo
 *
 * Class VariantJEXLContext
 *
 * implements the JEXML context for VariantContext; this saves us from
 * having to generate a JEXML context lookup map everytime we want to evaluate an expression.
 *
 * This is package protected, only classes in variantcontext should have access to it.
 *
 * // todo -- clean up to remove or better support genotype filtering 
 */

class VariantJEXLContext implements JexlContext {
    // our stored variant context
    private VariantContext vc;

    private interface AttributeGetter {
        public Object get(VariantContext vc);
    }

    private static Map<String, AttributeGetter> x = new HashMap<String, AttributeGetter>();

    static {
        x.put("vc",   new AttributeGetter() { public Object get(VariantContext vc) { return vc; }});
        x.put("CHROM",   new AttributeGetter() { public Object get(VariantContext vc) { return vc.getChr(); }});
        x.put("POS",     new AttributeGetter() { public Object get(VariantContext vc) { return vc.getStart(); }});
        x.put("TYPE",    new AttributeGetter() { public Object get(VariantContext vc) { return vc.getType().toString(); }});
        x.put("QUAL",    new AttributeGetter() { public Object get(VariantContext vc) { return 10 * vc.getNegLog10PError(); }});
        x.put("ALLELES", new AttributeGetter() { public Object get(VariantContext vc) { return vc.getAlleles(); }});
        x.put("N_ALLELES", new AttributeGetter() { public Object get(VariantContext vc) { return vc.getNAlleles(); }});
        x.put("FILTER",    new AttributeGetter() { public Object get(VariantContext vc) { return vc.isFiltered() ? "1" : "0"; }});

//        x.put("GT",        new AttributeGetter() { public Object get(VariantContext vc) { return g.getGenotypeString(); }});
        x.put("homRefCount",  new AttributeGetter() { public Object get(VariantContext vc) { return vc.getHomRefCount(); }});
        x.put("hetCount",     new AttributeGetter() { public Object get(VariantContext vc) { return vc.getHetCount(); }});
        x.put("homVarCount",  new AttributeGetter() { public Object get(VariantContext vc) { return vc.getHomVarCount(); }});
    }

    public VariantJEXLContext(VariantContext vc) {
        this.vc = vc;
    }

    public Object get(String name) {
        Object result = null;
        if ( x.containsKey(name) ) { // dynamic resolution of name -> value via map
            result = x.get(name).get(vc);
        } else if ( vc.hasAttribute(name)) {
            result = vc.getAttribute(name);
        } else if ( vc.getFilters().contains(name) ) {
            result = "1";
        }

        //System.out.printf("dynamic lookup %s => %s%n", name, result);

        return result;
    }

    public boolean has(String name) {
        return get(name) != null;
    }

    public void	set(String name, Object value) {
        throw new UnsupportedOperationException("remove() not supported on a VariantJEXLContext");
    }
}




/**
 * this is an implementation of a Map of JexlVCMatchExp to true or false values.  It lazy initializes each value
 * as requested to save as much processing time as possible.
 *
 * Compatible with JEXL 1.1 (this code will be easier if we move to 2.0, all of the functionality can go into the
 * JexlContext's get()
 * 
 */

class JEXLMap implements Map<VariantContextUtils.JexlVCMatchExp, Boolean> {
    // our variant context and/or Genotype
    private final VariantContext vc;
    private final Genotype g;

    // our context
    private JexlContext jContext = null;

    // our mapping from JEXLVCMatchExp to Booleans, which will be set to NULL for previously uncached JexlVCMatchExp
    private Map<VariantContextUtils.JexlVCMatchExp,Boolean> jexl;


    public JEXLMap(Collection<VariantContextUtils.JexlVCMatchExp> jexlCollection, VariantContext vc, Genotype g) {
        this.vc = vc;
        this.g = g;
        initialize(jexlCollection);
    }

    public JEXLMap(Collection<VariantContextUtils.JexlVCMatchExp> jexlCollection, VariantContext vc) {
        this(jexlCollection, vc, null);
    }

    private void initialize(Collection<VariantContextUtils.JexlVCMatchExp> jexlCollection) {
        jexl = new HashMap<VariantContextUtils.JexlVCMatchExp,Boolean>();
        for (VariantContextUtils.JexlVCMatchExp exp: jexlCollection) {
            jexl.put(exp, null);
        }
    }

    /**
     * create the internal JexlContext, only when required.  This code is where new JEXL context variables
     * should get added.
     *
     */
    private void createContext() {
        if ( g == null ) {
            // todo -- remove dependancy on g to the entire system
            jContext = new VariantJEXLContext(vc);
        } else {
            //
            // this whole branch is here just to support G jexl operations
            //
            Map<String, Object> infoMap = new HashMap<String, Object>();

            if ( vc != null ) {
                // create a mapping of what we know about the variant context, its Chromosome, positions, etc.
                infoMap.put("CHROM", vc.getChr());
                infoMap.put("POS", vc.getStart());
                infoMap.put("TYPE", vc.getType().toString());
                infoMap.put("QUAL", String.valueOf(vc.getPhredScaledQual()));

                // add alleles
                infoMap.put("ALLELES", Utils.join(";", vc.getAlleles()));
                infoMap.put("N_ALLELES", String.valueOf(vc.getNAlleles()));

                // add attributes
                addAttributesToMap(infoMap, vc.getAttributes());

                // add filter fields
                infoMap.put("FILTER", vc.isFiltered() ? "1" : "0");
                for ( Object filterCode : vc.getFilters() ) {
                    infoMap.put(String.valueOf(filterCode), "1");
                }

                // add genotype-specific fields
                // TODO -- implement me when we figure out a good way to represent this
                //    for ( Genotype g : vc.getGenotypes().values() ) {
                //        String prefix = g.getSampleName() + ".";
                //        addAttributesToMap(infoMap, g.getAttributes(), prefix);
                //        infoMap.put(prefix + "GT", g.getGenotypeString());
                //    }

                // add specific genotype if one is provided
                infoMap.put(VCFConstants.GENOTYPE_KEY, g.getGenotypeString());
                infoMap.put("isHomRef", g.isHomRef() ? "1" : "0");
                infoMap.put("isHet", g.isHet() ? "1" : "0");
                infoMap.put("isHomVar", g.isHomVar() ? "1" : "0");
                infoMap.put(VCFConstants.GENOTYPE_QUALITY_KEY, new Double(g.getPhredScaledQual()));
                for ( Map.Entry<String, Object> e : g.getAttributes().entrySet() ) {
                    if ( e.getValue() != null && !e.getValue().equals(VCFConstants.MISSING_VALUE_v4) )
                        infoMap.put(e.getKey(), e.getValue());
                }
            }

            // create the internal context that we can evaluate expressions against

            jContext = new MapContext(infoMap);
        }
    }

    /**
     * @return the size of the internal data structure
     */
    public int size() {
        return jexl.size();
    }

    /**
     * @return true if we're empty
     */
    public boolean isEmpty() { return this.jexl.isEmpty(); }

    /**
     * do we contain the specified key
     * @param o the key
     * @return true if we have a value for that key
     */
    public boolean containsKey(Object o) { return jexl.containsKey(o); }

    public Boolean get(Object o) {
        // if we've already determined the value, return it
        if (jexl.containsKey(o) && jexl.get(o) != null) return jexl.get(o);

        // try and cast the expression
        VariantContextUtils.JexlVCMatchExp e = (VariantContextUtils.JexlVCMatchExp) o;
        evaluateExpression(e);
        return jexl.get(e);
    }

    /**
     * get the keyset of map
     * @return a set of keys of type JexlVCMatchExp
     */
    public Set<VariantContextUtils.JexlVCMatchExp> keySet() {
        return jexl.keySet();
    }

    /**
     * get all the values of the map.  This is an expensive call, since it evaluates all keys that haven't
     * been evaluated yet.  This is fine if you truely want all the keys, but if you only want a portion, or  know
     * the keys you want, you would be better off using get() to get them by name.
     * @return a collection of boolean values, representing the results of all the variants evaluated
     */
    public Collection<Boolean> values() {
        // this is an expensive call
        for (VariantContextUtils.JexlVCMatchExp exp : jexl.keySet())
            if (jexl.get(exp) == null)
                evaluateExpression(exp);
        return jexl.values();
    }

    /**
     * evaulate a JexlVCMatchExp's expression, given the current context (and setup the context if it's null)
     * @param exp the JexlVCMatchExp to evaluate
     */
    private void evaluateExpression(VariantContextUtils.JexlVCMatchExp exp) {
        // if the context is null, we need to create it to evaluate the JEXL expression
        if (this.jContext == null) createContext();
        try {
            jexl.put (exp, (Boolean) exp.exp.evaluate(jContext));
        } catch (Exception e) {
            throw new UserException.CommandLineException(String.format("Invalid JEXL expression detected for %s with message %s", exp.name, e.getMessage()));
        }
    }

    /**
     * helper function: adds the list of attributes to the information map we're building
     * @param infoMap the map
     * @param attributes the attributes
     */
    private static void addAttributesToMap(Map<String, Object> infoMap, Map<String, ?> attributes ) {
        for (Map.Entry<String, ?> e : attributes.entrySet()) {
            infoMap.put(e.getKey(), String.valueOf(e.getValue()));
        }
    }

    public Boolean put(VariantContextUtils.JexlVCMatchExp jexlVCMatchExp, Boolean aBoolean) {
        return jexl.put(jexlVCMatchExp,aBoolean);
    }

    public void putAll(Map<? extends VariantContextUtils.JexlVCMatchExp, ? extends Boolean> map) {
        jexl.putAll(map);
    }

    // //////////////////////////////////////////////////////////////////////////////////////
    // The Following are unsupported at the moment
    // //////////////////////////////////////////////////////////////////////////////////////

    // this doesn't make much sense to implement, boolean doesn't offer too much variety to deal
    // with evaluating every key in the internal map.
    public boolean containsValue(Object o) {
        throw new UnsupportedOperationException("containsValue() not supported on a JEXLMap");
    }

    // this doesn't make much sense
    public Boolean remove(Object o) {
        throw new UnsupportedOperationException("remove() not supported on a JEXLMap");
    }


    public Set<Entry<VariantContextUtils.JexlVCMatchExp, Boolean>> entrySet() {
        throw new UnsupportedOperationException("clear() not supported on a JEXLMap");
    }

    // nope
    public void clear() {
        throw new UnsupportedOperationException("clear() not supported on a JEXLMap");
    }
}
