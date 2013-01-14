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


import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * This class encompasses all the basic information about a genotype.  It is immutable.
 *
 * @author Mark DePristo
 */
public class Genotype {

    public final static String PHASED_ALLELE_SEPARATOR = "|";
    public final static String UNPHASED_ALLELE_SEPARATOR = "/";

    protected InferredGeneticContext commonInfo;
    public final static double NO_NEG_LOG_10PERROR = InferredGeneticContext.NO_NEG_LOG_10PERROR;
    protected List<Allele> alleles = null; // new ArrayList<Allele>();
    protected Type type = null;

    protected boolean isPhased = false;
    protected boolean filtersWereAppliedToContext;

    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, Set<String> filters, Map<String, ?> attributes, boolean isPhased) {
        this(sampleName, alleles, negLog10PError, filters, attributes, isPhased, null);
    }

    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, Set<String> filters, Map<String, ?> attributes, boolean isPhased, double[] log10Likelihoods) {
        if ( alleles != null )
            this.alleles = Collections.unmodifiableList(alleles);
        commonInfo = new InferredGeneticContext(sampleName, negLog10PError, filters, attributes);
        if ( log10Likelihoods != null )
            commonInfo.putAttribute(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods(log10Likelihoods));
        filtersWereAppliedToContext = filters != null;
        this.isPhased = isPhased;
        validate();
    }

    /**
     * Creates a new Genotype for sampleName with genotype according to alleles.
     * @param sampleName
     * @param alleles
     * @param negLog10PError the confidence in these alleles
     * @param log10Likelihoods a log10 likelihoods for each of the genotype combinations possible for alleles, in the standard VCF ordering, or null if not known
     */
    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, double[] log10Likelihoods) {
        this(sampleName, alleles, negLog10PError, null, null, false, log10Likelihoods);
    }

    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError) {
        this(sampleName, alleles, negLog10PError, null, null, false);
    }

    public Genotype(String sampleName, List<Allele> alleles) {
        this(sampleName, alleles, NO_NEG_LOG_10PERROR, null, null, false);
    }


    // ---------------------------------------------------------------------------------------------------------
    //
    // Partial-cloning routines (because Genotype is immutable).
    //
    // ---------------------------------------------------------------------------------------------------------

    public static Genotype modifyName(Genotype g, String name) {
        return new Genotype(name, g.getAlleles(), g.getNegLog10PError(), g.filtersWereApplied() ? g.getFilters() : null, g.getAttributes(), g.isPhased());
    }

    public static Genotype modifyAttributes(Genotype g, Map<String, Object> attributes) {
        return new Genotype(g.getSampleName(), g.getAlleles(), g.getNegLog10PError(), g.filtersWereApplied() ? g.getFilters() : null, attributes, g.isPhased());
    }

    public static Genotype modifyAlleles(Genotype g, List<Allele> alleles) {
        return new Genotype(g.getSampleName(), alleles, g.getNegLog10PError(), g.filtersWereApplied() ? g.getFilters() : null, g.getAttributes(), g.isPhased());
    }

    /**
     * @return the alleles for this genotype
     */
    public List<Allele> getAlleles() {
        return alleles;
    }

    public List<Allele> getAlleles(Allele allele) {
        if ( getType() == Type.UNAVAILABLE )
            throw new ReviewedStingException("Requesting alleles for an UNAVAILABLE genotype");

        List<Allele> al = new ArrayList<Allele>();
        for ( Allele a : alleles )
            if ( a.equals(allele) )
                al.add(a);

        return Collections.unmodifiableList(al);
    }

    public Allele getAllele(int i) {
        if ( getType() == Type.UNAVAILABLE )
            throw new ReviewedStingException("Requesting alleles for an UNAVAILABLE genotype");
        return alleles.get(i);
    }

    public boolean isPhased() { return isPhased; }

    /**
     * @return the ploidy of this genotype
     */
    public int getPloidy() {
        if ( alleles == null )
            throw new ReviewedStingException("Requesting ploidy for an UNAVAILABLE genotype");
        return alleles.size();
    }

    public enum Type {
        NO_CALL,
        HOM_REF,
        HET,
        HOM_VAR,
        UNAVAILABLE,
        MIXED  // no-call and call in the same genotype
    }

    public Type getType() {
        if ( type == null ) {
            type = determineType();
        }
        return type;
    }

    protected Type determineType() {
        if ( alleles == null )
            return Type.UNAVAILABLE;

        boolean sawNoCall = false, sawMultipleAlleles = false;
        Allele observedAllele = null;

        for ( Allele allele : alleles ) {
            if ( allele.isNoCall() )
                sawNoCall = true;
            else if ( observedAllele == null )
                observedAllele = allele;
            else if ( !allele.equals(observedAllele) )
                sawMultipleAlleles = true;
        }

        if ( sawNoCall ) {
            if ( observedAllele == null )
                return Type.NO_CALL;
            return Type.MIXED;
        }

        if ( observedAllele == null )
            throw new ReviewedStingException("BUG: there are no alleles present in this genotype but the alleles list is not null");

        return sawMultipleAlleles ? Type.HET : observedAllele.isReference() ? Type.HOM_REF : Type.HOM_VAR;
    }

    /**
     * @return true if all observed alleles are the same (regardless of whether they are ref or alt); if any alleles are no-calls, this method will return false.
     */
    public boolean isHom()    { return isHomRef() || isHomVar(); }

    /**
     * @return true if all observed alleles are ref; if any alleles are no-calls, this method will return false.
     */
    public boolean isHomRef() { return getType() == Type.HOM_REF; }

    /**
     * @return true if all observed alleles are alt; if any alleles are no-calls, this method will return false.
     */
    public boolean isHomVar() { return getType() == Type.HOM_VAR; }
    
    /**
     * @return true if we're het (observed alleles differ); if the ploidy is less than 2 or if any alleles are no-calls, this method will return false.
     */
    public boolean isHet() { return getType() == Type.HET; }

    /**
     * @return true if this genotype is not actually a genotype but a "no call" (e.g. './.' in VCF); if any alleles are not no-calls (even if some are), this method will return false.
     */
    public boolean isNoCall() { return getType() == Type.NO_CALL; }

    /**
     * @return true if this genotype is comprised of any alleles that are not no-calls (even if some are).
     */
    public boolean isCalled() { return getType() != Type.NO_CALL && getType() != Type.UNAVAILABLE; }

    /**
     * @return true if this genotype is comprised of both calls and no-calls.
     */
    public boolean isMixed() { return getType() == Type.MIXED; }

    /**
     * @return true if the type of this genotype is set.
     */
    public boolean isAvailable() { return getType() != Type.UNAVAILABLE; }

    //
    // Useful methods for getting genotype likelihoods for a genotype object, if present
    //
    public boolean hasLikelihoods() {
        return (hasAttribute(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY) && !getAttribute(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY).equals(VCFConstants.MISSING_VALUE_v4)) ||
                (hasAttribute(VCFConstants.GENOTYPE_LIKELIHOODS_KEY) && !getAttribute(VCFConstants.GENOTYPE_LIKELIHOODS_KEY).equals(VCFConstants.MISSING_VALUE_v4));
    }
    
    public GenotypeLikelihoods getLikelihoods() {
        GenotypeLikelihoods x = getLikelihoods(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, true);
        if ( x != null )
            return x;
        else {
            x = getLikelihoods(VCFConstants.GENOTYPE_LIKELIHOODS_KEY, false);
            if ( x != null ) return x;
            else
                throw new IllegalStateException("BUG: genotype likelihood field in " + this.getSampleName() + " sample are not either a string or a genotype likelihood class!");
        }
    }

    private GenotypeLikelihoods getLikelihoods(String key, boolean asPL) {
        Object x = getAttribute(key);
        if ( x instanceof String ) {
            if ( asPL )
                return GenotypeLikelihoods.fromPLField((String)x);
            else
                return GenotypeLikelihoods.fromGLField((String)x);
        }
        else if ( x instanceof GenotypeLikelihoods ) return (GenotypeLikelihoods)x;
        else return null;
    }

    public void validate() {
        if ( alleles == null ) return;
        if ( alleles.size() == 0) throw new IllegalArgumentException("BUG: alleles cannot be of size 0");

        // int nNoCalls = 0;
        for ( Allele allele : alleles ) {
            if ( allele == null )
                throw new IllegalArgumentException("BUG: allele cannot be null in Genotype");
            // nNoCalls += allele.isNoCall() ? 1 : 0;
        }

        // Technically, the spec does allow for the below case so this is not an illegal state
        //if ( nNoCalls > 0 && nNoCalls != alleles.size() )
        //    throw new IllegalArgumentException("BUG: alleles include some No Calls and some Calls, an illegal state " + this);
    }

    public String getGenotypeString() {
        return getGenotypeString(true);
    }

    public String getGenotypeString(boolean ignoreRefState) {
        if ( alleles == null )
            return null;

        // Notes:
        // 1. Make sure to use the appropriate separator depending on whether the genotype is phased
        // 2. If ignoreRefState is true, then we want just the bases of the Alleles (ignoring the '*' indicating a ref Allele)
        // 3. So that everything is deterministic with regards to integration tests, we sort Alleles (when the genotype isn't phased, of course)
        return ParsingUtils.join(isPhased() ? PHASED_ALLELE_SEPARATOR : UNPHASED_ALLELE_SEPARATOR,
                ignoreRefState ? getAlleleStrings() : (isPhased() ? getAlleles() : ParsingUtils.sortList(getAlleles())));
    }

    private List<String> getAlleleStrings() {
        List<String> al = new ArrayList<String>();
        for ( Allele a : alleles )
            al.add(a.getBaseString());

        return al;
    }

    public String toString() {
        int Q = (int)Math.round(getPhredScaledQual());
        return String.format("[%s %s Q%s %s]", getSampleName(), getGenotypeString(false),
                Q == -10 ? "." : String.format("%2d",Q), sortedString(getAttributes()));
    }

    public String toBriefString() {
        return String.format("%s:Q%.2f", getGenotypeString(false), getPhredScaledQual());
    }

    public boolean sameGenotype(Genotype other) {
        return sameGenotype(other, true);
    }

    public boolean sameGenotype(Genotype other, boolean ignorePhase) {
        if (getPloidy() != other.getPloidy())
            return false; // gotta have the same number of allele to be equal

        // By default, compare the elements in the lists of alleles, element-by-element
        Collection<Allele> thisAlleles = this.getAlleles();
        Collection<Allele> otherAlleles = other.getAlleles();

        if (ignorePhase) { // do not care about order, only identity of Alleles
            thisAlleles = new TreeSet<Allele>(thisAlleles);   //implemented Allele.compareTo()
            otherAlleles = new TreeSet<Allele>(otherAlleles);
        }

        return thisAlleles.equals(otherAlleles);
    }

    /**
     * a utility method for generating sorted strings from a map key set.
     * @param c the map
     * @param <T> the key type
     * @param <V> the value type
     * @return a sting, enclosed in {}, with comma seperated key value pairs in order of the keys
     */
    private static <T extends Comparable<T>, V> String sortedString(Map<T, V> c) {
        // NOTE -- THIS IS COPIED FROM GATK UTILS TO ALLOW US TO KEEP A SEPARATION BETWEEN THE GATK AND VCF CODECS
        List<T> t = new ArrayList<T>(c.keySet());
        Collections.sort(t);

        List<String> pairs = new ArrayList<String>();
        for (T k : t) {
            pairs.add(k + "=" + c.get(k));
        }

        return "{" + ParsingUtils.join(", ", pairs.toArray(new String[pairs.size()])) + "}";
    }


    // ---------------------------------------------------------------------------------------------------------
    // 
    // get routines to access context info fields
    //
    // ---------------------------------------------------------------------------------------------------------
    public String getSampleName()       { return commonInfo.getName(); }
    public Set<String> getFilters()     { return commonInfo.getFilters(); }
    public boolean isFiltered()         { return commonInfo.isFiltered(); }
    public boolean isNotFiltered()      { return commonInfo.isNotFiltered(); }
    public boolean filtersWereApplied() { return filtersWereAppliedToContext; }
    public boolean hasNegLog10PError()  { return commonInfo.hasNegLog10PError(); }
    public double getNegLog10PError()   { return commonInfo.getNegLog10PError(); }
    public double getPhredScaledQual()  { return commonInfo.getPhredScaledQual(); }

    public Map<String, Object> getAttributes()  { return commonInfo.getAttributes(); }
    public boolean hasAttribute(String key)     { return commonInfo.hasAttribute(key); }
    public Object getAttribute(String key)      { return commonInfo.getAttribute(key); }
    
    public Object getAttribute(String key, Object defaultValue) {
        return commonInfo.getAttribute(key, defaultValue); 
    }

    public String getAttributeAsString(String key, String defaultValue)   { return commonInfo.getAttributeAsString(key, defaultValue); }
    public int getAttributeAsInt(String key, int defaultValue)            { return commonInfo.getAttributeAsInt(key, defaultValue); }
    public double getAttributeAsDouble(String key, double  defaultValue)  { return commonInfo.getAttributeAsDouble(key, defaultValue); }
    public boolean getAttributeAsBoolean(String key, boolean  defaultValue)  { return commonInfo.getAttributeAsBoolean(key, defaultValue); }
}