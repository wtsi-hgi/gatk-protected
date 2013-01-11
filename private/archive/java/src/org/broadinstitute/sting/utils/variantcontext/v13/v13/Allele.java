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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Immutable representation of an allele
 *
 * Types of alleles:
 *
 * Ref: a t C g a // C is the reference base
 *
 *    : a t G g a // C base is a G in some individuals
 *
 *    : a t - g a // C base is deleted w.r.t. the reference
 *
 *    : a t CAg a // A base is inserted w.r.t. the reference sequence
 *
 * In these cases, where are the alleles?
 *
 * SNP polymorphism of C/G  -> { C , G } -> C is the reference allele
 * 1 base deletion of C     -> { C , - } -> C is the reference allele
 * 1 base insertion of A    -> { - ; A } -> Null is the reference allele
 *
 * Suppose I see a the following in the population:
 *
 * Ref: a t C g a // C is the reference base
 *    : a t G g a // C base is a G in some individuals
 *    : a t - g a // C base is deleted w.r.t. the reference
 *
 * How do I represent this?  There are three segregating alleles:
 *
 *  { C , G , - }
 *
 * Now suppose I have this more complex example:
 *
 * Ref: a t C g a // C is the reference base
 *    : a t - g a
 *    : a t - - a
 *    : a t CAg a
 *
 * There are actually four segregating alleles:
 *
 *   { C g , - g, - -, and CAg } over bases 2-4
 *
 * However, the molecular equivalence explicitly listed above is usually discarded, so the actual
 * segregating alleles are:
 *
 *   { C g, g, -, C a g }
 *
 * Critically, it should be possible to apply an allele to a reference sequence to create the
 * correct haplotype sequence:
 *
 * Allele + reference => haplotype
 *
 * For convenience, we are going to create Alleles where the GenomeLoc of the allele is stored outside of the
 * Allele object itself.  So there's an idea of an A/C polymorphism independent of it's surrounding context.
 *
 * Given list of alleles it's possible to determine the "type" of the variation
 *
 *      A / C @ loc => SNP with
 *      - / A => INDEL
 *
 * If you know where allele is the reference, you can determine whether the variant is an insertion or deletion.
 *
 * Alelle also supports is concept of a NO_CALL allele.  This Allele represents a haplotype that couldn't be
 * determined. This is usually represented by a '.' allele.
 *
 * Note that Alleles store all bases as bytes, in **UPPER CASE**.  So 'atc' == 'ATC' from the perspective of an
 * Allele.

 * @author ebanks, depristo
 */
class Allele implements Comparable<Allele> {
    private static final byte[] EMPTY_ALLELE_BASES = new byte[0];

    private boolean isRef = false;
    private boolean isNull = false;
    private boolean isNoCall = false;
    private boolean isSymbolic = false;

    private byte[] bases = null;

    public final static String NULL_ALLELE_STRING = "-";
    public final static String NO_CALL_STRING = ".";
    /** A generic static NO_CALL allele for use */

    // no public way to create an allele
    private Allele(byte[] bases, boolean isRef) {
        // standardize our representation of null allele and bases
        if ( wouldBeNullAllele(bases) ) {
            bases = EMPTY_ALLELE_BASES;
            isNull = true;
        } else if ( wouldBeNoCallAllele(bases) ) {
            bases = EMPTY_ALLELE_BASES;
            isNoCall = true;
            if ( isRef ) throw new IllegalArgumentException("Cannot tag a NoCall allele as the reference allele");
        } else if ( wouldBeSymbolicAllele(bases) ) {
            isSymbolic = true;
            if ( isRef ) throw new IllegalArgumentException("Cannot tag a symbolic allele as the reference allele");
        }
//        else
//            bases = new String(bases).toUpperCase().getBytes(); // todo -- slow performance

        this.isRef = isRef;
        this.bases = bases;

        if ( ! acceptableAlleleBases(bases) )
            throw new IllegalArgumentException("Unexpected base in allele bases \'" + new String(bases)+"\'");
    }

    private Allele(String bases, boolean isRef) {
        this(bases.getBytes(), isRef);
    }


    private final static Allele REF_A = new Allele("A", true);
    private final static Allele ALT_A = new Allele("A", false);
    private final static Allele REF_C = new Allele("C", true);
    private final static Allele ALT_C = new Allele("C", false);
    private final static Allele REF_G = new Allele("G", true);
    private final static Allele ALT_G = new Allele("G", false);
    private final static Allele REF_T = new Allele("T", true);
    private final static Allele ALT_T = new Allele("T", false);
    private final static Allele REF_N = new Allele("N", true);
    private final static Allele ALT_N = new Allele("N", false);
    private final static Allele REF_NULL = new Allele(NULL_ALLELE_STRING, true);
    private final static Allele ALT_NULL = new Allele(NULL_ALLELE_STRING, false);
    public final static Allele NO_CALL = new Allele(NO_CALL_STRING, false);

    // ---------------------------------------------------------------------------------------------------------
    //
    // creation routines
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Create a new Allele that includes bases and if tagged as the reference allele if isRef == true.  If bases
     * == '-', a Null allele is created.  If bases ==  '.', a no call Allele is created.
     *
     * @param bases the DNA sequence of this variation, '-', of '.'
     * @param isRef should we make this a reference allele?
     * @throws IllegalArgumentException if bases contains illegal characters or is otherwise malformated
     */
    public static Allele create(byte[] bases, boolean isRef) {
        if ( bases == null )
            throw new IllegalArgumentException("create: the Allele base string cannot be null; use new Allele() or new Allele(\"\") to create a Null allele");

        if ( bases.length == 1 ) {
            // optimization to return a static constant Allele for each single base object
            switch (bases[0]) {
                case '.':
                    if ( isRef ) throw new IllegalArgumentException("Cannot tag a NoCall allele as the reference allele");
                    return NO_CALL;
                case '-': return isRef ? REF_NULL : ALT_NULL;
                case 'A': case 'a' : return isRef ? REF_A : ALT_A;
                case 'C': case 'c' : return isRef ? REF_C : ALT_C;
                case 'G': case 'g' : return isRef ? REF_G : ALT_G;
                case 'T': case 't' : return isRef ? REF_T : ALT_T;
                case 'N': case 'n' : return isRef ? REF_N : ALT_N;
                default: throw new IllegalArgumentException("Illegal base: " + (char)bases[0]);
            }
        } else {
            return new Allele(bases, isRef);
        }
    }

    public static Allele create(byte base, boolean isRef) {
//    public Allele(byte base, boolean isRef) {
        return create( new byte[]{ base }, isRef);
    }

    public static Allele create(byte base) {
        return create( base, false );
    }

    public static Allele extend(Allele left, byte[] right) {
        if (left.isSymbolic())
            throw new IllegalArgumentException("Cannot extend a symbolic allele");
        byte[] bases = null;
        if ( left.length() == 0 )
            bases = right;
        else {
            bases = new byte[left.length() + right.length];
            System.arraycopy(left.getBases(), 0, bases, 0, left.length());
            System.arraycopy(right, 0, bases, left.length(), right.length);
        }

        return create(bases, left.isReference());
    }

    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent the null allele
     */
    public static boolean wouldBeNullAllele(byte[] bases) {
        return (bases.length == 1 && bases[0] == '-') || bases.length == 0;
    }

    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent the NO_CALL allele
     */
    public static boolean wouldBeNoCallAllele(byte[] bases) {
        return bases.length == 1 && bases[0] == '.';
    }

    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent a symbolic allele
     */
    public static boolean wouldBeSymbolicAllele(byte[] bases) {
        return bases.length > 2 && bases[0] == '<' && bases[bases.length-1] == '>';
    }

    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent the well formatted allele
     */
    public static boolean acceptableAlleleBases(String bases) {
        return acceptableAlleleBases(bases.getBytes());
    }

    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent the well formatted allele
     */
    public static boolean acceptableAlleleBases(byte[] bases) {
        if ( wouldBeNullAllele(bases) || wouldBeNoCallAllele(bases) || wouldBeSymbolicAllele(bases) )
            return true;

        for ( int i = 0; i < bases.length; i++ ) {
            switch (bases[i]) {
                case 'A': case 'C': case 'G': case 'T': case 'N' : case 'a': case 'c': case 'g': case 't': case 'n' :
                    break;
                default:
                    return false;
            }
        }

        return true;
    }

    /**
     * @see Allele(byte[], boolean)
     *
     * @param bases  bases representing an allele
     * @param isRef  is this the reference allele?
     */
    public static Allele create(String bases, boolean isRef) {
    //public Allele(String bases, boolean isRef) {
        return create(bases.getBytes(), isRef);
    }


    /**
     * Creates a non-Ref allele.  @see Allele(byte[], boolean) for full information
     *
     * @param bases  bases representing an allele
     */
    public static Allele create(String bases) {
        return create(bases, false);
    }

    /**
     * Creates a non-Ref allele.  @see Allele(byte[], boolean) for full information
     *
     * @param bases  bases representing an allele
     */
    public static Allele create(byte[] bases) {
        return create(bases, false);
        //this(bases, false);
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // accessor routines
    //
    // ---------------------------------------------------------------------------------------------------------

    //Returns true if this is the null allele
    public boolean isNull()             { return isNull; }
    // Returns true if this is not the null allele
    public boolean isNonNull()          { return ! isNull(); }

    // Returns true if this is the NO_CALL allele
    public boolean isNoCall()           { return isNoCall; }
    // Returns true if this is not the NO_CALL allele
    public boolean isCalled()           { return ! isNoCall(); }

    // Returns true if this Allele is the reference allele
    public boolean isReference()        { return isRef; }
    // Returns true if this Allele is not the reference allele
    public boolean isNonReference()     { return ! isReference(); }

    // Returns true if this Allele is symbolic (i.e. no well-defined base sequence)
    public boolean isSymbolic()         { return isSymbolic; }

    // Returns a nice string representation of this object
    public String toString() {
        return (isNull() ? NULL_ALLELE_STRING : ( isNoCall() ? NO_CALL_STRING : getDisplayString() )) + (isReference() ? "*" : "");
    }

    /**
     * Return the DNA bases segregating in this allele.  Note this isn't reference polarized,
     * so the Null allele is represented by a vector of length 0
     *
     * @return the segregating bases
     */
    public byte[] getBases() { return isSymbolic ? EMPTY_ALLELE_BASES : bases; }

    /**
     * Return the DNA bases segregating in this allele in String format.
     * This is useful, because toString() adds a '*' to reference alleles and getBases() returns garbage when you call toString() on it.
     *
     * @return the segregating bases
     */
    public String getBaseString() { return new String(getBases()); }

    /**
     * Return the printed representation of this allele.
     * Same as getBaseString(), except for symbolic alleles.
     * For symbolic alleles, the base string is empty while the display string contains <TAG>.
     *
     * @return the allele string representation
     */
    public String getDisplayString() { return new String(bases); }

    /**
     * @param other  the other allele
     *
     * @return true if these alleles are equal
     */
    public boolean equals(Object other) {
        return ( ! (other instanceof Allele) ? false : equals((Allele)other, false) );
    }

    /**
     * @return hash code
     */
    public int hashCode() {
        int hash = 1;
        for (int i = 0; i < bases.length; i++)
            hash += (i+1) * bases[i];
        return hash;
    }

    /**
     * Returns true if this and other are equal.  If ignoreRefState is true, then doesn't require both alleles has the
     * same ref tag
     *
     * @param other            allele to compare to
     * @param ignoreRefState   if true, ignore ref state in comparison
     * @return true if this and other are equal
     */
    public boolean equals(Allele other, boolean ignoreRefState) {
        return this == other || (isRef == other.isRef || ignoreRefState) && isNull == other.isNull && isNoCall == other.isNoCall && (bases == other.bases || Arrays.equals(bases, other.bases));
    }

    /**
     * @param test  bases to test against
     *
     * @return  true if this Alelle contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
     */
    public boolean basesMatch(byte[] test) { return !isSymbolic && (bases == test || Arrays.equals(bases, test)); }

    /**
     * @param test  bases to test against
     *
     * @return  true if this Alelle contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
     */
    public boolean basesMatch(String test) { return basesMatch(test.toUpperCase().getBytes()); }

    /**
     * @param test  allele to test against
     *
     * @return  true if this Alelle contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
     */
    public boolean basesMatch(Allele test) { return basesMatch(test.getBases()); }

    /**
     * @return the length of this allele.  Null and NO_CALL alleles have 0 length.
     */
    public int length() {
        return isSymbolic ? 0 : bases.length;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // useful static functions
    //
    // ---------------------------------------------------------------------------------------------------------

    public static Allele getMatchingAllele(Collection<Allele> allAlleles, String alleleBases) {
        return getMatchingAllele(allAlleles, alleleBases.getBytes());
    }

    public static Allele getMatchingAllele(Collection<Allele> allAlleles, byte[] alleleBases) {
        for ( Allele a : allAlleles ) {
            if ( a.basesMatch(alleleBases) ) {
                return a;
            }
        }

        if ( wouldBeNoCallAllele(alleleBases) )
            return NO_CALL;
        else
            return null;    // couldn't find anything
    }

    public static List<Allele> resolveAlleles(List<Allele> possibleAlleles, List<String> alleleStrings) {
        List<Allele> myAlleles = new ArrayList<Allele>(alleleStrings.size());

        for ( String alleleString : alleleStrings ) {
            Allele allele = getMatchingAllele(possibleAlleles, alleleString);

            if ( allele == null ) {
                if ( Allele.wouldBeNoCallAllele(alleleString.getBytes()) ) {
                    allele = create(alleleString);
                } else {
                    throw new IllegalArgumentException("Allele " + alleleString + " not present in the list of alleles " + possibleAlleles);
                }
            }

            myAlleles.add(allele);
        }

        return myAlleles;
    }

    public int compareTo(Allele other) {
        if ( isReference() && other.isNonReference() )
            return -1;
        else if ( isNonReference() && other.isReference() ) 
            return 1;
        else
            return getBaseString().compareTo(other.getBaseString()); // todo -- potential performance issue
    }

    public static boolean oneIsPrefixOfOther(Allele a1, Allele a2) {
        if ( a1.isNull() || a2.isNull() )
            return true;

        if ( a2.length() >= a1.length() )
            return firstIsPrefixOfSecond(a1, a2);
        else
            return firstIsPrefixOfSecond(a2, a1);
    }

    private static boolean firstIsPrefixOfSecond(Allele a1, Allele a2) {
        String a1String = a1.getBaseString();
        return a2.getBaseString().substring(0, a1String.length()).equals(a1String);
    }
}
