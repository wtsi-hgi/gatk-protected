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

package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import com.mongodb.ReflectionDBObject;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.List;

/**
 * Genotype consistent with MongoDB
 *
 * Has a very limited subset of the functionality of a full Genotype object.
 *
 * Contains the following four fields:
 *
 * The genotype of NA12878 at this site, encoded as two ints, allele1 and allele2
 * An optional GQ value (int) and an optional DP (int) value, both of which
 * default to the value -1, meaning not present.
 *
 * Important information about the encoding:
 *
 * If allele1 and allele2 are -1, this means that the genotype is UNKNOWN (equivalent to the VCF encoding ./.)
 * If allele1 and allele2 are 0 or 1, this implies that NA12878 has a known genotype corresponding to the
 * VCF genotype encoding allele1/allele2.  So, if the site has ref/alt alleles and NA12878 is alt/alt here
 * then allele1 and allele2 should both be equal to 1.
 * If GQ == 0, then the genotype is considered discordant, but this should only be used by the consensus algorithm
 *
 * User: depristo
 * Date: 11/5/12
 * Time: 6:27 AM
 */
public class MongoGenotype extends ReflectionDBObject {
    private static final String SAMPLE_NAME = "NA12878";
    public final static Genotype NO_CALL = GenotypeBuilder.createMissing(SAMPLE_NAME, 2);
    public final static int DISCORDANT_GQ = 0;

    int allele1 = -1, allele2 = -1;
    int GQ = -1;
    int DP = -1;

    public static Genotype create(final VariantContext vc, int allele1, int allele2) {
        return new MongoGenotype(allele1, allele2).toGenotype(vc.getAlleles());
    }

    /**
     * Create a discordant Genotype based on input genotype g
     * @param g a genotype to make a discordant version of
     * @return a discordant version of g
     */
    public static Genotype createDiscordant(final Genotype g) {
        return new GenotypeBuilder(g.getSampleName(), g.getAlleles()).GQ(DISCORDANT_GQ).make();
    }

    /**
     * For MongoDB set() building approach.  Not for public consumption
     */
    public MongoGenotype() { }

    /**
     * Create a MongoGenotype from a VariantContext's list of alleles and a corresponding Genotype gt
     *
     * @param alleles list of alleles from VariantContext.getAlleles()
     * @param gt the Genotype we will use as the basis for this MongoGenotype
     */
    public MongoGenotype(final List<Allele> alleles, final Genotype gt) {
        if ( gt.getPloidy() == 0 ) {
            this.allele1 = this.allele2 = -1;
        } else if ( gt.getPloidy() != 2 ) {
            throw new IllegalArgumentException("Ploidy must be two for conversion to MongoGenotype " + gt);
        } else {
            this.allele1 = alleles.indexOf(gt.getAllele(0));
            this.allele2 = alleles.indexOf(gt.getAllele(1));
        }
        this.DP = gt.hasDP() ? gt.getDP() : -1;
        this.GQ = gt.hasGQ() ? gt.getGQ() : -1;
        validate();
    }

    /**
     * Create a simple MongoGenotype with alleles allele1 and allele2
     * @param allele1 the allele index of allele1
     * @param allele2 the allele index of allele2
     */
    public MongoGenotype(int allele1, int allele2) {
        this(allele1, allele2, -1, -1);
        validate();
    }

    /**
     * Full constructor: create a MongoGenotype with alleles allele1 and allele2, genotype quality and depth
     * @param allele1 the allele index of allele1
     * @param allele2 the allele index of allele2
     * @param GQ genotype quality, must be >= -1 (-1 means missing)
     * @param DP depth of the sequencing data supporting this case, must be >= -1 (-1 means missing)
     */
    public MongoGenotype(int allele1, int allele2, int GQ, int DP) {
        this.allele1 = allele1;
        this.allele2 = allele2;
        this.GQ = GQ;
        this.DP = DP;
    }

    /**
     * Convert this MongoGenotype to a VariantContext Genotype object
     * @param alleles the list of alleles from the VariantContext.getAlleles()
     * @return a Genotype corresponding to this MongoGenotype
     */
    public Genotype toGenotype(final List<Allele> alleles) {
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME);
        gb.alleles(Arrays.asList(getAllele(alleles, allele1), getAllele(alleles, allele2)));
        if ( DP != -1 ) gb.DP(DP);
        if ( GQ != -1 ) gb.GQ(GQ);
        return gb.make();
    }

    // -------------------------------------------------------------------------------------
    //
    // Polymorphic status information from this MongoDB
    //
    // -------------------------------------------------------------------------------------

    public PolymorphicStatus getPolymorphicStatus() {
        if ( isUnknown() ) return PolymorphicStatus.UNKNOWN;
        if ( isPolymorphic() ) return PolymorphicStatus.POLYMORPHIC;
        if ( isMonomorphic() ) return PolymorphicStatus.MONOMORPHIC;
        if ( isDiscordant() ) return PolymorphicStatus.DISCORDANT;
        throw new IllegalStateException("Expected polymorphic state " + this);
    }

    public boolean isUnknown() {
        return allele1 == -1;
    }

    public boolean isPolymorphic() {
        return (allele1 > 0 || allele2 > 0) && ! isDiscordant();
    }

    public boolean isMonomorphic() {
        return allele1 == 0 && allele2 == 0 && ! isDiscordant();
    }

    public boolean isDiscordant() {
        return GQ == DISCORDANT_GQ;
    }

    // -------------------------------------------------------------------------------------
    //
    // MongoDB getter / setters
    //
    // -------------------------------------------------------------------------------------

    public int getAllele1() {
        return allele1;
    }

    public void setAllele1(int allele1) {
        this.allele1 = allele1;
    }

    public int getAllele2() {
        return allele2;
    }

    public void setAllele2(int allele2) {
        this.allele2 = allele2;
    }

    public int getGQ() {
        return GQ;
    }

    public void setGQ(int GQ) {
        this.GQ = GQ;
    }

    public int getDP() {
        return DP;
    }

    public void setDP(int DP) {
        this.DP = DP;
    }

    private Allele getAllele(final List<Allele> alleles, final int i) {
        return i == -1 ? Allele.NO_CALL : alleles.get(i);
    }

    @Override
    public String toString() {
        return alleleString(allele1) + "/" + alleleString(allele2) +
                (isDiscordant() ? " DISCORDANT" : option(", GQ=", GQ)) +
                option(", DP=", DP);
    }

    private String option(final String prefix, int v) { return v == -1 ? "" : prefix + String.valueOf(v); }
    private String alleleString(int i) { return i == -1 ? "." : String.valueOf(i); }

    protected String validate() {
        if ( allele1 < -1 || allele1 > 1 ) return "allele1 " + allele1 + " not between -1 and 1";
        else if ( allele2 < -1 || allele2 > 1 ) return "allele2 " + allele2 + " not between -1 and 1";
        else if ( allele1 == -1 && allele2 != -1 ) return "Both allele1 and allele2 must be -1 if one is but saw " + allele1 + "/" + allele2;
        else if ( allele2 == -1 && allele1 != -1 ) return "Both allele1 and allele2 must be -1 if one is but saw " + allele1 + "/" + allele2;
        else if ( GQ < -1 ) return "GQ " + GQ + " < -1";
        else if ( DP < -1 ) return "DP " + DP + " < -1";
        return null;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        MongoGenotype that = (MongoGenotype) o;

        if (DP != that.DP) return false;
        if (GQ != that.GQ) return false;
        if (allele1 != that.allele1) return false;
        if (allele2 != that.allele2) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = allele1;
        result = 31 * result + allele2;
        result = 31 * result + GQ;
        result = 31 * result + DP;
        return result;
    }
}
