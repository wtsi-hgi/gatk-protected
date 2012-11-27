package org.broadinstitute.sting.gatk.walkers.na12878kb;

import com.mongodb.ReflectionDBObject;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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
        return allele1 == 0 && allele2 == 0;
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
                option(", GQ=", GQ) +
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
