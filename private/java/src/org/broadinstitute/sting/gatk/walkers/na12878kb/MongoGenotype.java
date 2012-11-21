package org.broadinstitute.sting.gatk.walkers.na12878kb;

import com.mongodb.ReflectionDBObject;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.List;

/**
 * Genotype consistent with Mongo
 *
 * Is a very limited subset of the functionality of a full genotype
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

    public MongoGenotype() { }

    public MongoGenotype(final List<Allele> alleles, final Genotype gt) {
        this.allele1 = alleles.indexOf(gt.getAllele(0));
        this.allele2 = alleles.indexOf(gt.getAllele(1));
        this.DP = gt.hasDP() ? gt.getDP() : -1;
        this.GQ = gt.hasGQ() ? gt.getGQ() : -1;
    }

    public MongoGenotype(int allele1, int allele2) {
        this(allele1, allele2, -1, -1);
    }

    public MongoGenotype(int allele1, int allele2, int GQ, int DP) {
        this.allele1 = allele1;
        this.allele2 = allele2;
        this.GQ = GQ;
        this.DP = DP;
    }

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

    public Genotype toGenotype(final List<Allele> alleles) {
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME);
        gb.alleles(Arrays.asList(getAllele(alleles, allele1), getAllele(alleles, allele2)));
        if ( DP != -1 ) gb.DP(DP);
        if ( GQ != -1 ) gb.GQ(GQ);
        return gb.make();
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
