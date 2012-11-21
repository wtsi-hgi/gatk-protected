package org.broadinstitute.sting.gatk.walkers.na12878kb;

import com.mongodb.ReflectionDBObject;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.util.*;

/**
 * VariantContext consistent with Mongo
 *
 * Is a very limited subset of the functionality of a full variant context
 *
 * User: depristo
 * Date: 11/5/12
 * Time: 6:27 AM
 */
public class MongoVariantContext extends ReflectionDBObject {
    private List<String> supportingCallsets = new ArrayList<String>(1);
    private String chr;
    private int start, stop;
    private String ref, alt;
    private Date date;
    private TruthStatus mongoType;
    private boolean reviewed;
    private MongoGenotype gt;

    public MongoVariantContext() { }

    protected static MongoVariantContext create(final String callSetName,
                                                final String chr,
                                                final int start,
                                                final String ref,
                                                final String alt,
                                                final boolean isReviewed) {
        return create(callSetName, chr, start, ref, alt, MongoGenotype.NO_CALL, isReviewed);
    }

    protected static MongoVariantContext create(final String callSetName,
                                                final String chr,
                                                final int start,
                                                final String ref,
                                                final String alt,
                                                final Genotype gt,
                                                final boolean isReviewed) {
        final int stop = start + ref.length() - 1;
        final VariantContextBuilder vcb = new VariantContextBuilder(callSetName, chr, start, stop, Arrays.asList(Allele.create(ref, true), Allele.create(alt)));
        return new MongoVariantContext(callSetName, vcb.make(), TruthStatus.TRUE_POSITIVE, new Date(), gt, isReviewed);
    }

    public static MongoVariantContext create(final List<String> supportingCallsets,
                                             final VariantContext vc,
                                             final TruthStatus type,
                                             final Date date,
                                             final Genotype gt,
                                             final boolean isReviewed) {
        return new MongoVariantContext(supportingCallsets, vc, type, date, gt, isReviewed);
    }

    public static MongoVariantContext create(final String callSetName,
                                             final VariantContext vc,
                                             final TruthStatus type,
                                             final Genotype gt) {
        return create(Arrays.asList(callSetName), vc, type, new Date(), gt, false);
    }

    public static MongoVariantContext createFromReview(final VariantContext vc) {
        final String callSet = parseReviewField(vc, "CallSetName");
        final TruthStatus type = TruthStatus.valueOf(parseReviewField(vc, "TruthStatus"));
        final Date date = new Date(Long.valueOf(parseReviewField(vc, "Date")));
        final Genotype gt = vc.hasGenotype("NA12878") ? vc.getGenotype("NA12878") : MongoGenotype.NO_CALL;
        return new MongoVariantContext(callSet, vc, type, date, gt, true);
    }

    protected MongoVariantContext(final String callset,
                                  final VariantContext vc,
                                  final TruthStatus type,
                                  final Date date,
                                  final Genotype gt,
                                  final boolean isReviewed) {
        this(Arrays.asList(callset), vc, type, date, gt, isReviewed);
    }

    protected MongoVariantContext(final List<String> supportingCallsets,
                                  final VariantContext vc,
                                  final TruthStatus type,
                                  final Date date,
                                  final Genotype gt,
                                  final boolean isReviewed) {
        if ( vc.getNAlleles() > 2 )
            throw new ReviewedStingException("MongoVariantContext only supports single alt allele, but saw " + vc);

        if ( vc.isSymbolic() )
            throw new ReviewedStingException("MongoVariantContext doesn't support symbolic alleles but got " + vc);

        this.supportingCallsets = supportingCallsets;
        this.chr = vc.getChr();
        this.start = vc.getStart();
        this.stop = vc.getEnd();
        this.ref = vc.getReference().getDisplayString();
        this.alt = vc.getAlternateAllele(0).getDisplayString();
        this.mongoType = type;
        this.date = date;
        this.reviewed = isReviewed;
        this.gt = new MongoGenotype(vc.getAlleles(), gt);
    }

    private List<Allele> getAlleles() {
        return Arrays.asList(getRefAllele(), getAltAllele());
    }

    public boolean isSingleCallset() {
        return supportingCallsets.size() == 1;
    }

    public String getCallSetName() {
        return Utils.join(",", getSupportingCallSets());
    }

    public List<String> getSupportingCallSets() {
        return supportingCallsets;
    }

    public void setSupportingCallSets(List<String> supportingCallSets) {
        this.supportingCallsets = supportingCallSets;
    }

    @Override
    public String toString() {
        return String.format("%s:%d-%d %s/%s %s/%s reviewed?=%b %s %s", getChr(), getStart(), getStop(), getRef(), getAlt(),
                getType(), getPolymorphicStatus(), isReviewed(), Utils.join(",", getSupportingCallSets()), getGt());
    }

    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getStop() {
        return stop;
    }

    public void setStop(int stop) {
        this.stop = stop;
    }

    public String getRef() {
        return ref;
    }

    public Allele getRefAllele() {
        return Allele.create(getRef(), true);
    }

    public void setRef(String ref) {
        this.ref = ref;
    }

    public String getAlt() {
        return alt;
    }

    public Allele getAltAllele() {
        return Allele.create(getAlt(), false);
    }

    public void setAlt(String alt) {
        this.alt = alt;
    }

    public VariantContext getVariantContext() {
        final VariantContextBuilder vcb = new VariantContextBuilder(getCallSetName(), chr, start, stop, getAlleles());

        if ( isReviewed() )
            addReviewInfoFields(vcb);

        vcb.genotypes(gt.toGenotype(getAlleles()));

        return vcb.make();
    }

    public TruthStatus getType() {
        return mongoType;
    }

    protected void addReviewInfoFields(final VariantContextBuilder vcb) {
        vcb.attribute("CallSetName", getCallSetName());
        vcb.attribute("TruthStatus", getType());
        vcb.attribute("PolymorphicStatus", getPolymorphicStatus());
        vcb.attribute("Date", getDate().getTime());
        //vcb.attribute("PhredConfidence", getPhredConfidence());
    }

    protected static Set<VCFHeaderLine> reviewHeaderLines() {
        final Set<VCFHeaderLine> lines = new HashSet<VCFHeaderLine>();

        lines.add(new VCFInfoHeaderLine("CallSetName", 1, VCFHeaderLineType.String, "Name of the review"));
        lines.add(new VCFInfoHeaderLine("TruthStatus", 1, VCFHeaderLineType.String, "What is the truth state of this call"));
        lines.add(new VCFInfoHeaderLine("PolymorphicStatus", 1, VCFHeaderLineType.String, "Is this call polymorphic in NA12878"));
        lines.add(new VCFInfoHeaderLine("Date", 1, VCFHeaderLineType.String, "Date/time as a long of this review"));
        lines.add(new VCFInfoHeaderLine("PhredConfidence", 1, VCFHeaderLineType.Integer, "Phred-scaled confidence in this review"));

        return lines;
    }

    private static String parseReviewField(final VariantContext vc, final String key) {
        final String value = vc.getAttributeAsString(key, null);
        if ( value == null )
            throw new UserException.BadInput("Missing key " + key + " from " + vc);
        return value;
    }

    public String getMongoType() {
        return mongoType.toString();
    }

    public void setMongoType(String mongoType) {
        this.mongoType = TruthStatus.valueOf(mongoType);
    }

    public PolymorphicStatus getPolymorphicStatus() {
        return getGt().getPolymorphicStatus();
    }

    public boolean isReviewed() {
        return reviewed;
    }

    public boolean getReviewed() {
        return reviewed;
    }

    public void setReviewed(boolean reviewed) {
        this.reviewed = reviewed;
    }

    public GenomeLoc getLocation(final GenomeLocParser parser) {
        return parser.createGenomeLoc(getChr(), getStart(), getStop());
    }

    public Date getDate() {
        return date;
    }

    public void setDate(Date date) {
        this.date = date;
    }

    public MongoGenotype getGt() {
        return gt;
    }

    public void setGt(MongoGenotype gt) {
        this.gt = gt;
    }

    public boolean isPolymorphic() { return getGt().isPolymorphic(); }
    public boolean isMonomorphic() { return getGt().isMonomorphic(); }
    public boolean isDiscordant() { return getGt().isDiscordant(); }
    public boolean isUnknown() { return getGt().isUnknown(); }

    public boolean matches(final VariantContext vc) {
        return vc.isBiallelic() &&
                getChr().equals(vc.getChr()) &&
                getStart() == vc.getStart() &&
                getStop() == vc.getEnd() &&
                getRef().equals(vc.getReference().getBaseString()) &&
                getAlt().equals(vc.getAlternateAllele(0).getBaseString());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        MongoVariantContext that = (MongoVariantContext) o;

        if (reviewed != that.reviewed) return false;
        if (start != that.start) return false;
        if (stop != that.stop) return false;
        if (alt != null ? !alt.equals(that.alt) : that.alt != null) return false;
        if (chr != null ? !chr.equals(that.chr) : that.chr != null) return false;
        if (date != null ? !date.equals(that.date) : that.date != null) return false;
        if (gt != null ? !gt.equals(that.gt) : that.gt != null) return false;
        if (mongoType != that.mongoType) return false;
        if (ref != null ? !ref.equals(that.ref) : that.ref != null) return false;
        if (supportingCallsets != null ? !supportingCallsets.equals(that.supportingCallsets) : that.supportingCallsets != null)
            return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = supportingCallsets != null ? supportingCallsets.hashCode() : 0;
        result = 31 * result + (chr != null ? chr.hashCode() : 0);
        result = 31 * result + start;
        result = 31 * result + stop;
        result = 31 * result + (ref != null ? ref.hashCode() : 0);
        result = 31 * result + (alt != null ? alt.hashCode() : 0);
        result = 31 * result + (date != null ? date.hashCode() : 0);
        result = 31 * result + (mongoType != null ? mongoType.hashCode() : 0);
        result = 31 * result + (reviewed ? 1 : 0);
        result = 31 * result + (gt != null ? gt.hashCode() : 0);
        return result;
    }
}
