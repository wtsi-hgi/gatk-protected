package org.broadinstitute.sting.gatk.walkers.na12878kb;

import com.mongodb.ReflectionDBObject;
import org.broadinstitute.sting.gatk.walkers.na12878kb.errors.MongoVariantContextException;
import org.broadinstitute.sting.utils.BaseUtils;
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
public class MongoVariantContext extends ReflectionDBObject implements Cloneable {
    /**
     * A list of at least 1 String describing the call set (i.e., OMNIPoly or CEU_best_practices).
     *
     * May contains multiple strings if the site is present in multiple callsets in the consensus
     */
    private List<String> supportingCallsets = new ArrayList<String>(1);

    /**
     * The chromosomal position of the site, following b37 conventions (i.e., chromosomes are named 1, 2, etc).
     */
    private String chr;

    /**
     * The start and stop position of the site (exactly the VCF convention)
     */
    private int start, stop;

    /**
     * Ref and alt alleles following the VCF convention.
     *
     * Must be in all caps
     *
     * For SNPs, ref is the base of the genome reference, and alt is the other allele, such as ref='A' and alt = 'C'
     * For indels, ref contains at a minimum the base in the genome reference.  For deletions ref contains additionall
     * the deleted reference bases, while alt is just the reference base again.  It's swapped for insertions.  For
     * example, a 2 bp deletion looks like ref=ACT alt=A, while for the symmetric 2 bp insertion you see
     * ref=A alt=ACT.
     *
     * Note that there can be only one alt allele per site.  VCF records with multiple alt alleles
     * must be split into separate records.
     */
    private String ref, alt;

    /**
     * The truth status of this site, can adopt the following values
     *
     * TRUE_POSITIVE -- a confirmed true positive call
     * FALSE_POSITIVE -- a confirmed false positive call
     * UNKNOWN -- nothing is known about the truth of this site
     * SUSPECT -- a likely false positive, but we aren't sure (for reviews)
     * DISCORDANT -- for the consensus only
     */
    private TruthStatus mongoType;

    /**
     * The genotype of NA12878 at this site (see MongoGenotype for more information)
     */
    private MongoGenotype gt;

    /**
     * The date this site was created.  Can be missing, in which case the date is assumed to be now
     */
    private Date date = new Date();

    /**
     * (Optional) reviewed status.  If true, indicates that this site represent information from a
     * manual review of data, and will be interpreted as more informative than a non-reviewed site.
     * If missing, assumed to *not* to be a reviewed site
     */
    private boolean reviewed = false;

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

    public MongoVariantContext(List<String> supportingCallsets, String chr, int start, int stop, String ref, String alt, TruthStatus mongoType, MongoGenotype gt, Date date, boolean reviewed) {
        this.supportingCallsets = supportingCallsets;
        this.chr = chr;
        this.start = start;
        this.stop = stop;
        this.ref = ref;
        this.alt = alt;
        this.mongoType = mongoType;
        this.gt = gt;
        this.date = date;
        this.reviewed = reviewed;
    }

    @Override
    protected MongoVariantContext clone() throws CloneNotSupportedException {
        return new MongoVariantContext(supportingCallsets, chr, start, stop, ref, alt, mongoType, gt, date, reviewed);
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
        try {
            return String.format("%s:%d-%d %s/%s %s/%s reviewed?=%b %s %s", getChr(), getStart(), getStop(), getRef(), getAlt(),
                    getType(), getPolymorphicStatus(), isReviewed(), Utils.join(",", getSupportingCallSets()), getGt());
        } catch ( Exception e ) {
            return String.format("%s at %s:%d [malformed]", supportingCallsets, chr, start);
        }
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

    /**
     * Convert this MongoVariantContext into an honest-to-goodness VariantContext
     *
     * @return a validate VariantContext
     */
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

    /**
     * Is this MongoVariantContext a match to VariantContext vc
     *
     * A match means that this and vc have the same location and same ref / alt alleles
     *
     * @param vc a non-null VariantContext vc
     * @return true if this and vc match, false otherwise
     */
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

    /**
     * Make sure this MongoVariantContext is valid, throwing a MongoVariantContextException if not
     *
     * @throws org.broadinstitute.sting.gatk.walkers.na12878kb.errors.MongoVariantContextException if this is malformed
     * @param parser a GenomeLocParser so we know what contigs are allowed
     */
    protected void validate(final GenomeLocParser parser) {
        if ( supportingCallsets == null || supportingCallsets.size() == 0 )
            error("SupportingCallSets has a bad value %s", supportingCallsets);
        if ( supportingCallsets.indexOf(null) != -1 )
            error("SupportingCallSets contains a null element %s", supportingCallsets);
        if ( start < 1 ) error("Start = %d < 1", start);
        if ( start > stop ) error("Start %d > Stop %d", start, stop);
        if ( chr == null ) error("Chr is null");
        if ( ! parser.contigIsInDictionary(chr) ) error("Chr %s is not in the b37 dictionary", chr);
        if ( ref == null || ref.equals("") ) error("ref allele is null or empty string");
        if ( ! BaseUtils.isUpperCase(ref.getBytes()) ) error("ref allele must be all upper case but got " + ref);
        if ( ! Allele.acceptableAlleleBases(ref) ) error("ref allele contains unacceptable bases " + ref);
        if ( alt == null || alt.equals("") ) error("alt allele is null or empty string");
        if ( ! BaseUtils.isUpperCase(alt.getBytes()) ) error("alt allele must be all upper case but got " + alt);
        if ( ! Allele.acceptableAlleleBases(alt) ) error("alt allele contains unacceptable bases " + alt);
        if ( gt == null ) error("gt is null");
        final String gtBad = gt.validate();
        if ( gtBad != null ) error("gt %s is bad: %s", gt, gtBad);
    }

    private void error(final String message, final Object ... values) {
        throw new MongoVariantContextException(
                String.format("MongoVariantContext is bad: reason = %s at context %s",
                        String.format(message, values), this.toString()),
                this);
    }
}
