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
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.MongoVariantContextException;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;

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
        return create(callSetName, chr, start, ref, alt, TruthStatus.TRUE_POSITIVE, gt, isReviewed);
    }

    protected static MongoVariantContext create(final String callSetName,
                                                final String chr,
                                                final int start,
                                                final String ref,
                                                final String alt,
                                                final TruthStatus truthStatus,
                                                final Genotype gt,
                                                final boolean isReviewed) {
        final int stop = start + ref.length() - 1;
        final VariantContextBuilder vcb = new VariantContextBuilder(callSetName, chr, start, stop, Arrays.asList(Allele.create(ref, true), Allele.create(alt)));
        return new MongoVariantContext(callSetName, vcb.make(), truthStatus, new Date(), gt, isReviewed);
    }

    public static MongoVariantContext create(final List<String> supportingCallsets,
                                             final VariantContext vc,
                                             final TruthStatus truthStatus,
                                             final Date date,
                                             final Genotype gt,
                                             final boolean isReviewed) {
        return new MongoVariantContext(supportingCallsets, vc, truthStatus, date, gt, isReviewed);
    }

    public static MongoVariantContext create(final String callSetName,
                                             final VariantContext vc,
                                             final TruthStatus truthStatus,
                                             final Genotype gt) {
        return create(Arrays.asList(callSetName), vc, truthStatus, new Date(), gt, false);
    }

    public static MongoVariantContext createFromReview(final VariantContext vc) {
        final String callSet = parseReviewField(vc, "CallSetName");
        final TruthStatus truthStatus = TruthStatus.valueOf(parseReviewField(vc, "TruthStatus"));
        final Date date = new Date(Long.valueOf(parseReviewField(vc, "Date")));
        final Genotype gt = vc.hasGenotype("NA12878") ? vc.getGenotype("NA12878") : MongoGenotype.NO_CALL;
        return new MongoVariantContext(callSet, vc, truthStatus, date, gt, true);
    }

    protected MongoVariantContext(final String callset,
                                  final VariantContext vc,
                                  final TruthStatus truthStatus,
                                  final Date date,
                                  final Genotype gt,
                                  final boolean isReviewed) {
        this(Arrays.asList(callset), vc, truthStatus, date, gt, isReviewed);
    }

    protected MongoVariantContext(final List<String> supportingCallsets,
                                  final VariantContext vc,
                                  final TruthStatus truthStatus,
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
        this.mongoType = truthStatus;
        this.date = date;
        this.reviewed = isReviewed;
        this.gt = new MongoGenotype(vc.getAlleles(), gt);

        //TODO Once the public gatk includes BaseUtils.isUpperCase, we can enable this
        //validate(null);
    }

    protected MongoVariantContext(List<String> supportingCallsets, String chr, int start, int stop, String ref, String alt, TruthStatus truthStatus, MongoGenotype gt, Date date, boolean reviewed) {
        this.supportingCallsets = supportingCallsets;
        this.chr = chr;
        this.start = start;
        this.stop = stop;
        this.ref = ref;
        this.alt = alt;
        this.mongoType = truthStatus;
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
        addReviewInfoFields(vcb);
        vcb.genotypes(gt.toGenotype(getAlleles()));
        return vcb.make();
    }

    public TruthStatus getType() {
        return mongoType;
    }

    public void setTruth(TruthStatus status) {
        this.mongoType = status;
    }

    protected void addReviewInfoFields(final VariantContextBuilder vcb) {
        vcb.attribute("CallSetName", getCallSetName());
        vcb.attribute("TruthStatus", getType());
        vcb.attribute("PolymorphicStatus", getPolymorphicStatus());
        vcb.attribute("Date", getDate().getTime());
        vcb.attribute("Reviewed", isReviewed());
        //vcb.attribute("PhredConfidence", getPhredConfidence());
    }

    public static Set<VCFHeaderLine> reviewHeaderLines() {
        final Set<VCFHeaderLine> lines = new HashSet<VCFHeaderLine>();

        lines.add(new VCFInfoHeaderLine("CallSetName", 1, VCFHeaderLineType.String, "Name of the review"));
        lines.add(new VCFInfoHeaderLine("TruthStatus", 1, VCFHeaderLineType.String, "What is the truth state of this call"));
        lines.add(new VCFInfoHeaderLine("PolymorphicStatus", 1, VCFHeaderLineType.String, "Is this call polymorphic in NA12878"));
        lines.add(new VCFInfoHeaderLine("Date", 1, VCFHeaderLineType.String, "Date/time as a long of this review"));
        lines.add(new VCFInfoHeaderLine("PhredConfidence", 1, VCFHeaderLineType.Integer, "Phred-scaled confidence in this review"));
        lines.add(new VCFInfoHeaderLine("Reviewed", 0, VCFHeaderLineType.Flag, "Was this a manually reviewed record?"));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY));

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

    /**
     * Is this and that duplicate entries?
     *
     * Duplicate entry contain all of the same essential information but may different in the time
     * they were added.  Sometimes multiple entries are added (for example, clicked save review twice in
     * IGV.  Such identical records have all the same field values but the date may be different.
     *
     * @param that another fully filled in MongoVariantContext to test for duplicate
     * @return true if this and that are duplicates, false otherwise
     */
    public boolean isDuplicate(final MongoVariantContext that) {
        if ( that == null ) throw new IllegalArgumentException("that cannot be null");

        if (reviewed != that.reviewed) return false;
        if (start != that.start) return false;
        if (stop != that.stop) return false;
        if (!ref.equals(that.ref)) return false;
        if (!alt.equals(that.alt)) return false;
        if (!chr.equals(that.chr)) return false;
        // ignores date if (date != null ? !date.equals(that.date) : that.date != null) return false;
        if (!gt.equals(that.gt)) return false;
        if (mongoType != that.mongoType) return false;
        if (!supportingCallsets.equals(that.supportingCallsets)) return false;

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
     * @throws org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.MongoVariantContextException if this is malformed
     * @param parser a GenomeLocParser so we know what contigs are allowed, can be null
     */
    protected void validate(final GenomeLocParser parser) {
        if ( supportingCallsets == null || supportingCallsets.size() == 0 )
            error("SupportingCallSets has a bad value %s", supportingCallsets);
        if ( supportingCallsets.indexOf(null) != -1 )
            error("SupportingCallSets contains a null element %s", supportingCallsets);
        if ( start < 1 ) error("Start = %d < 1", start);
        if ( start > stop ) error("Start %d > Stop %d", start, stop);
        if ( chr == null ) error("Chr is null");
        if ( parser != null && ! parser.contigIsInDictionary(chr) ) error("Chr %s is not in the b37 dictionary", chr);
        if ( parser == null && chr.toLowerCase().startsWith("chr") ) error("MongoVariantContext %s uses the UCSC convention -- must use the b37 convention (i.e., 1 not chr1)", chr);
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
