package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.Arrays;
import java.util.Date;
import java.util.List;

public class ImportCallset extends NA12878DBWalker {
    @Argument(shortName = "isReviewer", required=false)
    protected boolean isReviewer = false;

    // for creating callsets
    @Argument(shortName = "callSetName", required=false)
    public String callSetName = "anonymous";

    @Argument(shortName = "assumedCallTruth", required = false)
    public TruthStatus assumedCallTruth = TruthStatus.UNKNOWN;

    @Argument(shortName = "howToTreatFilteredSites", required = true)
    public FilteredTreatment howToTreatFilteredSites;

    public enum FilteredTreatment {
        SKIP,
        FALSE_POSITIVE
    }

    @Argument(shortName = "howToTreatAC0", required = true)
    public AC0Treatment howToTreatAC0;

    public enum AC0Treatment {
        SKIP,
        MARK_AS_NON_POLYMORPHIC,
        FALSE_POSITIVE
    }

    /**
     * Variants from this VCF file are used by this tool as input.
     * The file must at least contain the standard VCF header lines, but
     * can be empty (i.e., no variants are contained in the file).
     */
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=false)
    public RodBinding<VariantContext> variants;

    private CallSet callSet = null;

    @Override
    public NA12878DBArgumentCollection.DBType getDefaultDB() {
        return NA12878DBArgumentCollection.DBType.DEV;
    }

    public void initialize() {
        super.initialize();

        callSet = db.getCallSet(callSetName);
        if ( callSet == null ) {
            callSet = new CallSet(callSetName, new Date(), isReviewer);
            logger.info("Adding new callset " + callSet);
            db.addCallset(callSet);
        }
        callSet = db.getCallSet(callSetName);

        if ( callSet == null )
            throw new ReviewedStingException("Callset is null for callSetName " + callSetName);

        logger.info("Using callset  " + callSet);
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) return 0;

        for ( VariantContext vc : tracker.getValues(variants, ref.getLocus()) ) {
            if ( ! vc.isVariant() ) {
                logger.info("Skipping unsupported non-variant site " + vc);
                continue;
            }

            if ( ! vc.isBiallelic() && vc.hasGenotypes() && vc.isIndel() ) {
                logger.info("Skipping unsupported multi-allelic indel site with genotypes " + vc);
                continue;
            }

            if ( vc.isSymbolic() ) {
                logger.info("Skipping unsupported symbolic variant " + vc);
                continue;
            }

            final List<VariantContext> biallelics = VariantContextUtils.splitVariantContextToBiallelics(vc);
            final boolean isMultiAllelic = biallelics.size() > 1;
            if ( isMultiAllelic )
                logger.info("Splitting original multi-allelic " + vc);

            for ( final VariantContext biallelic : biallelics ) {
                importCall(isMultiAllelic ? addNA12878Genotype(vc, biallelic) : biallelic);
            }
        }

        return 1;
    }

    /**
     * Add a semi-reasonable genotype for NA12878 to a derived biallelic VC based on its original
     * genotype in a multi-allelic VC.
     *
     * In some cases (hom-ref or hom-alt) this is easy.  If NA12878 is het for ref and one alt its also easy.  But
     * het-alt cannot really be handled currently, and gets no called.
     *
     * TODO Does not work when derived bi has had its alleles trimmed for indels
     * TODO -- could be handled by allele matching inside this code...
     *
     * @param multi a multi-allelic variant context
     * @param bi a derived bi-allelic variant context
     * @return a bi-allelic variant with a semi-reasonable NA12878 genotype
     */
    private VariantContext addNA12878Genotype(final VariantContext multi, final VariantContext bi) {
        if ( multi.getNAlleles() < 3 )
            throw new IllegalArgumentException("addNA12878Genotype expected the original variant context to be multi-allelic " + multi);

        if ( bi.getNAlleles() != 2 )
            throw new IllegalArgumentException("addNA12878Genotype expected the biallelic variant context to actually be biallelic " + bi);

        final Genotype na12878 = multi.getGenotype("NA12878");

        if ( na12878 == null )
            return bi; // no NA12878 genotype

        final Genotype genotype;
        final Allele gt1 = na12878.getAllele(0);
        final Allele gt2 = na12878.getAllele(1);

        if ( bi.getAlleles().contains(gt1) || bi.getAlleles().contains(gt2) ) {
            switch ( na12878.getType() ) {
                // could be het matching these alleles or het-alt (1/2 e.g.)
                case HET:
                    if ( bi.getAlleles().contains(gt1) && bi.getAlleles().contains(gt2) )
                        genotype = na12878;
                    else {
                        logger.warn("Bailing out on het-alt genotype " + na12878 + " at " + bi + " with original VCF record " + multi);
                        genotype = MongoGenotype.NO_CALL;
                    }
                    break;
                /// both cases are the same for the biallelic split
                case HOM_REF:
                case HOM_VAR:
                    genotype = na12878;
                    break;
                case MIXED:
                case NO_CALL:
                case UNAVAILABLE:
                    genotype = MongoGenotype.NO_CALL;
                default:
                    throw new IllegalStateException("Unexpected genotype type " + na12878);
            }
        } else {
            // na12878 doesn't contain this allele at all, best we can represent is UNKNOWN
            genotype = MongoGenotype.NO_CALL;
        }

        return new VariantContextBuilder(bi).genotypes(genotype).make();
    }

    /**
     * Actually import a bi-allelic variant context into the knowledge base
     * @param vc a bi-allelic variant context that may or may not contain a genotype for NA12878
     */
    private void importCall(final VariantContext vc) {
        if ( vc.getNAlleles() != 2 )
            throw new IllegalArgumentException("Can only import biallelic variant contexts but got " + vc);

        Genotype gt = vc.hasGenotype("NA12878") ? vc.getGenotype("NA12878") : MongoGenotype.NO_CALL;

        TruthStatus type;
        if ( vc.isFiltered() ) {
            if ( howToTreatFilteredSites == FilteredTreatment.SKIP ) {
                logger.info("Skipping filtered site " + vc);
                return;
            } else {
                // assumed FP
                type = TruthStatus.FALSE_POSITIVE;
                gt = monomorphic(vc);
            }
        } else if ( isAC0(vc) ) {
            gt = monomorphic(vc);
            switch (howToTreatAC0) {
                case SKIP:
                    return;
                case MARK_AS_NON_POLYMORPHIC:
                    type = assumedCallTruth;
                    break;
                case FALSE_POSITIVE:
                    type = TruthStatus.FALSE_POSITIVE;
                    break;
                default:
                    throw new IllegalStateException("Unexpected AC0 treatment " + howToTreatAC0);
            }
        } else {
            // PASS sites use the assumed call status type
            type = assumedCallTruth;
        }

        final MongoVariantContext mvc = MongoVariantContext.create(callSet.getName(), vc, type, gt);
        logger.info("adding call " + mvc);
        db.addCall(mvc);
    }

    /**
     * Create a monomorphic genotype for NA12878 based on the vc
     * @param vc a variant context
     * @return a monomorphic Genotype for NA12878
     */
    private Genotype monomorphic(final VariantContext vc) {
        return new GenotypeBuilder("NA12878").alleles(Arrays.asList(vc.getReference(), vc.getReference())).make();
    }

    /**
     * Is the VC monomorphic w.r.t. NA12878 at this location?
     *
     * @param vc a bi-allelic variant context that has a genotype for NA12878 or an AC field, or none
     * @return true if we believe the VC has allele count 0 for for the alt allele
     */
    private boolean isAC0(final VariantContext vc) {
        if ( ! vc.isBiallelic() )
            throw new IllegalArgumentException("VariantContext should be biallelic " + vc);

        return (vc.hasGenotype("NA12878") && vc.getGenotype("NA12878").isHomRef())
                || vc.getAttributeAsInt("AC", -1) == 0;
    }
}