package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.Date;

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

            if ( ! vc.isBiallelic() ) {
                logger.info("Skipping unsupported multi-allelic site " + vc);
                continue;
            }

            if ( vc.isSymbolic() ) {
                logger.info("Skipping unsupported symbolic variant " + vc);
                continue;
            }

            Genotype gt = vc.hasGenotype("NA12878") ? vc.getGenotype("NA12878") : MongoGenotype.NO_CALL;

            TruthStatus type;
            if ( vc.isFiltered() ) {
                if ( howToTreatFilteredSites == FilteredTreatment.SKIP ) {
                    logger.info("Skipping filtered site " + vc);
                    continue;
                } else {
                    // assumed FP
                    type = TruthStatus.FALSE_POSITIVE;
                    gt = monomorphic(vc);
                }
            } else if ( isAC0(vc) ) {
                gt = monomorphic(vc);
                switch (howToTreatAC0) {
                    case SKIP:
                        continue;
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

        return 1;
    }

    private Genotype monomorphic(final VariantContext vc) {
        return new GenotypeBuilder("NA12878").alleles(Arrays.asList(vc.getReference(), vc.getReference())).make();
    }

    private boolean isAC0(final VariantContext vc) {
        return (vc.hasGenotype("NA12878") && vc.getGenotype("NA12878").isHomRef())
                || vc.getAttributeAsInt("AC", -1) == 0;
    }
}