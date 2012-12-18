package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.CallSet;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.NA12878DBArgumentCollection;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.VariantContext;

public class ImportReviews extends NA12878DBWalker {
    /**
     * Variants from this VCF file are used by this tool as input.
     * The file must at least contain the standard VCF header lines, but
     * can be empty (i.e., no variants are contained in the file).
     */
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=false)
    public RodBinding<VariantContext> variants;

    @Override
    public NA12878DBArgumentCollection.DBType getDefaultDB() {
        return NA12878DBArgumentCollection.DBType.DEV;
    }

    @Override
    public void initialize() {
        super.initialize();

        for ( final VCFHeaderLine line : GATKVCFUtils.getHeaderFields(getToolkit()) ) {
            if ( CallSet.isVCFHeaderLine(line) )
                ensureCallset(new CallSet(line));
        }
    }

    /**
     * Should be in header
     */
    private void ensureCallset(final CallSet callSet) {
        if ( db.getCallSet(callSet.getName()) == null ) {
            logger.info("Adding new reviewer callset " + callSet);
            db.addCallset(callSet);
        }
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) return 0;

        for ( VariantContext vc : tracker.getValues(variants, ref.getLocus()) ) {
            if ( vc.filtersWereApplied() )
                throw new UserException.MalformedVCF("Reviewed VCF shouldn't have filters applied");

            if ( ! vc.isBiallelic() ) {
                logger.info("Skipping unsupported non-ref / multi-allelic variant " + vc);
                continue;
            }

            final MongoVariantContext mvc = MongoVariantContext.createFromReview(vc);
            logger.info("adding reviewed site " + mvc);
            db.addCall(mvc);
        }

        return 1;
    }
}