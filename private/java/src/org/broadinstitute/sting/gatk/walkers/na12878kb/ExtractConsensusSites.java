package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.NA12878DBArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.PolymorphicStatus;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.TruthStatus;
import org.broadinstitute.sting.utils.variantcontext.GenotypeType;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;

import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * Extract a VCF from the NA12878 Knowledge Base meeting criteria
 *
 * See @link http://gatkforums.broadinstitute.org/discussion/1848/using-the-na12878-knowledge-base for more information
 */
public class ExtractConsensusSites extends NA12878DBWalker {
    @Output
    public VariantContextWriter out;

    @Argument(fullName="maxSites", shortName = "maxSites", doc="Max. number of bad sites to write out", required=false)
    public int maxSites = 10000;

    @Argument(fullName="variantType", shortName = "variantType", doc="", required=false)
    public VariantContext.Type variantType = null;

    @Argument(fullName="includeCallset", shortName = "includeCallset", doc="", required=false)
    public Set<String> includeCallset = null;

    @Argument(fullName="excludeCallset", shortName = "excludeCallset", doc="", required=false)
    public Set<String> excludeCallset = null;

    @Argument(fullName="truthStatus", shortName = "truthStatus", doc="", required=false)
    public TruthStatus truthStatus = null;

    @Argument(fullName="polymorphicStatus", shortName = "polymorphicStatus", doc="", required=false)
    public PolymorphicStatus polymorphicStatus = null;

    @Argument(fullName="genotype", shortName = "genotype", doc="", required=false)
    public GenotypeType genotypeType = null;

    @Argument(fullName="uniqueToOneCallset", shortName = "uniqueToOneCallset", doc="", required=false)
    public boolean uniqueToOneCallset = false;

    private abstract class ShouldBeReviewed {
        public abstract boolean exclude(final MongoVariantContext mvc, final VariantContext vc);
    }

    private List<ShouldBeReviewed> criteria = new LinkedList<ShouldBeReviewed>();

    @Override
    public void initialize() {
        super.initialize();

        if ( variantType != null ) criteria.add(new ByType());
        if ( excludeCallset != null ) criteria.add(new ByExclude());
        if ( includeCallset != null ) criteria.add(new ByInclude());
        if ( truthStatus != null ) criteria.add(new ByTruthStatus());
        if ( uniqueToOneCallset ) criteria.add(new ByUniqueToOneCallset());
        if ( genotypeType != null ) criteria.add(new ByGenotype());
        if ( polymorphicStatus != null ) criteria.add(new ByPolymorphicStatus());

        out.writeHeader(db.makeStandardVCFHeader());
    }

    @Override public boolean isDone() { return true; }

    @Override
    public void onTraversalDone(Integer result) {
        int nWritten = 0;
        for ( final MongoVariantContext mvc : db.getConsensusSites(makeSiteSelector())) {
            final VariantContext vc = mvc.getVariantContext();
            if ( shouldBeReviewed(mvc, vc) ) {
                out.add(vc);
                if ( nWritten++ > maxSites)
                    break;
            }
        }

        super.onTraversalDone(result);
    }

    private boolean shouldBeReviewed(final MongoVariantContext mvc, final VariantContext vc) {
        for ( final ShouldBeReviewed shouldBeReviewed : criteria )
            if ( shouldBeReviewed.exclude(mvc, vc) )
                return false;

        return true;
    }

    @Override
    public NA12878DBArgumentCollection.DBType getDefaultDB() {
        return NA12878DBArgumentCollection.DBType.PRODUCTION;
    }

    private class ByType extends ShouldBeReviewed {
        @Override
        public boolean exclude(MongoVariantContext mvc, VariantContext vc) {
            return vc.getType() != variantType;
        }
    }

    private class ByInclude extends ShouldBeReviewed {
        @Override
        public boolean exclude(MongoVariantContext mvc, VariantContext vc) {
            for ( final String callset : mvc.getSupportingCallSets() ) {
                if ( includeCallset.contains(callset)) {
                    return false;
                }
            }
            return true;
        }
    }
    private class ByExclude extends ShouldBeReviewed {
        @Override
        public boolean exclude(MongoVariantContext mvc, VariantContext vc) {
            for ( final String callset : mvc.getSupportingCallSets() ) {
                if ( excludeCallset.contains(callset)) {
                    return true;
                }
            }
            return false;
        }
    }

    private class ByUniqueToOneCallset extends ShouldBeReviewed {
        @Override
        public boolean exclude(MongoVariantContext mvc, VariantContext vc) {
            return mvc.getSupportingCallSets().size() != 1;
        }
    }

    private class ByGenotype extends ShouldBeReviewed {
        @Override
        public boolean exclude(MongoVariantContext mvc, VariantContext vc) {
            return vc.getGenotype("NA12878").getType() != genotypeType;
        }
    }

    private class ByTruthStatus extends ShouldBeReviewed {
        @Override
        public boolean exclude(MongoVariantContext mvc, VariantContext vc) {
            return mvc.getType() != truthStatus;
        }
    }

    private class ByPolymorphicStatus extends ShouldBeReviewed {
        @Override
        public boolean exclude(MongoVariantContext mvc, VariantContext vc) {
            return mvc.getPolymorphicStatus() != polymorphicStatus;
        }
    }

}