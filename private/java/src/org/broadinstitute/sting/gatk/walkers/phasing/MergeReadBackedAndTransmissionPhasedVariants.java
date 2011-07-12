package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class MergeReadBackedAndTransmissionPhasedVariants extends RodWalker<Integer, Integer> {
    private List<VariantContext> pbtCache = new ArrayList<VariantContext>();
    private List<VariantContext> rbpCache = new ArrayList<VariantContext>();

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            Collection<VariantContext> pbts = tracker.getVariantContexts(ref, "pbt", null, ref.getLocus(), true, true);
            Collection<VariantContext> rbps = tracker.getVariantContexts(ref, "rbp", null, ref.getLocus(), true, true);

            VariantContext pbt = pbts.iterator().hasNext() ? pbts.iterator().next() : null;
            VariantContext rbp = rbps.iterator().hasNext() ? rbps.iterator().next() : null;

            if (pbt != null && rbp != null) {
                /*
                pbtCache.add(pbt);
                rbpCache.add(rbp);

                if (pbtCache.size() > 1 && !pbt.isFiltered()) {

                }
                */

                logger.info(pbt.getGenotype(0) + " " + rbp.getGenotype(0));
            }
        }

        return null;
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }
}
