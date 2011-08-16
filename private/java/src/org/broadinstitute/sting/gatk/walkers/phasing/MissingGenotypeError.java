package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.Collection;

public class MissingGenotypeError extends RodWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @Input(fullName="truth", shortName = "truth", doc="truth", required=false)
    public RodBinding<VariantContext> truth;

    @Input(fullName="missing", shortName = "missing", doc="missing", required=false)
    public RodBinding<VariantContext> missing;

    @Input(fullName="complete", shortName = "complete", doc="complete", required=false)
    public RodBinding<VariantContext> complete;

    private int numMarkers = 0;
    private int numErrors = 0;

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            Collection<VariantContext> truthVCs = tracker.getValues(truth, context.getLocation());
            Collection<VariantContext> missingVCs = tracker.getValues(missing, context.getLocation());
            Collection<VariantContext> completeVCs = tracker.getValues(complete, context.getLocation());

            VariantContext truthVC = truthVCs.iterator().next();
            VariantContext missingVC = missingVCs.iterator().next();
            VariantContext completeVC = completeVCs.iterator().next();

            Genotype truthG = truthVC.getGenotype("NA12878");
            Genotype missingG = missingVC.getGenotype("NA12878");
            Genotype completeG = completeVC.getGenotype("NA12878");

            if (truthG.isCalled() && missingG.isNoCall() && completeG.isCalled()) {
                numMarkers++;

                if (!completeG.sameGenotype(truthG)) {
                    numErrors++;
                }
            }
        }

        return null;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public Integer reduceInit() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void onTraversalDone(Integer sum) {
        out.printf("numMarkers=%d\tnumErrors=%d\tpctErrors=%f%n", numMarkers, numErrors, 100.0 * ((double) numErrors / (double) numMarkers));
    }
}
