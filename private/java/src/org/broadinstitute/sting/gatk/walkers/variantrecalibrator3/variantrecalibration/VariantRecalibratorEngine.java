package org.broadinstitute.sting.gatk.walkers.variantrecalibrator3.variantrecalibration;

import org.broadinstitute.sting.gatk.walkers.bqsr.RecalDataManager;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 2/23/12
 * Time: 4:34 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class VariantRecalibratorEngine<T extends RecalibrationModel> {
    // the unified argument collection
    final protected VariantRecalibratorArgumentCollection VRAC;

    public VariantRecalibratorEngine(VariantRecalibratorArgumentCollection args) {
        this.VRAC = args;
    }

    public abstract T getModel(boolean wantGoodModel);

    protected abstract T generateModel(List<VariantDatum> data, boolean isFinal);

    public abstract void evaluateData(List<VariantDatum> data, T positive, boolean contrastive);

    public abstract void calculateWorstPerformingAnnotation(List<VariantDatum> data);

    public abstract void trainClassifier(VariantDataManager data, boolean isFinal);
}
