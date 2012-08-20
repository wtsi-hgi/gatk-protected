package org.broadinstitute.sting.gatk.walkers.variantrecalibrator3.variantrecalibration;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/2/12
 * Time: 12:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class RandomForestRecalibrationEngine extends VariantRecalibratorEngine<RandomForestBridge> {

    RandomForestBridge model;

    public RandomForestRecalibrationEngine(VariantRecalibratorArgumentCollection args) {
        super(args);
    }

    public RandomForestRecalibrationEngine clone() {
        RandomForestRecalibrationEngine engine = new RandomForestRecalibrationEngine(super.VRAC);
        engine.model = model;
        return engine;
    }

    public void evaluateData(List<VariantDatum> data, RandomForestBridge forest, boolean isContrastive) {
        for ( VariantDatum d : data ) {
            if ( forest.useFinal ) {
                d.setFinalLod(forest.evaluateDatum(d));
            } else {
                d.setInitialLod(forest.evaluateDatum(d));
            }
        }
    }

    public void trainClassifier(VariantDataManager manager, boolean finalAttribs ) {
        model = new RandomForestBridge(finalAttribs,manager);
        model.initialize(manager.getData(),VRAC);
        this.evaluateData(manager.getData(),model,false);
    }

    public void calculateWorstPerformingAnnotation(List<VariantDatum> data) {
        // do nothing yet
    }

    public RandomForestBridge getModel(boolean ignore) {
        return model;
    }

    public RandomForestBridge generateModel(List<VariantDatum> data, boolean isContrastive) {
        throw new IllegalStateException("generateModel doesn't work for RFRE -- bad implementation");
    }
}
