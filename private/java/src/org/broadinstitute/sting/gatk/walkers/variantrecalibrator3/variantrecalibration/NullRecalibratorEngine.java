package org.broadinstitute.sting.gatk.walkers.variantrecalibrator3.variantrecalibration;

import org.apache.log4j.Logger;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/1/12
 * Time: 5:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class NullRecalibratorEngine extends VariantRecalibratorEngine<NullRecalibratorEngine.NullRecalibratorModel> {

    private NullRecalibratorModel model;
    final Logger logger = Logger.getLogger(NullRecalibratorEngine.class);

    public NullRecalibratorEngine(VariantRecalibratorArgumentCollection args) {
        super(args);
    }

    public NullRecalibratorEngine clone() {
        NullRecalibratorEngine engine = new NullRecalibratorEngine(super.VRAC);
        engine.model = model;

        return engine;
    }

    public void evaluateData(List<VariantDatum> data, NullRecalibratorModel model, boolean isContrastive) {
        for ( VariantDatum d : data ) {
            if ( model.useFinal ) {
                d.setFinalLod(model.evaluateDatum(d));
            } else {
                d.setInitialLod(0.0);
            }
        }
    }

    public void trainClassifier(VariantDataManager manager, boolean isFinal) {
        logger.debug("No secondary training engine provided, passing initial classification through");
        model = generateModel(manager.getTrainingData(),isFinal);
        evaluateData(manager.getData(),model,false);
    }

    public void calculateWorstPerformingAnnotation(List<VariantDatum> data) {

    }

    public NullRecalibratorModel getModel(boolean ignore) {
        return model;
    }

    public NullRecalibratorModel generateModel(List<VariantDatum> data, boolean isFinal) {
        return new NullRecalibratorModel(isFinal);
    }

    public class NullRecalibratorModel extends RecalibrationModel {

        public NullRecalibratorModel(boolean useFinal) {
            super(useFinal);
        }

        public void initialize(List<VariantDatum> data, VariantRecalibratorArgumentCollection vrac) {

        }

        public double evaluateDatum(VariantDatum datum) {
            return datum.getLod();
        }
    }
}
