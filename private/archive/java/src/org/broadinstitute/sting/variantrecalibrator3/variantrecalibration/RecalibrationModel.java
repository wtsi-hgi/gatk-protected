package org.broadinstitute.sting.gatk.walkers.variantrecalibrator3.variantrecalibration;

import org.apache.log4j.Logger;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 2/23/12
 * Time: 4:48 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class RecalibrationModel {

    protected final static Logger logger = Logger.getLogger(RecalibrationModel.class);
    public final boolean useFinal;

    public RecalibrationModel(boolean isFinal){
        useFinal = isFinal;
    }

    public abstract void initialize(List<VariantDatum> data, VariantRecalibratorArgumentCollection vrac);

    /**
     *
     * @param datum - attributes of a variant context
     * @return Log10 of probability that the datum represents a true variant site
     */
    public abstract double evaluateDatum(VariantDatum datum);
}
