/*
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.variantrecalibrator3.variantrecalibration;

import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.walkers.bqsr.RecalDataManager;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class GMMRecalibratorEngine extends VariantRecalibratorEngine<GaussianMixtureModel> {

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////

    protected final static Logger logger = Logger.getLogger(GMMRecalibratorEngine.class);
    public final static double MIN_ACCEPTABLE_LOD_SCORE = -20000.0;

    private final static double MIN_PROB_CONVERGENCE = 2E-2;
    private GaussianMixtureModel goodModel = null;
    private GaussianMixtureModel badModel = null;

    /////////////////////////////
    // Public Methods to interface with the Engine
    /////////////////////////////

    public GMMRecalibratorEngine(final VariantRecalibratorArgumentCollection VRAC) {
        super(VRAC);
    }

    public GMMRecalibratorEngine clone() {
        GMMRecalibratorEngine toRet = new GMMRecalibratorEngine(super.VRAC);
        toRet.goodModel = goodModel;
        toRet.badModel = badModel;
        return toRet;
    }

    @Deprecated
    protected GaussianMixtureModel generateModel(final List<VariantDatum> data, boolean isFinal) {
        throw new ReviewedStingException("This implementaiton is deprecated. Please use the one allowing ngaussians immediately specified.");
    }

    protected GaussianMixtureModel generateModel(final List<VariantDatum> data, boolean isFinal, int nGaussians) {
        double[] annots = isFinal ? data.get(0).finalAnnotations : data.get(0).initialAnnotations;
        return new GaussianMixtureModel( nGaussians, annots.length, VRAC.SHRINKAGE, VRAC.DIRICHLET_PARAMETER, VRAC.PRIOR_COUNTS,isFinal );
    }

    public void evaluateData( final List<VariantDatum> data, final GaussianMixtureModel model, final boolean evaluateContrastively) {
        if( !model.isModelReadyForEvaluation ) {
            try {
                model.precomputeDenominatorForEvaluation();
            } catch( Exception e ) {
                logger.debug(e);
                logger.debug(e.getClass());
                logger.debug(e.getMessage());
                model.failedToConverge = true;
                return;
            }
        }
        
        logger.info("Evaluating full set of " + data.size() + " variants...");
        for( final VariantDatum datum : data ) {
            final double thisLod = evaluateDatum( datum, model );
            if( Double.isNaN(thisLod) ) {
                model.failedToConverge = true;
                //return;
            }

            double lod = ( evaluateContrastively ?
                            ( Double.isInfinite(datum.getLod()) ? // positive model said negative infinity
                                    ( MIN_ACCEPTABLE_LOD_SCORE + GenomeAnalysisEngine.getRandomGenerator().nextDouble() * MIN_ACCEPTABLE_LOD_SCORE ) // Negative infinity lod values are possible when covariates are extremely far away from their tight Gaussians
                                    : datum.prior + datum.getLod() - thisLod) // contrastive evaluation: (prior + positive model - negative model)
                            : thisLod ); // positive model only so set the lod and return
            if ( model.useFinal )
                datum.setFinalLod(lod);
            else
                datum.setInitialLod(lod);
        }
    }

    public GaussianMixtureModel getModel(boolean wantGoodModel) {
        return wantGoodModel ? goodModel : badModel;
    }

    public void calculateWorstPerformingAnnotation( final List<VariantDatum> data) {
        for( final VariantDatum datum : data ) {
            int worstAnnotation = -1;
            double[] annotations = goodModel.useFinal ? datum.finalAnnotations : datum.initialAnnotations;
            double minProb = Double.MAX_VALUE;
            for( int iii = 0; iii < annotations.length; iii++ ) {
                final Double goodProbLog10 = goodModel.evaluateDatumInOneDimension(datum, iii);
                logger.debug(datum.toString());
                logger.debug(String.format("%s",badModel));
                final Double badProbLog10 = badModel != null ? badModel.evaluateDatumInOneDimension(datum, iii) : 0.0;
                if( goodProbLog10 != null && badProbLog10 != null ) {
                    final double prob = goodProbLog10 - badProbLog10;
                    if(prob < minProb) { minProb = prob; worstAnnotation = iii; }
                }
            }
            datum.worstAnnotation = worstAnnotation;
        }
    }


    /////////////////////////////
    // Private Methods used for generating a GaussianMixtureModel
    /////////////////////////////

    public void trainClassifier(VariantDataManager dataManager,boolean isFinal) {
        ExpandingArrayList<VariantDatum> trainingData = dataManager.getTrainingData();
        goodModel = generateModel(trainingData,isFinal,VRAC.POS_MAX_GAUSSIANS);
        this.variationalBayesExpectationMaximization(goodModel,trainingData);
        this.evaluateData(dataManager.getData(),goodModel,false);

        if ( goodModel.failedToConverge ) {
            throw new UserException("NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. Please consider raising the number of variants used to train the negative model (via --percentBadVariants 0.05, for example) or lowering the maximum number of Gaussians to use in the model (via --maxGaussians 4, for example)");
        }

        if ( ! VRAC.noNegativeModel ) {

            final ExpandingArrayList<VariantDatum> negativeTrainingData = dataManager.selectWorstVariants(VRAC.PERCENT_BAD_VARIANTS,VRAC.MIN_NUM_BAD_VARIANTS);
            badModel = this.generateModel(negativeTrainingData,isFinal,VRAC.NEG_MAX_GAUSSIANS);
            this.variationalBayesExpectationMaximization(badModel,negativeTrainingData);
            this.evaluateData(dataManager.getData(),badModel,true);

            // Detect if the negative model failed to converge because of too few points and/or too many Gaussians and try again
            while( badModel.failedToConverge && VRAC.NEG_MAX_GAUSSIANS > 1 ) {
                logger.info("Negative model failed to converge. Retrying...");
                VRAC.NEG_MAX_GAUSSIANS--;
                badModel = this.generateModel( negativeTrainingData,isFinal, VRAC.NEG_MAX_GAUSSIANS );
                this.evaluateData( dataManager.getData(), goodModel, false );
                this.evaluateData( dataManager.getData(), badModel, true );
            }

            if( badModel.failedToConverge ) {
                logger.warn("Bad model consistently failing to converge; using flat prior. NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. Please consider raising the number of variants used to train the negative model (via --percentBadVariants 0.05, for example) or lowering the maximum number of Gaussians to use in the model (via --maxGaussians 4, for example)");
                this.evaluateData( dataManager.getData(), goodModel, false );
            }
        }
    }

    private void variationalBayesExpectationMaximization( final GaussianMixtureModel model, final List<VariantDatum> data) {
        model.initializeRandomModel(data, VRAC.NUM_KMEANS_ITERATIONS);
        model.normalizePMixtureLog10();
        // The VBEM loop
        model.expectationStep( data );
        double currentChangeInMixtureCoefficients;
        int iteration = 0;
        logger.info("Finished iteration " + iteration + ".");
        while( iteration < VRAC.MAX_ITERATIONS ) {
            iteration++;
            model.maximizationStep( data );
            currentChangeInMixtureCoefficients = model.normalizePMixtureLog10();
            model.expectationStep( data );
            if( iteration % 5 == 0 ) { // cut down on the number of output lines so that users can read the warning messages
                logger.info("Finished iteration " + iteration + ". \tCurrent change in mixture coefficients = " + String.format("%.5f", currentChangeInMixtureCoefficients));
            }
            if( iteration > 2 && currentChangeInMixtureCoefficients < MIN_PROB_CONVERGENCE ) {
                logger.info("Convergence after " + iteration + " iterations!");
                break;
            }
        }

        model.evaluateFinalModelParameters( data );
    }

    /////////////////////////////
    // Private Methods used for evaluating data given a GaussianMixtureModel
    /////////////////////////////

    private double evaluateDatum( final VariantDatum datum, final GaussianMixtureModel model ) {
        return model.evaluateDatum( datum );
    }
}
