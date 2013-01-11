/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.variantrecalibrator3.variantrecalibration;

import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
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
