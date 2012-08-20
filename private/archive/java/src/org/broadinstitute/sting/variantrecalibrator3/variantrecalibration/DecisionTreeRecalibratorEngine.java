package org.broadinstitute.sting.gatk.walkers.variantrecalibrator3.variantrecalibration;

import org.apache.log4j.Logger;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 2/27/12
 * Time: 1:27 PM
 * To change this template use File | Settings | File Templates.
 */
public class DecisionTreeRecalibratorEngine extends VariantRecalibratorEngine<DecisionTree> {

    private DecisionTree tree;
    private GMMRecalibratorEngine gmmEngine;
    private final static Logger logger = Logger.getLogger(DecisionTreeRecalibratorEngine.class);

    public DecisionTreeRecalibratorEngine(VariantRecalibratorArgumentCollection vrac) {
        super(vrac);
        gmmEngine = new GMMRecalibratorEngine(vrac);
    }

    public DecisionTreeRecalibratorEngine clone() {
        DecisionTreeRecalibratorEngine engine = new DecisionTreeRecalibratorEngine(super.VRAC);
        engine.gmmEngine = gmmEngine;
        return engine;
    }

    public void evaluateData(List<VariantDatum> data, DecisionTree tree, boolean contrastive) {
        // no such thing as contrastive
        if ( contrastive ) {
            throw new IllegalArgumentException("Decision Tree or Random Forest cannot be contrastive");
        }
        for ( VariantDatum d : data ) {
            double lod = tree.evaluateDatum(d);
            if ( tree.useFinal )
                d.setFinalLod(lod);
            else
                d.setInitialLod(lod);
        }
    }

    public void calculateWorstPerformingAnnotation(List<VariantDatum> data) {
        // unclear what to do here, so do nothing atm
    }

    public DecisionTree getModel(boolean positive) {
        if ( ! positive )
            throw new IllegalArgumentException("Decision tree or random forest do not have negative models");

        return tree;
    }

    public void trainClassifier(VariantDataManager manager, boolean isFinal) {
        logger.info(String.format("Generating decision tree on training set + worst %.2f percentf of data",VRAC.PERCENT_BAD_VARIANTS_FOR_RF));
        tree = generateModel(manager.getTrainingData(),isFinal);
        boolean splitHappened = true;
        while (splitHappened) {
            logger.info("Splitting an additional leaf...");
            splitHappened = tree.splitNext();
            logger.info("Leaves: "+Integer.toString(DecisionTree.leaves(tree.getTree()).size()));
        }
        tree.debug();
        logger.info("Evaluating data using tree. Leaves: "+ Integer.toString(DecisionTree.leaves(tree.getTree()).size()));
        this.evaluateData(manager.getData(),tree,false);
    }

    public DecisionTree generateModel(List<VariantDatum> data,boolean isFinal) {
        DecisionTree newTree = new DecisionTree(isFinal);
        newTree.initialize(data,VRAC);
        return newTree;
    }
}
