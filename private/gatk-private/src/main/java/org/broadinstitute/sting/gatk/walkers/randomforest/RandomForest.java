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

package org.broadinstitute.sting.gatk.walkers.randomforest;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.PriorityQueue;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 11/14/13
 */

/*
Each tree is constructed using the following algorithm:

        Let the number of training cases be N, and the number of variables in the classifier be M.
        We are told the number m of input variables to be used to determine the decision at a node of the tree; m should be much less than M.
        Choose a training set for this tree by choosing n times with replacement from all N available training cases (i.e., take a bootstrap sample). Use the rest of the cases to estimate the error of the tree, by predicting their classes.
        For each node of the tree, randomly choose m variables on which to base the decision at that node. Calculate the best split based on these m variables in the training set.
        Each tree is fully grown and not pruned (as may be done in constructing a normal tree classifier).
*/

public class RandomForest {
    final PriorityQueue<RandomForestDecisionNode> forestClassifier;
    final float PERCENT_TRAINING = 0.01f;
    final float PERCENT_TESTING = 0.01f;
    final int NUM_TRAINING_VARIANTS;
    final int NUM_TESTING_VARIANTS;
    final int MAX_TREE_DEPTH;
    private double sumAccuracy;
    private double CONVERGENCE_THRESHOLD = 1E-4;
    private int MAX_NUM_ITERATIONS = 50000;

    protected final static Logger logger = Logger.getLogger(RandomForest.class);

    /**
     * Construct a random forest classifier
     * @param trainingData  the training data to use to make this random forest
     * @param masterKeySet  the list of all possible annotation values to use
     * @param numTrees      the number of iterations to run
     */
    public RandomForest( final List<RandomForestDatum> trainingData, final LinkedHashSet<String> masterKeySet, final int numTrees ) {
        forestClassifier = new PriorityQueue<>(numTrees);
        NUM_TRAINING_VARIANTS = (int) Math.ceil( PERCENT_TRAINING * (float) trainingData.size() );
        NUM_TESTING_VARIANTS = (int) Math.ceil( PERCENT_TESTING * (float) trainingData.size() );
        MAX_TREE_DEPTH = (int) Math.floor( Math.log((double) masterKeySet.size()) );
        logger.info(String.format("Creating %d random trees each with depth %d by drawing %d examples out of %d training variants", numTrees, MAX_TREE_DEPTH, NUM_TRAINING_VARIANTS, trainingData.size()));

        final List<RandomForestDatum> positiveTrainingDataSNP = subsetToSNPs(subsetToSpecificTrainingStatus(trainingData, true), true);
        final List<RandomForestDatum> negativeTrainingDataSNP = subsetToSNPs(subsetToSpecificTrainingStatus(trainingData, false), true);
        final List<RandomForestDatum> positiveTrainingDataIndel = subsetToSNPs(subsetToSpecificTrainingStatus(trainingData, true), false);
        final List<RandomForestDatum> negativeTrainingDataIndel = subsetToSNPs(subsetToSpecificTrainingStatus(trainingData, false), false);
        int iteration = 0;
        boolean converged = false;
        double previousMeanAccuracy = 0.0;
        while( !converged && iteration++ < MAX_NUM_ITERATIONS ) {
            // Draw a training set for this tree by choosing n times with replacement from all available training cases
            final List<RandomForestDatum> trainingSample = MathUtils.randomSample(positiveTrainingDataSNP, NUM_TRAINING_VARIANTS / 4);
            trainingSample.addAll(MathUtils.randomSample(negativeTrainingDataSNP, NUM_TRAINING_VARIANTS / 4));
            trainingSample.addAll(MathUtils.randomSample(positiveTrainingDataIndel, NUM_TRAINING_VARIANTS / 4));
            trainingSample.addAll(MathUtils.randomSample(negativeTrainingDataIndel, NUM_TRAINING_VARIANTS / 4));
            // Draw a testing set for this tree by choosing n times without replacement from all available training cases
            final List<RandomForestDatum> testingSample = MathUtils.randomSubset(positiveTrainingDataSNP, NUM_TESTING_VARIANTS / 4);
            testingSample.addAll(MathUtils.randomSubset(negativeTrainingDataSNP, NUM_TESTING_VARIANTS / 4));
            testingSample.addAll(MathUtils.randomSubset(positiveTrainingDataIndel, NUM_TESTING_VARIANTS / 4));
            testingSample.addAll(MathUtils.randomSubset(negativeTrainingDataIndel, NUM_TESTING_VARIANTS / 4));

            if( trainingSample.isEmpty() || testingSample.isEmpty() ) {
                throw new IllegalStateException("There seems to be too little training data. Found " + trainingData.size() + " training variants.");
            }

            // Build a random decision tree using this smaller training set
            final RandomForestDecisionNode tree = new RandomForestDecisionNode(trainingSample, 0, MAX_TREE_DEPTH);
            // Evaluate the error of this tree using the full training set by predicting their classes
            tree.accuracy = tree.estimateAccuracy(testingSample);

            if( forestClassifier.size() < numTrees ) {
                forestClassifier.add(tree);
            } else if ( tree.compareTo(forestClassifier.peek()) > 0 ) {
                forestClassifier.poll();
                forestClassifier.add(tree);
            }

            if( iteration > 1 && (iteration) % 400 == 0 ) {
                sumAccuracy = sumAccuracy(forestClassifier);
                final double meanAccuracy = sumAccuracy / numTrees;
                logger.info("Completed iteration " + iteration + ". Random forest classifier contains " + forestClassifier.size() + " trees with mean accuracy = " + meanAccuracy);
                if( forestClassifier.size() >= numTrees && meanAccuracy - previousMeanAccuracy < CONVERGENCE_THRESHOLD ) {
                    logger.info("Convergence!");
                    converged = true;
                }
                previousMeanAccuracy = meanAccuracy;
            }
        }
        sumAccuracy = sumAccuracy(forestClassifier);
    }

    /**
     * Used to construct a forest directly for unit testing purposes
     * @param forestClassifier  a pre-made random forest classifier
     */
    protected RandomForest( final PriorityQueue<RandomForestDecisionNode> forestClassifier ) {
        NUM_TRAINING_VARIANTS = 0;
        NUM_TESTING_VARIANTS = 0;
        MAX_TREE_DEPTH = 1;
        this.forestClassifier = new PriorityQueue<>(forestClassifier.size());
        this.forestClassifier.addAll(forestClassifier);
        for( final RandomForestDecisionNode tree : this.forestClassifier ) {
            sumAccuracy += Math.abs(tree.accuracy);
        }
    }

    private double sumAccuracy( final PriorityQueue<RandomForestDecisionNode> forestClassifier ) {
        double accuracy = 0.0;
        for( final RandomForestDecisionNode tree : forestClassifier ) {
            accuracy += Math.abs(tree.accuracy);
        }
        return accuracy;
    }

    /**
     * Apply the random forest classifier to test data point and return its score
     * @param rfd   The test data point
     * @return      a double value between -1.0 and 1.0
     */
    public double classifyDatum( final RandomForestDatum rfd ) {
        double score = 0.0;
        for( final RandomForestDecisionNode tree : forestClassifier ) {
            score += tree.accuracy * ( tree.classifyDatum(rfd) ? 1.0 : -1.0 );
        }
        return score / sumAccuracy;
    }

    /**
     * Subset the list of data points to only those with the specified training status
     * @param data      list of data
     * @param isGood    the desired training status
     * @return          non-null output list
     */
    protected static List<RandomForestDatum> subsetToSpecificTrainingStatus(final List<RandomForestDatum> data, final boolean isGood) {
        final List<RandomForestDatum> returnData = new ArrayList<>();
        for( final RandomForestDatum rfd : data ) {
            if( (isGood && rfd.isGood && !rfd.isBad ) || (!isGood && rfd.isBad && !rfd.isGood) ) {
                returnData.add(rfd);
            }
        }
        return returnData;
    }

    /**
     * Subset the list of data points to only those with the specified variant type status (SNP versus non-SNP)
     * @param data      list of data
     * @param isSNP     the desired variant type status
     * @return          non-null output list
     */
    protected static List<RandomForestDatum> subsetToSNPs( final List<RandomForestDatum> data, final boolean isSNP ) {
        final List<RandomForestDatum> returnData = new ArrayList<>();
        for( final RandomForestDatum rfd : data ) {
            if( rfd.type.equals(VariantContext.Type.SNP) == isSNP ) {
                returnData.add(rfd);
            }
        }
        return returnData;
    }

    /**
     * Subset the list of data points to only those which are part of the training data
     * @param data      list of data
     * @return          non-null output list
     */
    protected static List<RandomForestDatum> subsetToTrainingData(final List<RandomForestDatum> data) {
        final List<RandomForestDatum> testData = new ArrayList<>();
        for( final RandomForestDatum rfd : data ) {
            if( (rfd.isGood || rfd.isBad) && !(rfd.isGood && rfd.isBad) ) {
                testData.add(rfd);
            }
        }
        return testData;
    }

    /**
     * Subset the list of data points to only those which are part of the input data
     * @param data      list of data
     * @return          non-null output list
     */
    protected static List<RandomForestDatum> subsetToInputData(final List<RandomForestDatum> data) {
        final List<RandomForestDatum> testData = new ArrayList<>();
        for( final RandomForestDatum rfd : data ) {
            if( rfd.isInputSite ) {
                testData.add(rfd);
            }
        }
        return testData;
    }

    /**
     * Subset the list of data points to only those which are part of the truth data
     * @param data      list of data
     * @return          non-null output list
     */
    protected static List<RandomForestDatum> subsetToTruthData(final List<RandomForestDatum> data) {
        final List<RandomForestDatum> testData = new ArrayList<>();
        for( final RandomForestDatum rfd : data ) {
            if( rfd.isTruthSite ) {
                testData.add(rfd);
            }
        }
        return testData;
    }
}
