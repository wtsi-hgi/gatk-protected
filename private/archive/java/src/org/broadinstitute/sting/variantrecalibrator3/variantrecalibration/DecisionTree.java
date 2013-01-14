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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 2/27/12
 * Time: 1:27 PM
 * To change this template use File | Settings | File Templates.
 */
public class DecisionTree extends RecalibrationModel {

    private DecisionTreeNode tree;
    private Logger logger = Logger.getLogger(DecisionTree.class);

    public DecisionTree(boolean isFinal) {
        super(isFinal);
    }

    public void initialize(List<VariantDatum> data, VariantRecalibratorArgumentCollection args) {
        int nTrue = 0;
        for ( VariantDatum datum : data ) {
            if ( datum.atPositiveTrainingSite )
                ++nTrue;
        }
        VariantDatum[] dataAr = new VariantDatum[data.size()];
        data.toArray(dataAr);
        tree = new DecisionTreeNode(nTrue,dataAr,useFinal);
    }

    public DecisionTreeNode getTree() { return tree; }

    public boolean splitNext() {
        // todo -- grossly inefficient; should probably store the information ratios in a nested hash map to avoid recalculating
        DecisionTreeNode bestLeaf = null;
        DataSplit bestSplit = null;
        double bestGainRatio = 0;
        double bestGain = 0;
        int leafIdx = -1;
        int bestLeafIdx = 0;
        for ( DecisionTreeNode leaf : leaves(tree) ) {
            ++leafIdx;
            if ( ! leaf.canSplit() )
                continue;
            VariantDatum[] data = leaf.getData();
            int lim = useFinal ? data[0].finalAnnotations.length : data[0].initialAnnotations.length;
            for ( int idx = 0; idx < lim; idx++ ) {
                logger.debug("Checking leaf "+Integer.toString(leafIdx)+" annot "+Integer.toString(idx));
                // now generate a preliminary split
                DataSplit subSplit = splitOnAttribute(data,idx,leaf.getNTrue());
                if ( subSplit != null && subSplit.gainRatio > bestGainRatio ) {
                    bestLeaf = leaf;
                    bestLeafIdx = leafIdx;
                    bestSplit = subSplit;
                    bestGain = subSplit.gain;
                    bestGainRatio = subSplit.gainRatio;
                }
            }
        }

        if ( bestSplit == null || bestSplit.chiSquare < 12.0 ) {
            return false;
        }
        logger.debug("Best split on leaf "+Integer.toString(bestLeafIdx) + " with index " + Integer.toString(bestSplit.attribIdx) + "with gain "+Double.toString(bestGain));

        bestLeaf.split(bestSplit.attribIdx, bestSplit.threshold, bestSplit.leftTrue, bestSplit.totalTrue - bestSplit.leftTrue,
                bestSplit.left, bestSplit.right);
        return true;
    }

    public void debug() {
        for ( DecisionTreeNode leaf : leaves(tree) ) {
            logger.debug(leaf);
        }
    }

    private DataSplit splitOnAttribute(VariantDatum[] data, final int index, int nTrue) {
        // want to find an optimal position on which to split the data. Can do this exhaustively
        // by checking (data-2) positions, by sampling randomly, or by checking k quantiles
        Arrays.sort(data,new Comparator<VariantDatum>() {
            @Override
            public int compare(VariantDatum variantDatum, VariantDatum variantDatum1) {
                if ( useFinal )
                    return Double.compare(variantDatum.finalAnnotations[index],variantDatum1.finalAnnotations[index]);
                return Double.compare(variantDatum.initialAnnotations[index],variantDatum1.initialAnnotations[index]);
            }
        });
        // calculate gain ratios in one pass
        double trueEntropy = entropy(((double) nTrue)/data.length);
        int leftIdx = 0;
        int nTrueLeft = 0;
        int rightIdx = getRightIdx(data,index,leftIdx);
        int bestLeft = 0;
        int bestLeftNT = 0;
        double bestRatio = 0;
        double bestGain = 0;
        double fracIncr = 1.0/data.length;
        double fracLeft = 0.0;
        while ( rightIdx < data.length-100 ) {
            // update n true left
            for ( int j = rightIdx-1; j >= leftIdx; j-- ) {
                nTrueLeft += data[j].atPositiveTrainingSite ? 1 : 0;
            }
            fracLeft += ((double)rightIdx-leftIdx)*fracIncr;
            double pTL = ((double)nTrueLeft)/(1+leftIdx);
            double pTR = ((double) nTrue - nTrueLeft)/(data.length-rightIdx);
            double gain = trueEntropy - fracLeft*entropy(pTL) - (1.0-fracLeft)*entropy(pTR);
            double ratio = gain/entropy(fracLeft);
            //logger.debug(String.format("Ratio: %.2f, Gain: %.2f",ratio,gain));
            if ( ratio > bestRatio && ! Double.isNaN(ratio) && leftIdx > 100 ) {
                bestLeft = leftIdx;
                bestRatio = ratio;
                bestLeftNT = nTrueLeft;
                bestGain = gain;
            }
            leftIdx = rightIdx;
            rightIdx = getRightIdx(data,index,leftIdx);
        }
        int bestRight = getRightIdx(data,index,bestLeft);
        double mdpt;
        if ( useFinal )
            mdpt = (data[bestRight].finalAnnotations[index]-data[bestLeft].finalAnnotations[index])/2;
        else
            mdpt = (data[bestRight].initialAnnotations[index]-data[bestLeft].initialAnnotations[index])/2;
        fracLeft = ((double) bestRight)/data.length;
        double leftPosExp = fracLeft*nTrue;
        double leftNegExp = fracLeft*(data.length-nTrue);
        double rightPosExp = (1.0-fracLeft)*nTrue;
        double rightNegExp = (1.0-fracLeft)*(data.length-nTrue);
        double chiSquare = Math.pow(leftPosExp-bestLeftNT,2)/leftPosExp +
                        Math.pow(leftNegExp-bestRight+bestLeftNT,2)/leftNegExp +
                        Math.pow(rightPosExp-(nTrue-bestLeftNT),2)/rightPosExp +
                        Math.pow(rightNegExp-((data.length-bestRight)-(nTrue-bestLeftNT)),2)/rightNegExp;
        return new DataSplit(index,mdpt,Arrays.copyOfRange(data,0,bestRight),Arrays.copyOfRange(data,bestRight,data.length),bestLeftNT,nTrue,bestGain,bestRatio,chiSquare);
    }

    protected static double entropy(double p) {
        if ( Double.compare(p,0.0) == 0 ) {
            return 0.0;
        }

        return -p*log2(p) - (1-p)*log2(1 - p);
    }

    protected static final double LOG2 = Math.log10(2);

    protected static double log2(double p) {
        return Math.log10(p)/LOG2;
    }

    private int getRightIdx(VariantDatum[] data, int atIdx, int leftIdx) {
        int rightIdx = leftIdx+1;
        if ( useFinal ) {
            while (rightIdx<data.length && data[leftIdx].finalAnnotations[atIdx] == data[rightIdx].finalAnnotations[atIdx]) {
                ++rightIdx;
            }
        } else {
            while (rightIdx<data.length && data[leftIdx].initialAnnotations[atIdx] == data[rightIdx].initialAnnotations[atIdx]) {
                ++rightIdx;
            }
        }
        return rightIdx;
    }

    public double evaluateDatum(VariantDatum point) {
        return Math.log10(tree.getConditionalProbability(point))-Math.log10(1.0-tree.getConditionalProbability(point));
    }

    protected static List<DecisionTreeNode> leaves(DecisionTreeNode root) {
        ArrayList<DecisionTreeNode> flattened = new ArrayList<DecisionTreeNode>(256);
        addChildren(flattened,root);
        return flattened;
    }

    private static void addChildren(List<DecisionTreeNode> accum, DecisionTreeNode node) {
        if ( node.leftChild == null ) {
            accum.add(node);
        } else {
            addChildren(accum,node.leftChild);
            addChildren(accum,node.rightChild);
        }
    }

    protected class DataSplit {
        final int attribIdx;
        final double threshold;
        final VariantDatum[] left;
        final VariantDatum[] right;
        final double gainRatio;
        final double chiSquare;
        final int totalTrue;
        final int leftTrue;
        final double gain;

        public DataSplit(int idx, double thresh, VariantDatum[] leftAr, VariantDatum[] rightAr, int nTrueLeft, int nTrue,
                         double gain, double ratio, double chi) {
            attribIdx = idx;
            threshold = thresh;
            left = leftAr;
            right = rightAr;
            gainRatio = ratio;
            chiSquare = chi;
            totalTrue = nTrue;
            leftTrue = nTrueLeft;
            this.gain = gain;
        }
    }

    protected class DecisionTreeNode  {
        private int DATA_INDEX;
        private double DATA_THRESHOLD;
        private final double CONDITIONAL_PTRUE;
        private final int nTrue;
        private DecisionTreeNode leftChild = null;
        private DecisionTreeNode rightChild = null;
        private VariantDatum[] data;
        private boolean useFinal;

        public DecisionTreeNode(int nTrue, VariantDatum[] data, boolean uFinal) {
            CONDITIONAL_PTRUE = ((double)nTrue)/data.length;
            this.nTrue = nTrue;
            this.data = data;
            this.useFinal = uFinal;
        }

        public int getNTrue() { return nTrue; }

        public void split(int idx, double thresh, int leftNTrue, int rightNTrue, VariantDatum[] leftData, VariantDatum[] rightData) {
            DATA_INDEX = idx;
            DATA_THRESHOLD = thresh;
            leftChild = new DecisionTreeNode(leftNTrue,leftData,useFinal);
            rightChild = new DecisionTreeNode(rightNTrue,rightData,useFinal);
            data = null;
        }

        public double getConditionalProbability(VariantDatum datum) {
            if ( leftChild == null )
                return CONDITIONAL_PTRUE;
            double[] annot = useFinal ? datum.finalAnnotations : datum.initialAnnotations;
            if ( annot[DATA_INDEX] < DATA_THRESHOLD )
                return leftChild.getConditionalProbability(datum);
            return rightChild.getConditionalProbability(datum);
        }

        public VariantDatum[] getData() {
            return data;
        }

        public void remove() { throw new IllegalArgumentException("Cannot use remove from iterator interface"); }
        public boolean canSplit() {
            if ( leftChild != null ) {
                throw new IllegalStateException("Should never be trying to split an internal node");
            }

            // todo -- this is a dumb heuristic, but <20 points suggests no power
            return data.length > 50;
        }
    }
}

