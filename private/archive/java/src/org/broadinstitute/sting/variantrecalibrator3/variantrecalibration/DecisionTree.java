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

