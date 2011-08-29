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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.Arrays;
import java.util.List;

public class LikelihoodCalculationEngine {

    private static final int MATCH_OFFSET = 0;
    private static final int X_OFFSET = 1;
    private static final int Y_OFFSET = 2;

    private static final int DIAG = 0;
    private static final int UP = 1;
    private static final int LEFT = 2;

    private static final int DIAG_GOTO_M = 0;
    private static final int DIAG_GOTO_X = 1;
    private static final int DIAG_GOTO_Y = 2;

    private static final int UP_GOTO_M = 4;
    private static final int UP_GOTO_X = 5;
    private static final int UP_GOTO_Y = 6;

    private static final int LEFT_GOTO_M = 8;
    private static final int LEFT_GOTO_X = 9;
    private static final int LEFT_GOTO_Y = 10;

    private static final int[] ACTIONS_M = {DIAG_GOTO_M, DIAG_GOTO_X, DIAG_GOTO_Y};
    private static final int[] ACTIONS_X = {UP_GOTO_M, UP_GOTO_X, UP_GOTO_Y};
    private static final int[] ACTIONS_Y = {LEFT_GOTO_M, LEFT_GOTO_X, LEFT_GOTO_Y};


    private final double logGapOpenProbability;
    private final double logGapContinuationProbability;

    private boolean DEBUG = false;

    private static final int MAX_CACHED_QUAL = 93;

    private static final double baseMatchArray[];
    private static final double baseMismatchArray[];

    private final static double LOG_ONE_HALF;
    private final static double END_GAP_COST;

    private static final int START_HRUN_GAP_IDX = 4;
    private static final int MAX_HRUN_GAP_IDX = 20;

    private static final double MIN_GAP_OPEN_PENALTY = 3.0;
    private static final double MIN_GAP_CONT_PENALTY = 1.0;
    private static final double GAP_PENALTY_HRUN_STEP = 1.0; // each increase in hrun decreases gap penalty by this.


    private boolean doViterbi = false;

    private final boolean useAffineGapModel = true;
    private boolean doContextDependentPenalties = false;

    private final double[] GAP_OPEN_PROB_TABLE;
    private final double[] GAP_CONT_PROB_TABLE;

    private boolean getGapPenaltiesFromFile = false;

    static {
        LOG_ONE_HALF = -Math.log10(2.0);
        END_GAP_COST = LOG_ONE_HALF;

        baseMatchArray = new double[MAX_CACHED_QUAL+1];
        baseMismatchArray = new double[MAX_CACHED_QUAL+1];
        for (int k=1; k <= MAX_CACHED_QUAL; k++) {
            double baseProb = Math.pow(10, -k/10.);


            baseMatchArray[k] =  Math.log10(1-baseProb);
            baseMismatchArray[k] = Math.log10(baseProb);
        }
    }

    public LikelihoodCalculationEngine(double indelGOP, double indelGCP, boolean deb, boolean doCDP, boolean dovit) {
        this(indelGOP, indelGCP, deb, doCDP);
        this.doViterbi = dovit;
    }

    public LikelihoodCalculationEngine(double indelGOP, double indelGCP, boolean deb, boolean doCDP) {

        this.logGapOpenProbability = -indelGOP/10.0; // QUAL to log prob
        this.logGapContinuationProbability = -indelGCP/10.0; // QUAL to log prob
        this.doContextDependentPenalties = doCDP;
        this.DEBUG = deb;


        // fill gap penalty table, affine naive model:
        this.GAP_CONT_PROB_TABLE = new double[MAX_HRUN_GAP_IDX];
        this.GAP_OPEN_PROB_TABLE = new double[MAX_HRUN_GAP_IDX];

        for (int i = 0; i < START_HRUN_GAP_IDX; i++) {
            GAP_OPEN_PROB_TABLE[i] = logGapOpenProbability;
            GAP_CONT_PROB_TABLE[i] = logGapContinuationProbability;
        }

        double gop = logGapOpenProbability;
        double gcp = logGapContinuationProbability;
        double step = GAP_PENALTY_HRUN_STEP/10.0;

        double maxGOP = -MIN_GAP_OPEN_PENALTY/10.0;  // phred to log prob
        double maxGCP = -MIN_GAP_CONT_PENALTY/10.0;  // phred to log prob

        for (int i=START_HRUN_GAP_IDX; i < MAX_HRUN_GAP_IDX; i++) {
            gop += step;
            if (gop > maxGOP)
                gop = maxGOP;

            gcp += step;
            if(gcp > maxGCP)
                gcp = maxGCP;
            GAP_OPEN_PROB_TABLE[i] = gop;
            GAP_CONT_PROB_TABLE[i] = gcp;
        }

    }

    public Pair<Haplotype, Haplotype> computeLikelihoods(final List<Haplotype> haplotypes, final List<SAMRecord> reads) {
        int numHaplotypes = haplotypes.size();
        double haplotypeLikehoodMatrix[][] = new double[numHaplotypes][numHaplotypes];
        double readLikelihoods[][] = new double[reads.size()][numHaplotypes];
        double gop[] = new double[haplotypes.get(0).bases.length];
        double gcp[] = new double[haplotypes.get(0).bases.length];

        Arrays.fill(gop, logGapOpenProbability); // this should eventually be derived from the data
        Arrays.fill(gcp, logGapContinuationProbability); // this should eventually be derived from the data

        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
                haplotypeLikehoodMatrix[iii][jjj] = 0.0;
            }
        }

        for( int iii = 0; iii < reads.size(); iii++ ) {
            final SAMRecord read = reads.get(iii);
            for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
                final Haplotype haplotype = haplotypes.get(jjj);
                readLikelihoods[iii][jjj] = computeReadLikelihoodGivenHaplotypeAffineGaps(haplotype.bases, read.getReadBases(), read.getBaseQualities(), gop, gcp);
            }
        }

        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = iii; jjj < numHaplotypes; jjj++ ) {
                for( int kkk = 0; kkk < reads.size(); kkk++ ) {

                    // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                    // First term is approximated by Jacobian log with table lookup.
                    if (Double.isInfinite(readLikelihoods[kkk][iii]) && Double.isInfinite(readLikelihoods[kkk][jjj])) {
                        continue;
                    }
                    haplotypeLikehoodMatrix[iii][jjj] += ( MathUtils.softMax(readLikelihoods[kkk][iii], readLikelihoods[kkk][jjj]) + LOG_ONE_HALF );
                }
            }
        }

        double maxElement = Double.NEGATIVE_INFINITY;
        int hap1 = -1;
        int hap2 = -1;
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = iii; jjj < numHaplotypes; jjj++ ) {
                if( haplotypeLikehoodMatrix[iii][jjj] > maxElement ) {
                    maxElement = haplotypeLikehoodMatrix[iii][jjj];
                    hap1 = iii;
                    hap2 = jjj;
                }
            }
        }
        return new Pair<Haplotype, Haplotype>(haplotypes.get(hap1), haplotypes.get(hap2));
    }

    private static double[] getHaplotypeLikelihoods(double[][] haplotypeLikehoodMatrix) {
        int hSize = haplotypeLikehoodMatrix.length;
        double[] genotypeLikelihoods = new double[hSize*(hSize+1)/2];

        int k=0;
        double maxElement = Double.NEGATIVE_INFINITY;
        for (int j=0; j < hSize; j++) {
            for (int i=0; i <= j; i++){
                genotypeLikelihoods[k++] = haplotypeLikehoodMatrix[i][j];
                if (haplotypeLikehoodMatrix[i][j] > maxElement)
                    maxElement = haplotypeLikehoodMatrix[i][j];
            }
        }

        // renormalize
        for (int i=0; i < genotypeLikelihoods.length; i++)
            genotypeLikelihoods[i] -= maxElement;

        return genotypeLikelihoods;
    }

    private double computeReadLikelihoodGivenHaplotypeAffineGaps(byte[] haplotypeBases, byte[] readBases, byte[] readQuals,
                                                                 double[] currentGOP, double[] currentGCP) {

        final int X_METRIC_LENGTH = readBases.length+1;
        final int Y_METRIC_LENGTH = haplotypeBases.length+1;

        // initialize path metric and traceback memories for likelihood computation
        double[][] matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        double[][] XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        double[][] YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        int[][] bestActionArrayM = new int[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        int[][] bestActionArrayX = new int[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        int[][] bestActionArrayY = new int[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        double c,d;
        matchMetricArray[0][0]= END_GAP_COST;//Double.NEGATIVE_INFINITY;

        for (int i=1; i < X_METRIC_LENGTH; i++) {
            //initialize first column
            matchMetricArray[i][0]  = Double.NEGATIVE_INFINITY;
            YMetricArray[i][0]      = Double.NEGATIVE_INFINITY;
            XMetricArray[i][0]      = END_GAP_COST*(i);//logGapOpenProbability + (i-1)*logGapContinuationProbability;

            bestActionArrayX[i][0] = bestActionArrayY[i][0] = bestActionArrayM[i][0] = UP_GOTO_X;
        }

        for (int j=1; j < Y_METRIC_LENGTH; j++) {
            // initialize first row
            matchMetricArray[0][j]  = Double.NEGATIVE_INFINITY;
            XMetricArray[0][j]      = Double.NEGATIVE_INFINITY;
            YMetricArray[0][j]      = END_GAP_COST*(j);//logGapOpenProbability + (j-1) * logGapContinuationProbability;

            bestActionArrayY[0][j] = bestActionArrayM[0][j] = bestActionArrayX[0][j] = LEFT_GOTO_Y;
        }

        for (int indI=1; indI < X_METRIC_LENGTH; indI++) {
            int im1 = indI-1;
            for (int indJ=1; indJ < Y_METRIC_LENGTH; indJ++) {
                int jm1 = indJ-1;
                byte x = readBases[im1];
                byte y = haplotypeBases[jm1];
                byte qual = readQuals[im1];

                double bestMetric = 0.0;
                int bestMetricIdx = 0;

                // compute metric for match/mismatch
                // workaround for reads whose bases quality = 0,
                if (qual < 1)
                    qual = 1;

                if (qual > MAX_CACHED_QUAL)
                    qual = MAX_CACHED_QUAL;

                double pBaseRead =  (x == y)? baseMatchArray[(int)qual]:baseMismatchArray[(int)qual];


                double[] metrics = new double[3];


                if (doViterbi) {
                    // update match array
                    metrics[MATCH_OFFSET] = matchMetricArray[im1][jm1] + pBaseRead;
                    metrics[X_OFFSET] = XMetricArray[im1][jm1] + pBaseRead;
                    metrics[Y_OFFSET] = YMetricArray[im1][jm1] + pBaseRead;

                    bestMetricIdx = MathUtils.maxElementIndex(metrics);
                    bestMetric = metrics[bestMetricIdx];
                }
                else
                    bestMetric = MathUtils.softMax(matchMetricArray[im1][jm1] + pBaseRead, XMetricArray[im1][jm1] + pBaseRead,
                            YMetricArray[im1][jm1] + pBaseRead);

                matchMetricArray[indI][indJ] = bestMetric;
                bestActionArrayM[indI][indJ] = ACTIONS_M[bestMetricIdx];

                // update X array
                // State X(i,j): X(1:i) aligned to a gap in Y(1:j).
                // When in last column of X, ie X(1:i) aligned to full Y, we don't want to penalize gaps

                //c = (indJ==Y_METRIC_LENGTH-1? END_GAP_COST: currentGOP[jm1]);
                //d = (indJ==Y_METRIC_LENGTH-1? END_GAP_COST: currentGCP[jm1]);
                if (getGapPenaltiesFromFile) {
                    c = currentGOP[im1];
                    d = logGapContinuationProbability;

                } else {
                    c = currentGOP[jm1];
                    d = currentGCP[jm1];
                }

                if (indJ == Y_METRIC_LENGTH-1)
                    c = d = END_GAP_COST;

                if (doViterbi) {
                    metrics[MATCH_OFFSET] = matchMetricArray[im1][indJ] + c;
                    metrics[X_OFFSET] = XMetricArray[im1][indJ] + d;
                    metrics[Y_OFFSET] = Double.NEGATIVE_INFINITY; //YMetricArray[indI-1][indJ] + logGapOpenProbability;

                    bestMetricIdx = MathUtils.maxElementIndex(metrics);
                    bestMetric = metrics[bestMetricIdx];
                }
                else
                    bestMetric = MathUtils.softMax(matchMetricArray[im1][indJ] + c, XMetricArray[im1][indJ] + d);

                XMetricArray[indI][indJ] = bestMetric;
                bestActionArrayX[indI][indJ] = ACTIONS_X[bestMetricIdx];

                // update Y array
                //c = (indI==X_METRIC_LENGTH-1? END_GAP_COST: currentGOP[jm1]);
                //d = (indI==X_METRIC_LENGTH-1? END_GAP_COST: currentGCP[jm1]);
                if (getGapPenaltiesFromFile) {
                    c = currentGOP[im1];
                    d = logGapContinuationProbability;
                }
                else {
                    c = currentGOP[jm1];
                    d = currentGCP[jm1];
                }

                if (indI == X_METRIC_LENGTH-1) {
                    c = d = END_GAP_COST;
                }

                if (doViterbi) {
                    metrics[MATCH_OFFSET] = matchMetricArray[indI][jm1] + c;
                    metrics[X_OFFSET] = Double.NEGATIVE_INFINITY; //XMetricArray[indI][indJ-1] + logGapOpenProbability;
                    metrics[Y_OFFSET] = YMetricArray[indI][jm1] + d;

                    bestMetricIdx = MathUtils.maxElementIndex(metrics);
                    bestMetric = metrics[bestMetricIdx];
                }
                else {
                    bestMetric = MathUtils.softMax(matchMetricArray[indI][jm1] + c, YMetricArray[indI][jm1] + d);
                }

                YMetricArray[indI][indJ] = bestMetric;
                bestActionArrayY[indI][indJ] = ACTIONS_Y[bestMetricIdx];

            }
        }

        double bestMetric;
        double metrics[] = new double[3];
        int bestTable=0, bestI=X_METRIC_LENGTH - 1, bestJ=Y_METRIC_LENGTH - 1;
        metrics[MATCH_OFFSET] = matchMetricArray[bestI][bestJ];
        metrics[X_OFFSET] = XMetricArray[bestI][bestJ];
        metrics[Y_OFFSET] = YMetricArray[bestI][bestJ];
        if (doViterbi) {
            bestTable = MathUtils.maxElementIndex(metrics);
            bestMetric = metrics[bestTable];
        }
        else {
            bestMetric = MathUtils.softMax(metrics);
        }

        // Do traceback (needed only for debugging!)
        if (DEBUG && doViterbi) {

            int bestAction;
            int i = bestI;
            int j = bestJ;


            System.out.println("Affine gap NW");


            String haplotypeString = new String (haplotypeBases);
            String readString = new String(readBases);


            while (i >0 || j >0) {
                if (bestTable == X_OFFSET) {
                    // insert gap in Y
                    haplotypeString = haplotypeString.substring(0,j)+"-"+haplotypeString.substring(j);
                    bestAction = bestActionArrayX[i][j];
                }
                else if (bestTable == Y_OFFSET) {
                    readString = readString.substring(0,i)+"-"+readString.substring(i);
                    bestAction = bestActionArrayY[i][j];

                }
                else {
                    bestAction = bestActionArrayM[i][j];
                }
                System.out.print(bestAction);


                // bestAction contains action to take at next step
                // encoding of bestAction: upper 2 bits = direction, lower 2 bits = next table

                // bestTable and nextDirection for next step
                bestTable = bestAction & 0x3;
                int nextDirection = bestAction >> 2;
                if (nextDirection == UP) {
                    i--;
                } else if (nextDirection == LEFT) {
                    j--;
                } else { //  if (nextDirection == DIAG)
                    i--; j--;
                }

            }




            System.out.println("\nAlignment: ");
            System.out.println("R:"+readString);
            System.out.println("H:"+haplotypeString);
            System.out.println();


        }
        if (DEBUG)
            System.out.format("Likelihood: %5.4f\n", bestMetric);

        return bestMetric;

    }

}