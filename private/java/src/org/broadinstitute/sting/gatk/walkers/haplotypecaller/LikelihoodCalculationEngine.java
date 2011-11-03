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
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.*;

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

    private static final double MIN_GAP_OPEN_PENALTY = 30.0;
    private static final double MIN_GAP_CONT_PENALTY = 10.0;
    private static final double GAP_PENALTY_HRUN_STEP = 1.0; // each increase in hrun decreases gap penalty by this.


    private boolean doViterbi = false;

    private final boolean useAffineGapModel = true;
    private boolean doContextDependentPenalties = false;

    private final double[] GAP_OPEN_PROB_TABLE;
    private final double[] GAP_CONT_PROB_TABLE;

    private boolean getGapPenaltiesFromFile = false;

    private final NestedHashMap kmerQualityTables;
    private final int CONTEXT_SIZE;

    public Double haplotypeLikehoodMatrix[][];

    static {
        LOG_ONE_HALF = -Math.log10(2.0);
        END_GAP_COST = LOG_ONE_HALF;

        baseMatchArray = new double[MAX_CACHED_QUAL+1];
        baseMismatchArray = new double[MAX_CACHED_QUAL+1];
        for (int k=1; k <= MAX_CACHED_QUAL; k++) {
            double baseProb = Math.pow(10.0, ((double) -k)/10.0);

            baseMatchArray[k] =  Math.log10(1.0-baseProb);
            baseMismatchArray[k] = Math.log10(baseProb);
        }
    }

    public LikelihoodCalculationEngine( double indelGOP, double indelGCP, boolean deb, boolean doCDP, boolean dovit, final NestedHashMap kmerQualityTables, final int contextSize ) {
        this(indelGOP, indelGCP, deb, doCDP, kmerQualityTables, contextSize);
        this.doViterbi = dovit;
    }

    public LikelihoodCalculationEngine( double indelGOP, double indelGCP, boolean deb, boolean doCDP, final NestedHashMap kmerQualityTables, final int contextSize ) {

        this.kmerQualityTables = kmerQualityTables;
        this.CONTEXT_SIZE = contextSize;

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

    public void computeLikelihoods( final ArrayList<Haplotype> haplotypes, final ArrayList<SAMRecord> reads ) {
        int numHaplotypes = haplotypes.size();
        double readLikelihoods[][] = new double[reads.size()][numHaplotypes];
        haplotypeLikehoodMatrix = new Double[numHaplotypes][numHaplotypes];

        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
                haplotypeLikehoodMatrix[iii][jjj] = 0.0;
            }
        }

        int maxHaplotypeLength = 0;
        for( final Haplotype h : haplotypes ) {
            int length = h.bases.length;
            if(length > maxHaplotypeLength) { maxHaplotypeLength = length; }
        }

        for( int iii = 0; iii < reads.size(); iii++ ) {
            final SAMRecord read = reads.get(iii);
            final String readGroup = read.getReadGroup().getReadGroupId();

            // initialize path metric and traceback memories for likelihood computation
            double[][] matchMetricArray = null, XMetricArray = null, YMetricArray = null;
            byte[] previousHaplotypeSeen = null;
            double[] previousGOP = null;
            int startIdx;

            for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
                final Haplotype haplotype = haplotypes.get(jjj);
                haplotype.extendHaplotype( maxHaplotypeLength );
                final byte[] haplotypeBases = haplotype.extendedBases;
                final double[] contextLogGapOpenProbabilities = new double[haplotypeBases.length];
                final double[] contextLogGapContinuationProbabilities = new double[haplotypeBases.length];

                if (matchMetricArray == null) {
                    final int X_METRIC_LENGTH = read.getReadLength()+1;
                    final int Y_METRIC_LENGTH = haplotypeBases.length+1;

                    matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
                    XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
                    YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
                }

                //Arrays.fill(contextLogGapOpenProbabilities, logGapOpenProbability); // this should eventually be derived from the data
                fillGapProbabilitiesFromQualityTables( readGroup, haplotypeBases, contextLogGapOpenProbabilities );
                Arrays.fill(contextLogGapContinuationProbabilities, logGapContinuationProbability); // this should eventually be derived from the data

                if (previousHaplotypeSeen == null) {
                    startIdx = 0;
                } else {
                    int s1 = computeFirstDifferingPosition(haplotypeBases, previousHaplotypeSeen);
                    int s2 = computeFirstDifferingPosition(contextLogGapOpenProbabilities, previousGOP);
                    startIdx = Math.min(s1,s2);
                }
                previousHaplotypeSeen = haplotypeBases.clone();
                previousGOP = contextLogGapOpenProbabilities.clone();


                readLikelihoods[iii][jjj] = computeReadLikelihoodGivenHaplotypeAffineGaps(haplotypeBases, read.getReadBases(), read.getBaseQualities(),
                        contextLogGapOpenProbabilities, contextLogGapContinuationProbabilities, startIdx, matchMetricArray, XMetricArray, YMetricArray);
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

                    final SAMRecord read = reads.get(kkk);
                    final int mappingLength = read.getAlignmentEnd() - read.getAlignmentStart() + 1;
                    final double mappingProb = 1.0 - Math.max(0.0, (76.0 - ((double)mappingLength)) / 76.0); //BUGBUG: 101!, needs to pull from the empirical read length distribution per read group

                    haplotypeLikehoodMatrix[iii][jjj] += (mappingProb*mappingProb) * ( MathUtils.softMax(readLikelihoods[kkk][iii], readLikelihoods[kkk][jjj]) + LOG_ONE_HALF ); // BUGBUG: needs to be a logged probability
                }
            }
        }

        double[] genotypeLikelihoods = new double[numHaplotypes*(numHaplotypes+1)/2];
        int k=0;
        for (int j=0; j < numHaplotypes; j++) {
            for (int i=0; i <= j; i++){
                genotypeLikelihoods[k++] = haplotypeLikehoodMatrix[i][j];
            }
        }

        genotypeLikelihoods = MathUtils.normalizeFromLog10(genotypeLikelihoods, false, true);

        k=0;
        for (int j=0; j < numHaplotypes; j++) {
            for (int i=0; i <= j; i++){
                haplotypeLikehoodMatrix[i][j] = genotypeLikelihoods[k++];
            }
        }


        for( int iii = 1; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj < iii; jjj++ ) {
                haplotypeLikehoodMatrix[iii][jjj] = haplotypeLikehoodMatrix[jjj][iii]; // fill in the symmetric lower triangular part of the matrix for convenience later
            }
        }
    }

    public Set<Haplotype> chooseBestHaplotypes( final ArrayList<Haplotype> haplotypes ) {

        // For now we choose the top two haplotypes by finding the max value of the pairwise matrix
        // in the future we could use AIC or some other criterion to select more haplotypes to best explain the read data

        final int numHaplotypes = haplotypes.size();
        final HashSet<Haplotype> returnHaplotypeSet = new HashSet<Haplotype>();
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

        returnHaplotypeSet.add(haplotypes.get(hap1));
        returnHaplotypeSet.add(haplotypes.get(hap2));

        return returnHaplotypeSet;
    }

    private void fillGapProbabilitiesFromQualityTables( final String readGroup, final byte[] refBytes, final double[] contextLogGapOpenProbabilities ) {

        final Object[] key = new Object[2];
        key[0] = readGroup;
        for(int i = 0; i < refBytes.length; i++) {

            Double gop = null;
            final String bases = ( i-CONTEXT_SIZE < 0 ? null : new String(Arrays.copyOfRange(refBytes,i-CONTEXT_SIZE,i)) );
            if( bases != null ) {
                key[1] = bases;
                gop = (Double) kmerQualityTables.get( key );
            }
            contextLogGapOpenProbabilities[i] = ( gop != null ? gop : -4.5 );
        }
    }

    private int computeFirstDifferingPosition(byte[] b1, byte[] b2) {
        if (b1.length != b2.length)
            return 0; // sanity check

        for (int i=0; i < b1.length; i++ ){
            if ( b1[i] != b2[i] )
                return i;
        }
        return 0; // sanity check
    }

    private int computeFirstDifferingPosition(double[] b1, double[] b2) {
        if (b1.length != b2.length)
            return 0; // sanity check

        for (int i=0; i < b1.length; i++ ){
            if ( b1[i] != b2[i] )
                return i;
        }
        return 0; // sanity check
    }

    private void updateCell(final int indI, final int indJ, final int X_METRIC_LENGTH, final int Y_METRIC_LENGTH, byte[] readBases, byte[] readQuals, byte[] haplotypeBases,
                            double[] currentGOP, double[] currentGCP,  double[][] matchMetricArray,  double[][] XMetricArray,  double[][] YMetricArray) {
        if (indI > 0 && indJ > 0) {
            final int im1 = indI -1;
            final int jm1 = indJ - 1;
            // update current point
            final byte x = readBases[im1];
            final byte y = haplotypeBases[jm1];
            final byte qual = readQuals[im1] < 1 ? 1 : (readQuals[im1] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : readQuals[im1]);

            final double pBaseRead =  (x == y)? baseMatchArray[(int)qual]:baseMismatchArray[(int)qual];

            matchMetricArray[indI][indJ] = MathUtils.softMax(matchMetricArray[im1][jm1] + pBaseRead, XMetricArray[im1][jm1] + pBaseRead,
                    YMetricArray[im1][jm1] + pBaseRead);

            final double c1 = indJ == Y_METRIC_LENGTH-1 ? END_GAP_COST : currentGOP[jm1];
            final double d1 = indJ == Y_METRIC_LENGTH-1 ? END_GAP_COST : currentGCP[jm1];

            XMetricArray[indI][indJ] = MathUtils.softMax(matchMetricArray[im1][indJ] + c1, XMetricArray[im1][indJ] + d1);

            // update Y array
            final double c2 = indI == X_METRIC_LENGTH-1 ? END_GAP_COST : currentGOP[jm1];
            final double d2 = indI == X_METRIC_LENGTH-1 ? END_GAP_COST : currentGCP[jm1];
            YMetricArray[indI][indJ] = MathUtils.softMax(matchMetricArray[indI][jm1] + c2, YMetricArray[indI][jm1] + d2);
        }
    }

    private double computeReadLikelihoodGivenHaplotypeAffineGaps(byte[] haplotypeBases, byte[] readBases, byte[] readQuals,
                                                                 double[] currentGOP, double[] currentGCP, int indToStart,
                                                                 double[][] matchMetricArray, double[][] XMetricArray, double[][] YMetricArray) {

        final boolean bandedLikelihoods = false;

        final int X_METRIC_LENGTH = readBases.length+1;
        final int Y_METRIC_LENGTH = haplotypeBases.length+1;

        if (indToStart == 0) {
            // default initialization for all arrays

            for (int i=0; i < X_METRIC_LENGTH; i++) {
                Arrays.fill(matchMetricArray[i],Double.NEGATIVE_INFINITY);
                Arrays.fill(YMetricArray[i],Double.NEGATIVE_INFINITY);
                Arrays.fill(XMetricArray[i],Double.NEGATIVE_INFINITY);
            }

            for (int i=1; i < X_METRIC_LENGTH; i++) {
                //initialize first column
                XMetricArray[i][0]      = END_GAP_COST*(i);
            }

            for (int j=1; j < Y_METRIC_LENGTH; j++) {
                // initialize first row
                YMetricArray[0][j]      = END_GAP_COST*(j);
            }
            matchMetricArray[0][0]= END_GAP_COST;//Double.NEGATIVE_INFINITY;
            XMetricArray[0][0]=  YMetricArray[0][0] = 0;
        }


        if (bandedLikelihoods) {
            final double DIAG_TOL = 40; // means that max - min element in diags have to be > this number for banding to take effect.

            final int numDiags = X_METRIC_LENGTH +  Y_METRIC_LENGTH -1;
            final int elemsInDiag = Math.min(X_METRIC_LENGTH, Y_METRIC_LENGTH);

            int idxWithMaxElement = 0;

            for (int  diag=indToStart; diag <  numDiags; diag++) {
                // compute default I and J start positions at edge of diagonals
                int indI = 0;
                int indJ = diag;
                if (diag >= Y_METRIC_LENGTH ) {
                    indI = diag-(Y_METRIC_LENGTH-1);
                    indJ = Y_METRIC_LENGTH-1;
                }

                // first pass: from max element to edge
                int idxLow =  idxWithMaxElement;

                // reset diag max value before starting
                double maxElementInDiag = Double.NEGATIVE_INFINITY;
                // set indI, indJ to correct values
                indI += idxLow;
                indJ -= idxLow;
                if (indI >= X_METRIC_LENGTH || indJ < 0) {
                    idxLow--;
                    indI--;
                    indJ++;
                }


                for (int el = idxLow; el < elemsInDiag; el++) {
                    updateCell(indI, indJ, X_METRIC_LENGTH, Y_METRIC_LENGTH, readBases, readQuals, haplotypeBases,
                            currentGOP, currentGCP,  matchMetricArray,  XMetricArray, YMetricArray);
                    // update max in diagonal
                    final double bestMetric = MathUtils.max(matchMetricArray[indI][indJ], XMetricArray[indI][indJ], YMetricArray[indI][indJ]);

                    // check if we've fallen off diagonal value by threshold
                    if (bestMetric > maxElementInDiag) {
                        maxElementInDiag = bestMetric;
                        idxWithMaxElement = el;
                    }
                    else if (bestMetric < maxElementInDiag - DIAG_TOL && idxWithMaxElement > 0)
                        break; // done w/current diagonal

                    indI++;
                    if (indI >=X_METRIC_LENGTH )
                        break;
                    indJ--;
                    if (indJ <= 0)
                        break;
                }
                if (idxLow > 0) {
                    // now do second part in opposite direction
                    indI = 0;
                    indJ = diag;
                    if (diag >= Y_METRIC_LENGTH ) {
                        indI = diag-(Y_METRIC_LENGTH-1);
                        indJ = Y_METRIC_LENGTH-1;
                    }

                    indI += idxLow-1;
                    indJ -= idxLow-1;
                    for (int el = idxLow-1; el >= 0; el--) {

                        updateCell(indI, indJ, X_METRIC_LENGTH, Y_METRIC_LENGTH, readBases, readQuals, haplotypeBases,
                                currentGOP, currentGCP,  matchMetricArray,  XMetricArray, YMetricArray);
                        // update max in diagonal
                        final double bestMetric = MathUtils.max(matchMetricArray[indI][indJ], XMetricArray[indI][indJ], YMetricArray[indI][indJ]);

                        // check if we've fallen off diagonal value by threshold
                        if (bestMetric > maxElementInDiag) {
                            maxElementInDiag = bestMetric;
                            idxWithMaxElement = el;
                        }
                        else if (bestMetric < maxElementInDiag - DIAG_TOL)
                            break; // done w/current diagonal

                        indJ++;
                        if (indJ >= Y_METRIC_LENGTH )
                            break;
                        indI--;
                        if (indI <= 0)
                            break;
                    }
                }
            }
        }
        else {
            // simplified rectangular version of update loop
            for (int indI=1; indI < X_METRIC_LENGTH; indI++) {
                for (int indJ=indToStart+1; indJ < Y_METRIC_LENGTH; indJ++) {
                    updateCell(indI, indJ, X_METRIC_LENGTH, Y_METRIC_LENGTH, readBases, readQuals, haplotypeBases,
                            currentGOP, currentGCP,  matchMetricArray,  XMetricArray, YMetricArray);

                }
            }
        }

        final int bestI = X_METRIC_LENGTH - 1, bestJ = Y_METRIC_LENGTH - 1;
        return MathUtils.softMax(matchMetricArray[bestI][bestJ], XMetricArray[bestI][bestJ], YMetricArray[bestI][bestJ]);
    }
}