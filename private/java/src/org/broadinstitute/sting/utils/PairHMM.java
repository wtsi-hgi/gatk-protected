/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.SAMUtils;

import java.util.Arrays;

/**
 * Util class for performing the pair HMM for local alignment. Figure 4.3 in Durbin 1998 book.
 * User: rpoplin
 * Date: 3/1/12
 */

public class PairHMM {
    private static final int MAX_CACHED_QUAL = (int)Byte.MAX_VALUE;
    private static final double matchPenalty[];
    private static final double mismatchPenalty[];
    private static final byte DEFAULT_GOP = (byte) 45;
    private static final byte DEFAULT_GCP = (byte) 10;

    static {
        matchPenalty = new double[MAX_CACHED_QUAL+1];
        mismatchPenalty = new double[MAX_CACHED_QUAL+1];
        for (int k=1; k <= MAX_CACHED_QUAL; k++) {
            double baseProb = Math.pow(10.0, ((double) -k)/10.0);

            matchPenalty[k] =  Math.log10(1.0-baseProb);
            mismatchPenalty[k] = Math.log10(baseProb);
        }
    }

    @Requires({"readBases.length == readQuals.length","readBases.length == insertionGOP.length","readBases.length == deletionGOP.length","readBases.length == overallGCP.length"})
    @Ensures({"result <= 0.0"})
    public double computeReadLikelihoodGivenHaplotype( final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals,
                                                       final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP ) {
        
        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions
        final int X_METRIC_LENGTH = readBases.length + 1;
        final int Y_METRIC_LENGTH = haplotypeBases.length + 1;
        
        // initial arrays to hold the probabilities of being in the match, insertion and deletion cases
        final double[][] matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        final double[][] XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        final double[][] YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        
        for( int iii=0; iii < X_METRIC_LENGTH; iii++ ) {
            Arrays.fill(matchMetricArray[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(XMetricArray[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(YMetricArray[iii], Double.NEGATIVE_INFINITY);
        }
        
        // the initial condition
        matchMetricArray[1][1] = 0.0; // Math.log10(1.0);

        // simplified rectangular version of update loop
        for( int iii=1; iii < X_METRIC_LENGTH; iii++ ) {
            for( int jjj=1; jjj < Y_METRIC_LENGTH; jjj++ ) {
                if( iii == 1 && jjj == 1 ) { continue; }
                updateCell(iii, jjj, haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP,
                           matchMetricArray, XMetricArray, YMetricArray);

            }
        }

        // final probability is the log10 sum of the last element in all three state arrays
        final int endI = X_METRIC_LENGTH - 1;
        final int endJ = Y_METRIC_LENGTH - 1;
        
        /*
        if(MathUtils.approximateLog10SumLog10(matchMetricArray[endI][endJ], XMetricArray[endI][endJ], YMetricArray[endI][endJ]) > 0.0) {
            System.out.println("ERROR! + " + String.format("%f,%f,%f",matchMetricArray[endI][endJ], XMetricArray[endI][endJ], YMetricArray[endI][endJ]) + " = " +MathUtils.approximateLog10SumLog10(matchMetricArray[endI][endJ], XMetricArray[endI][endJ], YMetricArray[endI][endJ]));
            double val = 0.0;
            for( int iii = 0; iii < 780009; iii++) {
                double x = 0.0;
                val += x;
            }
            System.out.println(val);
            System.out.println(new String(haplotypeBases));
            System.out.println(new String(readBases));
            System.out.println(SAMUtils.phredToFastq(readQuals));
            System.out.println(SAMUtils.phredToFastq(insertionGOP));
            System.out.println(SAMUtils.phredToFastq(deletionGOP));
            System.out.println(SAMUtils.phredToFastq(overallGCP));
            for( int iii=1; iii < X_METRIC_LENGTH; iii++ ) {
                for( int jjj=1; jjj < Y_METRIC_LENGTH; jjj++ ) {
                    System.out.print(String.format("%.2f,", matchMetricArray[iii][jjj]));
                }
                System.out.println();
            }
            System.out.println();System.out.println();System.out.println();System.out.println();
            for( int iii=1; iii < X_METRIC_LENGTH; iii++ ) {
                for( int jjj=1; jjj < Y_METRIC_LENGTH; jjj++ ) {
                    System.out.print(String.format("%.2f,", XMetricArray[iii][jjj]));
                }
                System.out.println();
            }
            System.out.println();System.out.println();System.out.println();System.out.println();
            for( int iii=1; iii < X_METRIC_LENGTH; iii++ ) {
                for( int jjj=1; jjj < Y_METRIC_LENGTH; jjj++ ) {
                    System.out.print(String.format("%.2f,", YMetricArray[iii][jjj]));
                }
                System.out.println();
            }
            System.out.println();
        }
        */
        return MathUtils.approximateLog10SumLog10(new double[]{matchMetricArray[endI][endJ], XMetricArray[endI][endJ], YMetricArray[endI][endJ]});
    }

    private void updateCell( final int indI, final int indJ, final byte[] haplotypeBases, final byte[] readBases,
                             final byte[] readQuals, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP,
                             final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray ) {

        // the read and haplotype indicies are offset by one because the state arrays have an extra column to hold the initial conditions
        final int im1 = indI - 1;
        final int jm1 = indJ - 1;

        // update the match array
        double pBaseReadLog10 = 0.0; // Math.log10(1.0);
        if( im1 > 0 && jm1 > 0 ) { // the emission probability is applied when leaving the state
            final byte x = readBases[im1-1];
            final byte y = haplotypeBases[jm1-1];
            final byte qual = ( readQuals[im1-1] < QualityUtils.MIN_USABLE_Q_SCORE ? QualityUtils.MIN_USABLE_Q_SCORE : (readQuals[im1-1] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : readQuals[im1-1]) );
            pBaseReadLog10 = ( x == y || x == (byte) 'N' || y == (byte) 'N' ? matchPenalty[(int)qual] : mismatchPenalty[(int)qual] );
        }
        final int qualIndexGOP = ( im1 == 0 ? DEFAULT_GOP + DEFAULT_GOP : ( insertionGOP[im1-1] + deletionGOP[im1-1] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : insertionGOP[im1-1] + deletionGOP[im1-1]) );
        final double d0 = matchPenalty[qualIndexGOP];
        final double e0 = ( im1 == 0 ? matchPenalty[(int)DEFAULT_GCP] : matchPenalty[(int)overallGCP[im1-1]] );
        matchMetricArray[indI][indJ] = pBaseReadLog10 + MathUtils.approximateLog10SumLog10(
                new double[]{matchMetricArray[indI-1][indJ-1] + d0, XMetricArray[indI-1][indJ-1] + e0, YMetricArray[indI-1][indJ-1] + e0});

        // update the X (insertion) array
        final double d1 = ( im1 == 0 ? mismatchPenalty[(int)DEFAULT_GOP] : mismatchPenalty[(int)insertionGOP[im1-1]] );
        final double e1 = ( im1 == 0 ? mismatchPenalty[(int)DEFAULT_GCP] : mismatchPenalty[(int)overallGCP[im1-1]] );
        final double qBaseReadLog10 = 0.0; // Math.log10(1.0) -- we don't have an estimate for this emission probability so assume q=1.0
        XMetricArray[indI][indJ] = qBaseReadLog10 + MathUtils.approximateLog10SumLog10(matchMetricArray[indI-1][indJ] + d1, XMetricArray[indI-1][indJ] + e1);

        // update the Y (deletion) array, with penalty of zero on the left and right flanks to allow for a local alignment within the haplotype
        final double d2 = ( im1 == 0 || im1 == readBases.length - 1 ? 0.0 : mismatchPenalty[(int)deletionGOP[im1-1]] );
        final double e2 = ( im1 == 0 || im1 == readBases.length - 1 ? 0.0 : mismatchPenalty[(int)overallGCP[im1-1]] );
        final double qBaseRefLog10 = 0.0; // Math.log10(1.0) -- we don't have an estimate for this emission probability so assume q=1.0
        YMetricArray[indI][indJ] = qBaseRefLog10 + MathUtils.approximateLog10SumLog10(matchMetricArray[indI][indJ-1] + d2, YMetricArray[indI][indJ-1] + e2);
    }
}

