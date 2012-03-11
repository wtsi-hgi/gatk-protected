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

import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.PairHMM;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

public class LikelihoodCalculationEngine {

    private final static double LOG_ONE_HALF;
    private final byte constantGCP;
    private final boolean DEBUG;
    private final PairHMM pairHMM;
    public Double haplotypeLikehoodMatrix[][];

    static {
        LOG_ONE_HALF = -Math.log10(2.0);
    }

    public LikelihoodCalculationEngine( final byte constantGCP, final boolean debug, final boolean noBanded ) {
        pairHMM = new PairHMM( noBanded );
        this.constantGCP = constantGCP;
        DEBUG = debug;
    }

    public void computeLikelihoods( final ArrayList<Haplotype> haplotypes, final ArrayList<GATKSAMRecord> reads ) {
        final int numHaplotypes = haplotypes.size();
        final double readLikelihoods[][] = new double[reads.size()][numHaplotypes];
        haplotypeLikehoodMatrix = new Double[numHaplotypes][numHaplotypes];

        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
                haplotypeLikehoodMatrix[iii][jjj] = 0.0;
            }
        }

        int maxHaplotypeLength = 0;
        for( final Haplotype h : haplotypes ) {
            int length = h.getBases().length;
            if(length > maxHaplotypeLength) { maxHaplotypeLength = length; }
        }

        for( int iii = 0; iii < reads.size(); iii++ ) {
            final GATKSAMRecord read = reads.get(iii);
            if( DEBUG ) { System.out.println(read.toString()); }
            final byte[] overallGCP = new byte[read.getReadLength()];
            Arrays.fill( overallGCP, constantGCP ); // Is there a way to derive empirical estimates for this from the data?
            for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {                
                readLikelihoods[iii][jjj] = pairHMM.computeReadLikelihoodGivenHaplotype(haplotypes.get(jjj).getBases(), read.getReadBases(),
                        read.getBaseQualities(), read.getBaseInsertionQualities(), read.getBaseDeletionQualities(), overallGCP);
                if( DEBUG ) { System.out.println(">> " + jjj + ":\t" + readLikelihoods[iii][jjj]); }
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

                    /*
                    final GATKSAMRecord read = reads.get(kkk);
                    final int mappingLength = read.getAlignmentEnd() - read.getAlignmentStart() + 1;
                    final double mappingProb = 1.0 - Math.max(0.0, (100.0 - ((double)mappingLength)) / 100.0); //BUGBUG: 101!, needs to pull from the empirical read length distribution per read group
                    haplotypeLikehoodMatrix[iii][jjj] += (mappingProb*mappingProb) * ( MathUtils.softMax(readLikelihoods[kkk][iii], readLikelihoods[kkk][jjj]) + LOG_ONE_HALF ); // BUGBUG: needs to be a logged probability
                    */
                    haplotypeLikehoodMatrix[iii][jjj] += MathUtils.approximateLog10SumLog10(readLikelihoods[kkk][iii], readLikelihoods[kkk][jjj]) + LOG_ONE_HALF;
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

    public ArrayList<Haplotype> chooseBestHaplotypes( final ArrayList<Haplotype> haplotypes ) {

        // For now we choose the top two haplotypes by finding the max value of the pairwise matrix
        // in the future we could use AIC or some other criterion to select more haplotypes to best explain the read data

        final int numHaplotypes = haplotypes.size();
        final ArrayList<Haplotype> bestHaplotypesList = new ArrayList<Haplotype>();
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

        bestHaplotypesList.add(haplotypes.get(hap1));
        if( hap1 != hap2 ) { bestHaplotypesList.add(haplotypes.get(hap2)); }

        return bestHaplotypesList;
    }
}