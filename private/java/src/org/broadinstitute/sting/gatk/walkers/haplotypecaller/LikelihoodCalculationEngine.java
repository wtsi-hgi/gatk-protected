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

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.PairHMM;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

public class LikelihoodCalculationEngine {

    private final static double LOG_ONE_HALF = -Math.log10(2.0);
    private final byte constantGCP;
    private final boolean DEBUG;
    private final PairHMM pairHMM;

    public LikelihoodCalculationEngine( final byte constantGCP, final boolean debug, final boolean noBanded ) {
        pairHMM = new PairHMM( noBanded );
        this.constantGCP = constantGCP;
        DEBUG = debug;
    }

    public void computeReadLikelihoods( final ArrayList<Haplotype> haplotypes, final HashMap<String, ArrayList<GATKSAMRecord>> perSampleReadList ) {
        final int numHaplotypes = haplotypes.size();

        int X_METRIC_LENGTH = 0;
        for( final String sample : perSampleReadList.keySet() ) {
            for( final GATKSAMRecord read : perSampleReadList.get(sample) ) {
                final int readLength = read.getReadLength();
                if( readLength > X_METRIC_LENGTH ) { X_METRIC_LENGTH = readLength; }
            }
        }
        int Y_METRIC_LENGTH = 0;
        for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
            final int haplotypeLength = haplotypes.get(jjj).getBases().length;
            if( haplotypeLength > Y_METRIC_LENGTH ) { Y_METRIC_LENGTH = haplotypeLength; }
        }
        X_METRIC_LENGTH++;
        Y_METRIC_LENGTH++;

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
        
        // for each sample's reads
        for( final String sample : perSampleReadList.keySet() ) {
            if( DEBUG ) { System.out.println("Evaluating sample " + sample + " with " + perSampleReadList.get( sample ).size() + " passing reads"); }
            // evaluate the likelihood of the reads given those haplotypes
            computeReadLikelihoods( haplotypes, perSampleReadList.get(sample), sample, matchMetricArray, XMetricArray, YMetricArray );
        }
    }

    private void computeReadLikelihoods( final ArrayList<Haplotype> haplotypes, final ArrayList<GATKSAMRecord> reads, final String sample,
                                         final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray ) {

        final int numHaplotypes = haplotypes.size();
        final int numReads = reads.size();
        final double[][] readLikelihoods = new double[numHaplotypes][numReads];
        
        for( int iii = 0; iii < numReads; iii++ ) {
            final GATKSAMRecord read = reads.get(iii);
            final byte[] overallGCP = new byte[read.getReadLength()];
            Arrays.fill( overallGCP, constantGCP ); // Is there a way to derive empirical estimates for this from the data?
            Haplotype previousHaplotypeSeen = null;
            for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {            
                final Haplotype haplotype = haplotypes.get(jjj);
                final int haplotypeStart = ( previousHaplotypeSeen == null ? 0 : computeFirstDifferingPosition(haplotype.getBases(), previousHaplotypeSeen.getBases()) );
                previousHaplotypeSeen = haplotype;
                readLikelihoods[jjj][iii] = pairHMM.computeReadLikelihoodGivenHaplotype(haplotype.getBases(), read.getReadBases(),
                        read.getBaseQualities(), read.getBaseInsertionQualities(), read.getBaseDeletionQualities(), overallGCP,
                        haplotypeStart, matchMetricArray, XMetricArray, YMetricArray);
            }
        }
        for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
            haplotypes.get(jjj).addReadLikelihoods( sample, readLikelihoods[jjj] );
        }
    }

    private static int computeFirstDifferingPosition( final byte[] b1, final byte[] b2 ) {
        for( int iii = 0; iii < b1.length && iii < b2.length; iii++ ){
            if( b1[iii] != b2[iii] ) {
                return iii;
            }
        }
        return b1.length;
    }

    @Requires({"haplotypes.size() > 0"})
    @Ensures({"result.length == result[0].length", "result.length == haplotypes.size()"})
    public static double[][] computeDiploidHaplotypeLikelihoods( final ArrayList<Haplotype> haplotypes, final String sample ) {
        // set up the default 1-to-1 haplotype mapping object
        final ArrayList<ArrayList<Integer>> haplotypeMapping = new ArrayList<ArrayList<Integer>>();
        final int numHaplotypes = haplotypes.size();
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            final ArrayList<Integer> list = new ArrayList<Integer>();
            list.add(iii);
            haplotypeMapping.add(list);
        }
        return computeDiploidHaplotypeLikelihoods( haplotypes, sample, haplotypeMapping );
    }
    
    @Requires({"haplotypes.size() > 0","haplotypeMapping.size() > 0","haplotypeMapping.size() <= haplotypes.size()"})
    @Ensures({"result.length == result[0].length", "result.length == haplotypeMapping.size()"})
    public static double[][] computeDiploidHaplotypeLikelihoods( final ArrayList<Haplotype> haplotypes, final String sample, final ArrayList<ArrayList<Integer>> haplotypeMapping ) {

        final int numHaplotypes = haplotypeMapping.size();
        final int numReads = haplotypes.get(0).getReadLikelihoods(sample).length; // BUGBUG: assume all haplotypes saw the same reads
        final double[][] haplotypeLikelihoodMatrix = new double[numHaplotypes][numHaplotypes];
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            Arrays.fill(haplotypeLikelihoodMatrix[iii], 0.0);
        }

        // compute the diploid haplotype likelihoods
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( final int iii_mapped : haplotypeMapping.get(iii) ) {
                final double[] readLikelihoods_iii = haplotypes.get(iii_mapped).getReadLikelihoods(sample);
                for( int jjj = 0; jjj <= iii; jjj++ ) {
                    for( final int jjj_mapped : haplotypeMapping.get(jjj) ) {
                        final double[] readLikelihoods_jjj = haplotypes.get(jjj_mapped).getReadLikelihoods(sample);
                        for( int kkk = 0; kkk < numReads; kkk++ ) {
                            // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                            // First term is approximated by Jacobian log with table lookup.
                            haplotypeLikelihoodMatrix[iii][jjj] += MathUtils.approximateLog10SumLog10(readLikelihoods_iii[kkk], readLikelihoods_jjj[kkk]) + LOG_ONE_HALF;
                        }
                    }
                }
            }
        }

        // normalize the diploid likelihoods matrix
        return normalizeDiploidLikelihoodMatrixFromLog10( haplotypeLikelihoodMatrix );        
    }

    @Requires({"likelihoodMatrix.length == likelihoodMatrix[0].length"})
    @Ensures({"result.length == result[0].length", "result.length == likelihoodMatrix.length"})
    protected static double[][] normalizeDiploidLikelihoodMatrixFromLog10( final double[][] likelihoodMatrix ) {
        final int numHaplotypes = likelihoodMatrix.length;
        double[] genotypeLikelihoods = new double[numHaplotypes*(numHaplotypes+1)/2];
        int index = 0;
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj <= iii; jjj++ ){
                genotypeLikelihoods[index++] = likelihoodMatrix[iii][jjj];
            }
        }
        genotypeLikelihoods = MathUtils.normalizeFromLog10(genotypeLikelihoods, false, true);
        index = 0;
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj <= iii; jjj++ ){
                likelihoodMatrix[iii][jjj] = genotypeLikelihoods[index++];
            }
        }
        return likelihoodMatrix;
    }

    @Requires({"haplotypes.size() > 0"})
    @Ensures({"result.size() <= haplotypes.size()"})
    public ArrayList<Haplotype> selectBestHaplotypes( final ArrayList<Haplotype> haplotypes ) {

        // For now we choose the top two haplotypes by finding the max value per sample of the diploid haplotype likelihood matrix
        // in the future we could use AIC or some other criterion to select more haplotypes to best explain the read data

        final int numHaplotypes = haplotypes.size();
        final Set<String> stringKeySet = haplotypes.get(0).getSampleKeySet(); // BUGBUG: assume all haplotypes saw the same samples
        final ArrayList<Integer> bestHaplotypesIndexList = new ArrayList<Integer>();
        bestHaplotypesIndexList.add(0); // always start with the reference haplotype
        
        for( final String sample : stringKeySet ) {
            final double[][] haplotypeLikelihoodMatrix = computeDiploidHaplotypeLikelihoods( haplotypes, sample );
                    
            double maxElement = Double.NEGATIVE_INFINITY;
            int hap1 = -1;
            int hap2 = -1;
            for( int iii = 0; iii < numHaplotypes; iii++ ) {
                for( int jjj = 0; jjj <= iii; jjj++ ) {
                    if( haplotypeLikelihoodMatrix[iii][jjj] > maxElement ) {
                        maxElement = haplotypeLikelihoodMatrix[iii][jjj];
                        hap1 = iii;
                        hap2 = jjj;
                    }
                }
            }
    
            if( !bestHaplotypesIndexList.contains(hap1) ) { bestHaplotypesIndexList.add(hap1); }
            if( !bestHaplotypesIndexList.contains(hap2) ) { bestHaplotypesIndexList.add(hap2); }
        }

        final ArrayList<Haplotype> bestHaplotypes = new ArrayList<Haplotype>();
        for( final int hIndex : bestHaplotypesIndexList ) {
            bestHaplotypes.add( haplotypes.get(hIndex) );
        }
        return bestHaplotypes;
    }
}