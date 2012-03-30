/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.PrintStream;
import java.util.*;

public class PoolAFCalculationModel extends AlleleFrequencyCalculationModel {

    // private final static boolean DEBUG = false;

    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6
    final protected PoolCallerUnifiedArgumentCollection UAC;

    protected PoolAFCalculationModel(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
        if (UAC instanceof PoolCallerUnifiedArgumentCollection)
            this.UAC = (PoolCallerUnifiedArgumentCollection)UAC;
        else
            this.UAC = new PoolCallerUnifiedArgumentCollection(); // dummy copy

    }

    public static class SumIterator {
        private int[] currentState;
        private final int[] finalState;
        private final int restrictSumTo;
        private final int dim;
        private boolean hasNext;
        private int linearIndex;
        
        public SumIterator(int[] finalState,int restrictSumTo) {
            this.finalState = finalState;
            this.dim = finalState.length;
            this.restrictSumTo = restrictSumTo;
            currentState = new int[dim];
            reset();
            if (restrictSumTo>0) {
                next();
                linearIndex = 0;
            }
        }

        public void next() {
            if (restrictSumTo > 0) {
                do {
                    hasNext = next(finalState, dim);
                    if (!hasNext)
                        break;
                }
                while(getCurrentSum() != restrictSumTo);

            }
            else
                hasNext = next(finalState, dim);

            if (hasNext)
                linearIndex++;
        }
        private boolean next(final int[] finalState, final int numValues) {
            int x = currentState[numValues-1]+1;

            if (x > finalState[numValues-1]) {
                // recurse into subvector
                currentState[numValues-1] = 0;
                if (numValues > 1) {
                    return next(finalState,numValues-1);
                }
                else
                    return false;
            }
            else
                currentState[numValues-1] = x;

            return true;
            
        }
        
        public void reset() {
            Arrays.fill(currentState,0); 
            hasNext = true;
            linearIndex = 0;
        }
        public int[] getCurrentVector() {
            return currentState;
        }
        public int getCurrentSum() {
            return (int)MathUtils.sum(currentState);
        }
        
        public int getLinearIndex() {
            return linearIndex;
        }
        
        public boolean hasNext() {
            return hasNext;
        }
    }
    public List<Allele> getLog10PNonRef(final VariantContext vc,
                                        final double[] log10AlleleFrequencyPriors,
                                        final AlleleFrequencyCalculationResult result) {

        GenotypesContext GLs = vc.getGenotypes();
        List<Allele> alleles = vc.getAlleles();

        // don't try to genotype too many alternate alleles
        if ( vc.getAlternateAlleles().size() > MAX_ALTERNATE_ALLELES_TO_GENOTYPE ) {
            logger.warn("this tool is currently set to genotype at most " + MAX_ALTERNATE_ALLELES_TO_GENOTYPE + " alternate alleles in a given context, but the context at " + vc.getChr() + ":" + vc.getStart() + " has " + (vc.getAlternateAlleles().size()) + " alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument");

            alleles = new ArrayList<Allele>(MAX_ALTERNATE_ALLELES_TO_GENOTYPE + 1);
            alleles.add(vc.getReference());
/*            alleles.addAll(chooseMostLikelyAlternateAlleles(vc, MAX_ALTERNATE_ALLELES_TO_GENOTYPE));

            // todo - this is wrong for pool calls!!
            GLs = UnifiedGenotyperEngine.subsetAlleles(vc, alleles, false);
  */      }

        simplePoolBiallelic(GLs, alleles.size() - 1, UAC.nSamplesPerPool, log10AlleleFrequencyPriors, result);

        return alleles;
    }


    private static final ArrayList<double[]> getGLs(GenotypesContext GLs) {
        ArrayList<double[]> genotypeLikelihoods = new ArrayList<double[]>(GLs.size());

        genotypeLikelihoods.add(new double[]{0.0,0.0,0.0}); // dummy
        for ( Genotype sample : GLs.iterateInSampleNameOrder() ) {
            if ( sample.hasLikelihoods() ) {
                double[] gls = sample.getLikelihoods().getAsVector();

                if ( MathUtils.sum(gls) < UnifiedGenotyperEngine.SUM_GL_THRESH_NOCALL )
                    genotypeLikelihoods.add(gls);
            }
        }

        return genotypeLikelihoods;
    }


    public static void simplePoolBiallelic(final GenotypesContext GLs,
                                               final int numAlternateAlleles,
                                               final int nSamplesPerPool,
                                               final double[] log10AlleleFrequencyPriors,
                                               final AlleleFrequencyCalculationResult result) {

        int[] alleleCounts = new int[numAlternateAlleles];
        final ArrayList<double[]> genotypeLikelihoods = getGLs(GLs);
        final int numSamples = genotypeLikelihoods.size()-1;
        final int numChr = 2*(numSamples*nSamplesPerPool);

        double[] combinedLikelihoods = new double[1];
        
        for (int p=1; p<genotypeLikelihoods.size(); p++) {
            combinedLikelihoods = combinePoolNaively(genotypeLikelihoods.get(p), combinedLikelihoods);    
        }

        final double log10Lof0 =combinedLikelihoods[0];
        result.setLog10LikelihoodOfAFzero(log10Lof0);
        result.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);
        
        for (int k=1; k < combinedLikelihoods.length; k++) {
            //Arrays.fill(alleleCounts,0);
            alleleCounts[0] = k;
            result.updateMLEifNeeded(combinedLikelihoods[k],alleleCounts);
            result.updateMAPifNeeded(combinedLikelihoods[k] + log10AlleleFrequencyPriors[k],alleleCounts);
        }
        //result.log10AlleleFrequencyLikelihoods[0] = MathUtils.vectorSum(combinedLikelihoods,log10AlleleFrequencyPriors);
    }


    /* For two pools of size m1 and m2, we can compute the combined likelihood as:
   Pr(D|AC=k) = Sum_{j=0}^k Pr(D|AC1=j) Pr(D|AC2=k-j) * choose(m1,j)*choose(m2,k-j)/choose(m1+m2,k)

    */

    public static double[] combinePoolNaively(double[] x, double[] y) {
     
        int m1 = x.length;
        int m2 = y.length;
        double[] result = new double[m1+m2-1];

        result[0] = x[0]+y[0];
        for (int k=1; k < result.length; k++) {
            double den = MathUtils.log10BinomialCoefficient(m1+m2,k);
            double[] acc = new double[k+1];
            Arrays.fill(acc,Double.NEGATIVE_INFINITY);

            for (int j=0; j <=k; j++) {
                double x1,y1;


                if (k-j>=0 && k-j < y.length)
                    y1 = y[k-j];
                else
                    continue;

                if (j < x.length)
                    x1 = x[j];
                else
                    continue;

                acc[j] = x1+ y1 + MathUtils.log10BinomialCoefficient(m1,j) + MathUtils.log10BinomialCoefficient(m2,k-j) - den;
            }    
            result[k] = MathUtils.log10sumLog10(acc);
        }
        

        return MathUtils.normalizeFromLog10(result,false, true);
    }

}
