/*
 * Copyright (c) 2010 The Broad Institute
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

import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;

import java.util.*;

public abstract class PoolGenotypeLikelihoods {
    protected final int numChromosomes;

    protected static final boolean VERBOSE = false;

    //
    // The fundamental data arrays associated with a Genotype Likelhoods object
    //
    protected double[] log10Likelihoods;

    protected final int nSamplesPerPool;
    protected final HashMap<String, ErrorModel> perLaneErrorModels;
    protected final int likelihoodDim;
    protected final boolean ignoreLaneInformation;

    protected final int nAlleles;
    protected final List<Allele> alleles;

    private static final int MAX_NUM_ALLELES_TO_CACHE = 20;
    private static final int MAX_NUM_SAMPLES_PER_POOL = 1000;
    // constructor with given logPL elements
    public PoolGenotypeLikelihoods(final List<Allele> alleles, final double[] logLikelihoods, final int ploidy,
                                   final HashMap<String, ErrorModel> perLaneErrorModels, final boolean ignoreLaneInformation) {
        this.alleles = alleles;
        this.nAlleles = alleles.size();
        numChromosomes = ploidy;
        nSamplesPerPool = numChromosomes/2;
        this.perLaneErrorModels = perLaneErrorModels;
        this.ignoreLaneInformation = ignoreLaneInformation;
        // check sizes
        if (nAlleles > MAX_NUM_ALLELES_TO_CACHE)
            throw new UserException("No support for this number of alleles");

        if (nSamplesPerPool > MAX_NUM_SAMPLES_PER_POOL)
            throw new UserException("No support for such large number of samples per pool");

        likelihoodDim = GenotypeLikelihoods.calculateNumLikelihoods(nAlleles, numChromosomes);

        if (logLikelihoods == null){
            log10Likelihoods = new double[likelihoodDim]; 
            Arrays.fill(log10Likelihoods,Double.NEGATIVE_INFINITY);
        } else {
            if (logLikelihoods.length != likelihoodDim)
                throw new ReviewedStingException("BUG: inconsistent parameters when creating PoolGenotypeLikelihoods object");

            log10Likelihoods = logLikelihoods; //.clone(); // is clone needed?
        }           
    }


    /**
     * Crucial inner class that handles addressing elements of pool likelihoods. We store likelihoods as a map
     * of form int[] -> double (to be more precise, IntArrayWrapper -> Double).
     * For a given ploidy (chromosome count) and number of alleles, we need a form to iterate deterministically
     * across all possible allele conformations.
     * Problem equivalent to listing in determistic order all possible ways in which N integers will sum to P,
     * where N is number of alleles and P is number of chromosomes.
     * There's an option to list all integers so that sum will be UP to P.
     * For example, with P=2,N=2, restrictSumTo = 2 iterator will produce
     * [2 0 ] [1 1] [ 0 2]
     *
     *
     */
    protected static class SumIterator {
        private int[] currentState;
        private final int[] finalState;
        private final int restrictSumTo;
        private final int dim;
        private boolean hasNext;
        private int linearIndex;
        private int currentSum;

        /**
         * Default constructor. Typical use case: restrictSumTo = -1 if there's no sum restriction, or will generate int[]
         * vectors so that all add to this value.
         *
         * @param finalState                    End state - typically we should set value to (P,P,P,...)
         * @param restrictSumTo                 See above
         */
        public SumIterator(final int[] finalState,final int restrictSumTo) {
            this.finalState = finalState;
            this.dim = finalState.length;
            this.restrictSumTo = restrictSumTo;
            currentState = new int[dim];
            reset();

        }

        /**
         * Shortcut constructor for common use case: iterator will produce 
         * all vectors of length numAlleles whose sum = numChromosomes
         * @param numAlleles              Number of alleles
         * @param numChromosomes          Ploidy
         */
        public SumIterator(final int numAlleles, final int numChromosomes) {
            this(getInitialStateVector(numAlleles,numChromosomes), numChromosomes);            
        }


        private static int[] getInitialStateVector(final int nAlleles, final int numChromosomes) {
            int[] initialState = new int[nAlleles];
            Arrays.fill(initialState,numChromosomes);
            return initialState;
        }
        
        public void next() {
            int initialDim = (restrictSumTo > 0)?1:0;
            hasNext = next(finalState, initialDim);
            if (hasNext)
                linearIndex++;
        }

        private boolean next(final int[] finalState, int initialDim) {
            boolean hasNextState = false;
            for (int currentDim=initialDim; currentDim < finalState.length; currentDim++) {
                final int x = currentState[currentDim]+1;

                if (x > finalState[currentDim] || (currentSum >= restrictSumTo && initialDim > 0)) {
                    // update vector sum, and reset position
                    currentSum -= currentState[currentDim];
                    currentState[currentDim] = 0;
                    if (currentDim >= dim-1) {
                        hasNextState = false;
                        break;
                    }
                }
                else {
                    currentState[currentDim] = x;
                    hasNextState = true;
                    currentSum++;
                    break;
                }
            }
            if (initialDim > 0) {
                currentState[0] = restrictSumTo - currentSum;
            }
            return hasNextState;
        }

        public void reset() {
            Arrays.fill(currentState, 0);
            if (restrictSumTo > 0)
                currentState[0] = restrictSumTo;
            hasNext = true;
            linearIndex = 0;
            currentSum = 0;
        }
        public int[] getCurrentVector() {
            return currentState;
        }
      /*  public int getCurrentSum() {
            return currentSum;
        }
        */
        public int getLinearIndex() {
            return linearIndex;
        }

        public boolean hasNext() {
            return hasNext;
        }
    }

    public List<Allele> getAlleles() { return alleles;}

    /**
     * Returns an array of log10 likelihoods for each genotype conformation, with ordering determined by SumIterator class.
     *
     * @return likelihoods array
     */
    public double[] getLikelihoods() {
        return log10Likelihoods;
    }





    /**
     * Set particular element of logPL vector
     * @param idx          index of allele count conformation to modify
     * @param pl                    Likelihood to associate with map
     */
    public void setLogPLs(final int idx, final double pl) {
            log10Likelihoods[idx] = pl;
    }

    public void renormalize() {
        log10Likelihoods = MathUtils.normalizeFromLog10(log10Likelihoods,false,true);
    }
    /** Compute most likely AC conformation based on currently stored PL's - just loop through log PL map and output max value
     *
     * @return vector with most likely allele count, ordered according to this object's alleles
     */
    public Pair<int[],Double> getMostLikelyACCount() {

        int[] mlInd = null;
        double maxVal = Double.NEGATIVE_INFINITY;

        final SumIterator iterator = new SumIterator(alleles.size(),numChromosomes);

        int idx = 0;
        while (iterator.hasNext()) {
            double pl = log10Likelihoods[idx++];
            if (pl > maxVal) {
                maxVal = pl;
                mlInd = iterator.getCurrentVector().clone();
                
            }
            iterator.next();
        }
        if (VERBOSE) {
            System.out.print("MLAC: ");
            outputVectorAsString(mlInd);
            System.out.println();
        }
        return new Pair<int[], Double>(mlInd,maxVal);
    }

    /**
     * Given set of alleles with corresponding vector of likelihoods, subset to a new set of alleles
     *
     * @param oldLikelihoods        Vector of PL's corresponding to original alleles
     * @param numChromosomes        Ploidy (number of chromosomes describing PL's)
     * @param originalAlleles       List of original alleles
     * @param allelesToSubset       Alleles to subset
     * @return                      Vector of new PL's, ordered accorrding to SumIterator's ordering
     */
    public static double[] subsetToAlleles(final double[] oldLikelihoods, final int numChromosomes,
                                                   final List<Allele> originalAlleles, final List<Allele> allelesToSubset) {

        int newPLSize = PoolGenotypeLikelihoods.getNumLikelihoodElements(allelesToSubset.size(), numChromosomes);
        double[] newPLs = new double[newPLSize];


        int idx = 0;
        // First fill boolean array stating whether each original allele is present in new mapping
        final boolean [] allelePresent = new boolean[originalAlleles.size()];
        for ( Allele allele : originalAlleles )
            allelePresent[idx++] = allelesToSubset.contains(allele);


        // compute mapping from old idx to new idx
        // This might be needed in case new allele set is not ordered in the same way as old set
        // Example. Original alleles: {T*,C,G,A}. New alleles: {G,C}. Permutation key = [2,1]

        int[] permutationKey = new int[allelesToSubset.size()];
        for (int k=0; k < allelesToSubset.size(); k++)
            // for each allele to subset, find corresponding index in original allele list
            permutationKey[k] = originalAlleles.indexOf(allelesToSubset.get(k));


        if (VERBOSE) {
            System.out.print("permutationKey:");
            outputVectorAsString(permutationKey);    
        }

        final SumIterator iterator = new SumIterator(originalAlleles.size(),numChromosomes);

        while (iterator.hasNext()) {
            // for each entry in logPL table, associated originally with allele count stored in vec[],
            // see if this allele count conformation will be present in new logPL table.
            // For entry to be present, elements in dimensions not present in requested allele list have to have count = 0
            int[] pVec = iterator.getCurrentVector();
            double pl = oldLikelihoods[iterator.getLinearIndex()];
            
            boolean keyPresent = true;
            for (int k=0; k < allelePresent.length; k++)
                if ( pVec[k]>0 && !allelePresent[k] )
                    keyPresent = false;

            if (keyPresent) {// skip to next entry in logPLs if this conformation is not present in subset

                final int[] newCount = new int[allelesToSubset.size()];
    
                // map from old allele mapping count to new allele mapping
                // In pseudo-Matlab notation: newCount = vec[permutationKey] for permutationKey vector
                for (idx = 0; idx < newCount.length; idx++)
                    newCount[idx] =  pVec[permutationKey[idx]];
    
                // get corresponding index from new count
                int outputIdx = PoolGenotypeLikelihoods.getLinearIndex(newCount, allelesToSubset.size(), numChromosomes);
                newPLs[outputIdx] = pl;
                if (VERBOSE) {
                    System.out.print("Old Key:");
                    outputVectorAsString(pVec);
                    System.out.print("New Key:");
                    outputVectorAsString(newCount);
                }
            }
            iterator.next();
        }

        return  newPLs;
    }

    public static int getLinearIndex(int[] vectorIdx, int numAlleles, int numChromosomes) {

        // brain-dead implementation.
        // BIG to-do, ideally, should compute closed form formula for this
        final SumIterator iterator = new SumIterator(numAlleles,numChromosomes);

 
        while (iterator.hasNext()) {
            int[] vec = iterator.getCurrentVector();
            if (Arrays.equals(vec,vectorIdx))
                return iterator.getLinearIndex();
            iterator.next();
            
        }


        return -1;
        
    }

    /**
     * Given a scalar index, what's the alelle count conformation corresponding to it?
     * @param nAlleles                    Number of alleles
     * @param numChromosomes              Ploidy
     * @param PLindex                     Index to query
     * @return                            Allele count conformation, according to iteration order from SumIterator
     */
    public static int[] getAlleleCountFromPLIndex(final int nAlleles, final int numChromosomes, final int PLindex) {

        // todo - another brain-dead inefficient implementation, can do much better by computing in closed form
        final SumIterator iterator = new SumIterator(nAlleles,numChromosomes);
        while (iterator.hasNext()) {
            final int[] plVec = iterator.getCurrentVector();
            if (iterator.getLinearIndex() == PLindex)
                return plVec;

            iterator.next();
        }

        return null;

    }

    /*
    * a cache of the PL ivector sizes as a function of # of alleles and pool sizes
    */
    
    public static int getNumLikelihoodElements(int numAlleles, int ploidy) {
        return GenotypeLikelihoodVectorSizes[numAlleles][ploidy];
    }

    private final static int[][] GenotypeLikelihoodVectorSizes = fillGLVectorSizeCache(MAX_NUM_ALLELES_TO_CACHE, 2*MAX_NUM_SAMPLES_PER_POOL);

    private static int[][] fillGLVectorSizeCache(int maxAlleles, int maxPloidy) {
        
        int[][] cache = new int[maxAlleles][maxPloidy];
        for (int numAlleles=1; numAlleles < maxAlleles; numAlleles++) {
            for (int ploidy=0; ploidy < maxPloidy; ploidy++) {

                if (numAlleles == 1)
                    cache[numAlleles][ploidy] = 1;
                else if (ploidy == 1)
                    cache[numAlleles][ploidy] = numAlleles;
                else {
                    int acc =0;
                    for (int k=0; k <= ploidy; k++ )
                        acc += cache[numAlleles-1][ploidy-k];

                    cache[numAlleles][ploidy] = acc;
                }
            }
        }
        return cache;
    }

    /**
     * Return a string representation of this object in a moderately usable form
     *
     * @return string representation
     */
    public String toString() {
        StringBuilder s = new StringBuilder(1000);

        s.append("Alleles:");
        for (Allele a: this.alleles){
            s.append(a.getDisplayString());
            s.append(",");
        }
        s.append("\nGLs:\n");
        SumIterator iterator = new SumIterator(nAlleles,numChromosomes);
        while (iterator.hasNext()) {
            s.append("Count [");
            StringBuilder b = new StringBuilder(iterator.getCurrentVector().length*2);
            for (int it:iterator.getCurrentVector()) {
                b.append(it);
                b.append(",");
            }
            s.append(b.toString());
            s.append(String.format("] GL=%4.3f\n",this.getLikelihoods()[iterator.getLinearIndex()]) );
            iterator.next();
        }
        return s.toString();
    }


    // small helper routines to dump primitive array to screen
    public static void outputVectorAsString(int[] array) {
        for (int d:array) 
            System.out.format("%d ",d);
        System.out.println();
    }
}
