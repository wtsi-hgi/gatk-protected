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

import org.broadinstitute.sting.gatk.walkers.poolcaller.ErrorModel;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;

import java.util.*;

public abstract class PoolGenotypeLikelihoods {
    protected final int numChromosomes;

    protected static final boolean VERBOSE = false;

    //
    // The fundamental data arrays associated with a Genotype Likelhoods object
    //
    protected double[] log10Likelihoods = null;
    protected double[] log10Posteriors = null;
    protected final double[] genotypeZeros;

    protected PoolGenotypePriors priors = null;

    protected final int nSamplesPerPool;
    protected final HashMap<String, ErrorModel> perLaneErrorModels;

    protected final boolean ignoreLaneInformation;

    protected Map<IntArrayWrapper,Double> logPLs;
//    protected ArrayList<IntArrayWrapper> alleleConformationList;
    
    protected final int nAlleles;
    protected final List<Allele> alleles;

    public PoolGenotypeLikelihoods(final List<Allele> alleles, final PoolGenotypePriors priors,
                                   final HashMap<String, ErrorModel> perLaneErrorModels, final boolean ignoreLaneInformation) {
        this.priors = priors;
        this.alleles = alleles;
        this.nAlleles = alleles.size();
        nSamplesPerPool = priors.getNSamplesPerPool();
        numChromosomes = 2*nSamplesPerPool;
        genotypeZeros = new double[priors.getPriors().length];
        this.perLaneErrorModels = perLaneErrorModels;
        this.ignoreLaneInformation = ignoreLaneInformation;
        initializePLs();

    }

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
        if (logLikelihoods.length != GenotypeLikelihoods.calculateNumLikelihoods(nAlleles,ploidy))
            throw new ReviewedStingException("BUG: inconsistent parameters when creating PoolGenotypeLikelihoods object");
        logPLs = initializePLsFromVector(nAlleles, numChromosomes, logLikelihoods);
        genotypeZeros = new double[logLikelihoods.length];
    }

    /**
     * Wrapper class because we can't directly implement a hash of int[] -> double:
     * In java, an array of primitive types is not equal to another one when an element-wise comparison is made,
     * so this just encapsulates the array and adds the equals method to do element-wise comparison.
     *
     */
    protected final static class IntArrayWrapper /*implements Comparable<Object> */
    {
        private final int[] data;

        public IntArrayWrapper(int[] data)
        {
            if (data == null)
            {
                throw new NullPointerException();
            }
            this.data = data.clone();
        }

        @Override
        public boolean equals(Object other)
        {
            if (!(other instanceof IntArrayWrapper))
            {
                return false;
            }
            return Arrays.equals(data, ((IntArrayWrapper)other).data);
        }

        @Override
        public int hashCode()
        {
            return Arrays.hashCode(data);
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
            if (restrictSumTo>0) {
                next();
                linearIndex = 0;
            }
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
            if (restrictSumTo > 0) {
                do {
                    hasNext = next(finalState, 0);
                    if (!hasNext)
                        break;
                }
                while(getCurrentSum() != restrictSumTo);

            }
            else
                hasNext = next(finalState, 0);

            if (hasNext)
                linearIndex++;
        }
        private boolean next(final int[] finalState, final int initialDim) {
            final int x = currentState[initialDim]+1;

            if (x > finalState[initialDim]) {
                // recurse into subvector
                currentState[initialDim] = 0;
                if (initialDim < dim-1) {
                    return next(finalState,initialDim+1);
                }
                else
                    return false;
            }
            else
                currentState[initialDim] = x;
            
            return true;

        }

        public void reset() {
            Arrays.fill(currentState, 0);
            hasNext = true;
            linearIndex = 0;
        }
        public int[] getCurrentVector() {
            return currentState;
        }
        public int getCurrentSum() {
            return (int) MathUtils.sum(currentState);
        }

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
        final SumIterator iterator = new SumIterator(nAlleles,numChromosomes);
        int idx = 0;
        final int n = GenotypeLikelihoods.calculateNumLikelihoods(nAlleles, numChromosomes);
        double[] logPLVector = new double[n];
        while (iterator.hasNext()) {
            logPLVector[idx++] = logPLs.get(new IntArrayWrapper(iterator.getCurrentVector()));
            iterator.next();

        }
        
        // sanity check. Did we traverse the whole logPL map?
        if (idx != n)
            throw new ReviewedStingException("BUG: inconsistent size of logPLs");
        return logPLVector;
    }


    /**
     * Initialize hash map with correct keys and set all likelihoods to -infinity 
     */
    private void initializePLs() {
        final int vecDim = GenotypeLikelihoods.calculateNumLikelihoods(nAlleles, numChromosomes);
        double[] ss = new double[vecDim];
        Arrays.fill(ss,Double.NEGATIVE_INFINITY);
        logPLs = initializePLsFromVector(nAlleles, numChromosomes, ss);
    }

    /**
     * Given a particular vector of PL's, construct a hash map of int[]-> double with these PL values.
     * @param nAlleles               Number of alleles
     * @param numChromosomes         Ploidy
     * @param logVec                 Vector of likelihoods to fill
     * @return                       Map filled with logVec
     */
    public static Map<IntArrayWrapper,Double> initializePLsFromVector(final int nAlleles, final int numChromosomes, 
                                                                      final double[] logVec) {
        final SumIterator iterator = new SumIterator(nAlleles,numChromosomes);
        Map<IntArrayWrapper,Double> logPLs = new HashMap<IntArrayWrapper,Double>();

        int k=0;
        while (iterator.hasNext()) {
            IntArrayWrapper key = new IntArrayWrapper( iterator.getCurrentVector());
            //alleleConformationList.add(k,key);
            logPLs.put(key,logVec[k++]);

            iterator.next();

        }
        return logPLs;
    }

    /**
     * Given a scalar index, what's the alelle count conformation corresponding to it?
     * @param nAlleles                    Number of alleles
     * @param numChromosomes              Ploidy
     * @param PLindex                     Index to query
     * @return                            Allele count conformation, according to iteration order from SumIterator
     */
    public static int[] getAlleleCountFromPLIndex(final int nAlleles, final int numChromosomes, final int PLindex) {
        final SumIterator iterator = new SumIterator(nAlleles,numChromosomes);
        while (iterator.hasNext()) {
            final int[] plVec = iterator.getCurrentVector();
            if (iterator.getLinearIndex() == PLindex)
                return plVec;

            iterator.next();
        }

        return null;

    }

    /**
     * Helper routine. Given Map of int[] -> double in canonical ordering, unwrap it to linear vector,
     * with ordering given by SumIterator.
     * @param nAlleles                          Number of Alleles
     * @param numChromosomes                    Ploidy (number of chromosomes)
     * @param logPLs                            Map of IntArrayWrapper to Double containing actual likelihood values
     * @return                                  Likelihood vector
     */
    public static double[] getPLVectorFromPLObject(final int nAlleles, final int numChromosomes, 
                                                   final Map<IntArrayWrapper,Double> logPLs) {

        final SumIterator iterator = new SumIterator(nAlleles,numChromosomes);
        int k=0;
        double[] logVec = new double[GenotypeLikelihoods.calculateNumLikelihoods(nAlleles, numChromosomes)];
        while (iterator.hasNext()) {
            logVec[k++] = logPLs.get(new IntArrayWrapper(iterator.getCurrentVector()));
            iterator.next();
        }

        return logVec;
    }
    /**
     * Returns an array of priors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return log10 prior as a double array
     */
    public double[] getPriors() {
        return priors.getPriors();
    }

    /**
     * Set particular element of logPL map
     * @param alleleCounts          allele count conformation to modify
     * @param pl                    Likelihood to associate with map
     */
    public void setLogPLs(final int[] alleleCounts, final double pl) {
        if (alleleCounts.length != nAlleles || !logPLs.containsKey(new IntArrayWrapper(alleleCounts)))
            throw new ReviewedStingException("BUG: trying to call setLogPLs with wrong key");
        
        logPLs.put(new IntArrayWrapper(alleleCounts), pl);
        
    }

    /**
     * Query likelihood map 
     * @param ac            Allele count conformation to query
     * @return              Likelihood of that conformation
     */
    public double getLogPLofAC(final int[] ac) {
        return logPLs.get(new IntArrayWrapper(ac));
    }

    /** Compute most likely AC conformation based on currently stored PL's - just loop through log PL map and output max value
     *
     * @return vector with most likely allele count, ordered according to this object's alleles
     */
    public int[] getMostLikelyACCount() {

        final int[] mlInd = new int[nAlleles];
        double maxVal = Double.NEGATIVE_INFINITY;
        IntArrayWrapper mlList = new IntArrayWrapper(mlInd);
        
        for (Map.Entry<IntArrayWrapper,Double> elt : logPLs.entrySet()) {
            double val = elt.getValue();
            if (val > maxVal) {
                maxVal = val;
                mlList = elt.getKey();
                
            }
        }
        if (VERBOSE) {
            System.out.print("MLAC: ");
            for (int k:mlInd)
                System.out.format("%d ",k);
            System.out.println();
        }
        return mlList.data;
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
        final Map<IntArrayWrapper,Double> oldPLs = initializePLsFromVector(originalAlleles.size(),numChromosomes, oldLikelihoods);

        Map<IntArrayWrapper,Double> newPLs = new HashMap<IntArrayWrapper,Double>();


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

        for (final Map.Entry<IntArrayWrapper,Double> elt : oldPLs.entrySet()) {
            final IntArrayWrapper vec = elt.getKey();
            // for each entry in logPL table, associated originally with allele count stored in vec[],
            // see if this allele count conformation will be present in new logPL table.
            // For entry to be present, elements in dimensions not present in requested allele list have to have count = 0
            //
            boolean keyPresent = true;
            for (int k=0; k < allelePresent.length; k++)
                if (!allelePresent[k] && vec.data[k]>0)
                    keyPresent = false;

            if (!keyPresent) continue; // skip to next entry in logPLs if this conformation is not present in subset

            final IntArrayWrapper newCount = new IntArrayWrapper(new int[allelesToSubset.size()]);

            // map from old allele mapping count to new allele mapping
            // In pseudo-Matlab notation: newCount = vec[permutationKey] for permutationKey vector
            for (idx = 0; idx < newCount.data.length; idx++)
                newCount.data[idx] =  vec.data[permutationKey[idx]];

            // compute new key: all dimensions not present in
            newPLs.put(newCount, elt.getValue());
            if (VERBOSE) {
                System.out.print("Old Key:");
                outputVectorAsString(vec.data);
                System.out.print("New Key:");
                outputVectorAsString(newCount.data);
            }
        }

        return getPLVectorFromPLObject(allelesToSubset.size(), numChromosomes, newPLs);
    }

    
    // small helper routines to dump primitive array to screen
    public static void outputVectorAsString(int[] array) {
        for (int d:array) 
            System.out.format("%d ",d);
        System.out.println();
    }
    public static void outputVectorAsString(double[] array) {
        for (double d:array)
            System.out.format("%4.3f ",d);
        System.out.println();
    }
}
