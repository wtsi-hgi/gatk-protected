package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.walkers.poolcaller.ErrorModel;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 4/1/12
 * Time: 7:28 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class PoolGenotypeLikelihoods {
    protected final int numChromosomes;

    protected boolean VERBOSE = false;

    //
    // The fundamental data arrays associated with a Genotype Likelhoods object
    //
    protected double[] log10Likelihoods = null;
    protected double[] log10Posteriors = null;
    protected final double[] genotypeZeros;
//    private int[] initialState;

    protected PoolGenotypePriors priors = null;

    protected final int nSamplesPerPool;
    protected final HashMap<String, ErrorModel> perLaneErrorModels;

    protected final boolean ignoreLaneInformation;

    protected Map<IntArrayWrapper,Double> logPLs;
    protected ArrayList<IntArrayWrapper> alleleConformationList;
    
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

    protected static class SumIterator {
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
     * Returns an array of log10 likelihoods for each genotype
     * @return likelihoods array
     */
    public double[] getLikelihoods() {
        SumIterator iterator = new SumIterator(getInitialStateVector(nAlleles,numChromosomes),numChromosomes);
        int idx = 0;
        int n = GenotypeLikelihoods.calculateNumLikelihoods(nAlleles, numChromosomes);
        double[] logPLVector = new double[n];
        while (iterator.hasNext()) {
            logPLVector[idx++] = logPLs.get(new IntArrayWrapper(iterator.getCurrentVector()));
            iterator.next();

        }
        if (idx != n)
            throw new ReviewedStingException("BUG: inconsistent size of logPLs");
        return logPLVector;
    }

    protected static int[] getInitialStateVector(int nAlleles, int numChromosomes) {
        int[] initialState = new int[nAlleles];
        Arrays.fill(initialState,numChromosomes);
        return initialState;
    }

    private void initializePLs() {
        int vecDim = GenotypeLikelihoods.calculateNumLikelihoods(nAlleles, numChromosomes);
        double[] ss = new double[vecDim];
        Arrays.fill(ss,Double.NEGATIVE_INFINITY);
        //alleleConformationList = new ArrayList<IntArrayWrapper>(vecDim);
        logPLs = initializePLsFromVector(nAlleles, numChromosomes, ss/*, alleleConformationList*/);
    }

    public static Map<IntArrayWrapper,Double> initializePLsFromVector(final int nAlleles, final int numChromosomes, 
                                                                      final double[] logVec/*, ArrayList<IntArrayWrapper> alleleConformationList*/) {
        SumIterator iterator = new SumIterator(getInitialStateVector(nAlleles,numChromosomes),numChromosomes);
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

    public static List<Integer> asIntegerList(int[] vec) {
        Integer[] ivec = new Integer[vec.length];
        int idx = 0;
        for (int k:vec)
            ivec[idx++] = k;
        return Arrays.asList(ivec);

    }
    
    public static int[] asIntArray(List<Integer> vec) {
        int[] ivec = new int[vec.size()];
        
        int idx = 0;
        for (Integer k : vec)
            ivec[idx++] = k;
        
        return ivec;
    }
    
    public static int[] getAlleleCountFromPLIndex(int nAlleles, int numChromosomes, int PLindex) {
        SumIterator iterator = new SumIterator(getInitialStateVector(nAlleles,numChromosomes),numChromosomes);
        while (iterator.hasNext()) {
            int[] plVec = iterator.getCurrentVector();
            if (iterator.getLinearIndex() == PLindex)
                return plVec;

            iterator.next();
        }

        return null;

    }
    public static double[] getPLVectorFromPLObject(int nAlleles, int numChromosomes, Map<IntArrayWrapper,Double> logPLs) {

        SumIterator iterator = new SumIterator(getInitialStateVector(nAlleles,numChromosomes),numChromosomes);
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

    public void setLogPLs(int[] alleleCounts, double pl) {
        if (alleleCounts.length != nAlleles || !logPLs.containsKey(new IntArrayWrapper(alleleCounts)))
            throw new ReviewedStingException("BUG: trying to call setLogPLs with wrong key");
        
        logPLs.put(new IntArrayWrapper(alleleCounts), pl);
        
    }

    public double getLogPLofAC(final int[] ac) {
        return logPLs.get(new IntArrayWrapper(ac));
    }

    public int[] getMostLikelyACCount() {

        int[] mlInd = new int[nAlleles];
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
    
    public static double[] subsetToAlleles(final double[] oldLikelihoods, final int numChromosomes,
                                                   final List<Allele> originalAlleles, final List<Allele> allelesToSubset) {
        Map<IntArrayWrapper,Double> oldPLs = initializePLsFromVector(originalAlleles.size(),numChromosomes, oldLikelihoods);

        Map<IntArrayWrapper,Double> newPLs = new HashMap<IntArrayWrapper,Double>();


        // compute mapping from old idx to new idx
        int idx = 0;
        int[] permutationKey = new int[allelesToSubset.size()];
        final boolean [] allelePresent = new boolean[originalAlleles.size()];
        for ( Allele allele : originalAlleles )
            allelePresent[idx++] = allelesToSubset.contains(allele);

        // 
        for (int k=0; k < permutationKey.length; k++) {
            // for each allele to subset, find corresponding index in original allele list
            for (int a=0; a < originalAlleles.size(); a++) {
                if (allelesToSubset.get(k) == originalAlleles.get(a)) {
                    permutationKey[k] = a;
                    break;
                }
            }
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

            for (idx = 0; idx < newCount.data.length; idx++)
                newCount.data[idx] =  vec.data[permutationKey[idx]];

            // compute new key: all dimensions not present in
            newPLs.put(newCount, elt.getValue());
        }

        return getPLVectorFromPLObject(allelesToSubset.size(), numChromosomes, newPLs);
    }

}
