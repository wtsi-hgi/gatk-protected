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
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.PrintStream;
import java.util.*;

public class PoolAFCalculationModel extends AlleleFrequencyCalculationModel {

    public static final String MAXIMUM_LIKELIHOOD_AC_KEY = "MLAC";
    public static final String MAXIMUM_LIKELIHOOD_AF_KEY= "MLAF";
    static final int MAX_LENGTH_FOR_POOL_PL_LOGGING = 10; // if PL vectors longer than this # of elements, don't log them
    final protected PoolCallerUnifiedArgumentCollection UAC;

    private final int ploidy;
    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6
    private final static boolean USE_NAIVE_COMBINER = true;
    private final static boolean VERBOSE = false;

    protected PoolAFCalculationModel(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
        if (UAC instanceof PoolCallerUnifiedArgumentCollection) {
            this.UAC = (PoolCallerUnifiedArgumentCollection)UAC;
        }
        else {
            this.UAC = new PoolCallerUnifiedArgumentCollection(); // dummy copy
        }
        ploidy = 2*this.UAC.nSamplesPerPool;


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
            alleles.addAll(chooseMostLikelyAlternateAlleles(vc, MAX_ALTERNATE_ALLELES_TO_GENOTYPE, ploidy));


            GLs = subsetAlleles(vc, alleles, false, ploidy);
        }

        combineSinglePools(GLs, alleles.size() - 1, UAC.nSamplesPerPool, log10AlleleFrequencyPriors, result);

        return alleles;
    }


    static class CombinedPoolLikelihoods {
        private List<Pair<ExactACset,Double>> alleleCountSetList;
        private double maxLikelihood;

        
        //final ExactACset set = new ExactACset(1, new ExactACcounts(outputIterator.getCurrentAltVector()));
        public CombinedPoolLikelihoods(int numAltAlleles) {
            alleleCountSetList = new ArrayList<Pair<ExactACset,Double>>();
            maxLikelihood = Double.NEGATIVE_INFINITY;
        }

        public void add(double likelihood, ExactACset set, int idx) {
            if (idx >= alleleCountSetList.size())
                alleleCountSetList.add(idx,new Pair<ExactACset, Double>(set,likelihood));
            else
                alleleCountSetList.set(idx,new Pair<ExactACset, Double>(set,likelihood));

            if (likelihood > maxLikelihood )
                maxLikelihood = likelihood;

        }

        public double[] getLikelihoodsAsVector(boolean normalize) {
            int k=0;
            double fac = (normalize?maxLikelihood:0.0);
            double[] likelihoods = new double[alleleCountSetList.size()];
            for (Pair<ExactACset, Double> pair : alleleCountSetList) {
                likelihoods[k++] = pair.second - fac;

            }

            return likelihoods;
        }

        public int getLength() {
            return alleleCountSetList.size();
        }
        //public void addLikelihoodList(double[] likelihoods, )
    }

    private static List<Allele> chooseMostLikelyAlternateAlleles(VariantContext vc, int numAllelesToChoose, int ploidy) {
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final LikelihoodSum[] likelihoodSums = new LikelihoodSum[numOriginalAltAlleles];
        for ( int i = 0; i < numOriginalAltAlleles; i++ )
            likelihoodSums[i] = new LikelihoodSum(vc.getAlternateAllele(i));

        // based on the GLs, find the alternate alleles with the most probability; sum the GLs for the most likely genotype
        final ArrayList<double[]> GLs = getGLs(vc.getGenotypes());
        for ( final double[] likelihoods : GLs ) {

            final int PLindexOfBestGL = MathUtils.maxElementIndex(likelihoods);
            final int[] acCount = PoolGenotypeLikelihoods.getAlleleCountFromPLIndex(1+numOriginalAltAlleles,ploidy,PLindexOfBestGL);
            // by convention, first count coming from getAlleleCountFromPLIndex comes from reference allele
            for (int k=1; k < acCount.length;k++) {
                if (acCount[k] > 0)
                    likelihoodSums[k-1].sum += likelihoods[PLindexOfBestGL];

            }
        }

        // sort them by probability mass and choose the best ones
        Collections.sort(Arrays.asList(likelihoodSums));
        final ArrayList<Allele> bestAlleles = new ArrayList<Allele>(numAllelesToChoose);
        for ( int i = 0; i < numAllelesToChoose; i++ )
            bestAlleles.add(likelihoodSums[i].allele);

        final ArrayList<Allele> orderedBestAlleles = new ArrayList<Allele>(numAllelesToChoose);
        for ( Allele allele : vc.getAlternateAlleles() ) {
            if ( bestAlleles.contains(allele) )
                orderedBestAlleles.add(allele);
        }

        return orderedBestAlleles;
    }


    /**
     * Simple non-optimized version that combines GLs from several pools and produces global AF distribution.
     * @param GLs                              Inputs genotypes context with per-pool GLs
     * @param numAlternateAlleles              Number of alternate alleles
     * @param nSamplesPerPool                  Number of samples per pool
     * @param log10AlleleFrequencyPriors       Frequency priors
     * @param result                           object to fill with output values
     */
    protected static void combineSinglePools(final GenotypesContext GLs,
                                               final int numAlternateAlleles,
                                               final int nSamplesPerPool,
                                               final double[] log10AlleleFrequencyPriors,
                                               final AlleleFrequencyCalculationResult result) {

        int[] alleleCounts = new int[numAlternateAlleles];
        final ArrayList<double[]> genotypeLikelihoods = getGLs(GLs);


        int ploidy1 = 0;
        int ploidy2 = 2*nSamplesPerPool;
        // Combine each pool incrementally - likelihoods will be renormalized at each step
        if (numAlternateAlleles > 1) {

            CombinedPoolLikelihoods combinedPoolLikelihoods = new CombinedPoolLikelihoods(numAlternateAlleles);
            combinedPoolLikelihoods.add(0.0, new ExactACset(1, new ExactACcounts(new int[numAlternateAlleles])),0);
            for (int p=1; p<genotypeLikelihoods.size(); p++) {
                result.reset();
                if (USE_NAIVE_COMBINER)
                    combineMultiallelicPoolNaively(combinedPoolLikelihoods, genotypeLikelihoods.get(p), ploidy1, ploidy2,
                            numAlternateAlleles+1, log10AlleleFrequencyPriors, result);
                else
                    fastCombineMultiallelicPool(combinedPoolLikelihoods, genotypeLikelihoods.get(p), ploidy1, ploidy2,
                            numAlternateAlleles + 1, log10AlleleFrequencyPriors, result);
                ploidy1 = ploidy2 + ploidy1; // total number of chromosomes in combinedLikelihoods
            }
        }
        else {
            double[] combinedLikelihoods = new double[1];
            // simple biallelic logic, for efficiency
            for (int p=1; p<genotypeLikelihoods.size(); p++)
                combinedLikelihoods = combineBiallelicPoolNaively(genotypeLikelihoods.get(p), combinedLikelihoods);
            final double log10Lof0 =combinedLikelihoods[0];
            result.setLog10LikelihoodOfAFzero(log10Lof0);
            result.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);

            for (int k=1; k < combinedLikelihoods.length; k++) {
                alleleCounts[0] = k;
                result.updateMLEifNeeded(combinedLikelihoods[k],alleleCounts);
                result.updateMAPifNeeded(combinedLikelihoods[k] + log10AlleleFrequencyPriors[k],alleleCounts);
            }
        }
    }


    public static void fastCombineMultiallelicPool(CombinedPoolLikelihoods xx, double[] yy, int ploidy1, int ploidy2, int numAlleles,
                                                          final double[] log10AlleleFrequencyPriors,
                                                          final AlleleFrequencyCalculationResult result) {

        /*
        maxElement =-Infinity
        for each element in originalPool, with vector originalVector,
         x = xx.likelihood + corr(x)
         for each element in yy (with allele vector currentVector,
           y = yy.likelihood + corr(y)
           double likelihood = computeLofK(set, x,y, log10AlleleFrequencyPriors, numAlleles, ploidy1, ploidy2, result);
           if (likelihood > maxElement - THRESH) {


           }

         */
    }

    /**
     * Naive combiner of two multiallelic pools - number of alt alleles must be the same.
     * Math is generalization of biallelic combiner.
     *
     * For vector K representing an allele count conformation,
     * Pr(D | AC = K) = Sum_G Pr(D|AC1 = G) Pr (D|AC2=K-G) * F(G,K)
     * where F(G,K) = choose(m1,[g0 g1 ...])*choose(m2,[...]) / choose(m1+m2,[k1 k2 ...])
     * @param originalPool                    First log-likelihood pool GL vector
     * @param yy                    Second pool GL vector
     * @param ploidy1               Ploidy of first pool (# of chromosomes in it)
     * @param ploidy2               Ploidy of second pool
     * @param numAlleles            Number of alleles
     * @param log10AlleleFrequencyPriors Array of biallelic priors
     * @param result                Af calculation result object                  
     */
    public static void combineMultiallelicPoolNaively(CombinedPoolLikelihoods originalPool, double[] yy, int ploidy1, int ploidy2, int numAlleles,
                                                          final double[] log10AlleleFrequencyPriors,
                                                          final AlleleFrequencyCalculationResult result) {

        final int dim1 = GenotypeLikelihoods.calculateNumLikelihoods(numAlleles, ploidy1);
        final int dim2 = GenotypeLikelihoods.calculateNumLikelihoods(numAlleles, ploidy2);

        if (dim1 != originalPool.getLength() || dim2 != yy.length)
            throw new ReviewedStingException("BUG: Inconsistent vector length");

        if (ploidy2 == 0)
            return;

        final int newPloidy = ploidy1 + ploidy2;

        // Say L1(K) = Pr(D|AC1=K) * choose(m1,K)
        // and L2(K) = Pr(D|AC2=K) * choose(m2,K)
        PoolGenotypeLikelihoods.SumIterator firstIterator = new PoolGenotypeLikelihoods.SumIterator(numAlleles,ploidy1);
        final double[] x = originalPool.getLikelihoodsAsVector(true);
        while(firstIterator.hasNext()) {
            x[firstIterator.getLinearIndex()] += MathUtils.log10MultinomialCoefficient(ploidy1,firstIterator.getCurrentVector());
            firstIterator.next();
        }

        PoolGenotypeLikelihoods.SumIterator secondIterator = new PoolGenotypeLikelihoods.SumIterator(numAlleles,ploidy2);
        final double[] y = yy.clone();
        while(secondIterator.hasNext()) {
            y[secondIterator.getLinearIndex()] += MathUtils.log10MultinomialCoefficient(ploidy2,secondIterator.getCurrentVector());
            secondIterator.next();
        }

        // initialize output to -log10(choose(m1+m2,[k1 k2...])
        final int outputDim = GenotypeLikelihoods.calculateNumLikelihoods(numAlleles, newPloidy);
        final PoolGenotypeLikelihoods.SumIterator outputIterator = new PoolGenotypeLikelihoods.SumIterator(numAlleles,newPloidy);


        // Now, result(K) =  logSum_G (L1(G)+L2(K-G)) where G are all possible vectors that sum UP to K
        while(outputIterator.hasNext()) {
            final ExactACset set = new ExactACset(1, new ExactACcounts(outputIterator.getCurrentAltVector()));
            double likelihood = computeLofK(set, x,y, log10AlleleFrequencyPriors, numAlleles, ploidy1, ploidy2, result);

            originalPool.add(likelihood, set, outputIterator.getLinearIndex());
            outputIterator.next();
        }
        
    }

    private static double computeLofK(final ExactACset set,
                                    final double[] firstGL, final double[] secondGL,
                                    final double[] log10AlleleFrequencyPriors,
                                    final int numAlleles, final int ploidy1, final int ploidy2,
                                    final AlleleFrequencyCalculationResult result) {

        final int newPloidy = ploidy1 + ploidy2;

        final int totalK = set.getACsum();

        // special case for k = 0 over all k
        if ( totalK == 0 ) {
            final double log10Lof0 = firstGL[HOM_REF_INDEX] + secondGL[HOM_REF_INDEX];
            set.log10Likelihoods[0] = log10Lof0;
                    
            result.setLog10LikelihoodOfAFzero(log10Lof0);
            result.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);

        }   else {

            // initialize result with denominator
            // ExactACset holds by convention only the conformation of alt alleles, since the sum of all allele count is just the ploidy.
            // To compute n!/k1!k2!k3!... we need to compute first n!/(k2!k3!...) and then further divide by k1! where k1=ploidy-sum_k_i

            int[] currentCount = new int[1+set.ACcounts.getCounts().length];
            for (int k=0; k < set.ACcounts.getCounts().length; k++)
                currentCount[k+1] = set.ACcounts.getCounts()[k];

            currentCount[0] = newPloidy - set.getACsum();
            double denom =  -MathUtils.log10MultinomialCoefficient(newPloidy, currentCount);

            // for current conformation, get all possible ways to break vector K into two components G1 and G2
            final PoolGenotypeLikelihoods.SumIterator innerIterator = new PoolGenotypeLikelihoods.SumIterator(numAlleles,ploidy1);
            set.log10Likelihoods[0] = Double.NEGATIVE_INFINITY;
            double innerMax = Double.NEGATIVE_INFINITY;
            while (innerIterator.hasNext()) {
                // check if breaking current conformation into g1 and g2 is feasible.
                final int[] g2 = MathUtils.vectorDiff(currentCount, innerIterator.getCurrentVector());
                if (g2[MathUtils.minElementIndex(g2)] >= 0) {
                    final int idx2 = PoolGenotypeLikelihoods.getLinearIndex(g2,numAlleles,ploidy2);
                    final double sum = firstGL[innerIterator.getLinearIndex()] + secondGL[idx2];

                    set.log10Likelihoods[0] = MathUtils.approximateLog10SumLog10(set.log10Likelihoods[0], sum);
                }
                else {
                    int k=0;
                }
                innerIterator.next();
            }

            set.log10Likelihoods[0] += denom;
        }

        double log10LofK = set.log10Likelihoods[0];

        // update the MLE if necessary
        result.updateMLEifNeeded(log10LofK, set.ACcounts.counts);

        // apply the priors over each alternate allele
        for ( final int ACcount : set.ACcounts.getCounts() ) {
            if ( ACcount > 0 )
                log10LofK += log10AlleleFrequencyPriors[ACcount];
        }
        result.updateMAPifNeeded(log10LofK, set.ACcounts.counts);
        
        return log10LofK;
    }


    /**
     * Naive combiner of two biallelic pools (of arbitrary size).
     * For two pools of size m1 and m2, we can compute the combined likelihood as:
     *   Pr(D|AC=k) = Sum_{j=0}^k Pr(D|AC1=j) Pr(D|AC2=k-j) * choose(m1,j)*choose(m2,k-j)/choose(m1+m2,k)
     * @param xx               GL vector, x[k] = Pr(AC_i = k) for alt allele i
     * @param yy               Second GL vector
     * @return                Combined likelihood vector
     */
    public static double[] combineBiallelicPoolNaively(final double[] xx, final double[] yy) {
     
        final int m1 = xx.length;
        final int m2 = yy.length;
        final double[] result = new double[m1+m2-1];

        /** Pre-fill result array and incorporate weights into input vectors
        * Say L1(k) = Pr(D|AC1=k) * choose(m1,k)
        * and L2(k) = Pr(D|AC2=k) * choose(m2,k)
        * equation reduces to
        * Pr(D|AC=k) = 1/choose(m1+m2,k) * Sum_{j=0}^k L1(k) L2(k-j)
        * which is just plain convolution of L1 and L2 (with pre-existing vector)
           */

        // intialize result vector to log10(1/choose(m1+m2,k)).

        for (int k=0; k < result.length; k++)
            result[k] = - MathUtils.log10BinomialCoefficient(m1+m2,k);


        final double[] x = xx.clone();
        final double[] y = yy.clone();

        for (int k=0; k < x.length; k++)
            x[k] += MathUtils.log10BinomialCoefficient(m1,k);

        for (int k=0; k < y.length; k++)
            y[k] +=  MathUtils.log10BinomialCoefficient(m2,k);


        result[0] += x[0]+y[0];

        double maxElement = result[0];
        int maxElementIdx = 0;
        for (int k=1; k < result.length; k++) {
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

                if (Double.isInfinite(x1) || Double.isInfinite(y1))
                    continue;
                acc[j] = x1+ y1;
            }    
            result[k] += MathUtils.log10sumLog10(acc);
            maxElementIdx = k;
            double maxDiff = result[k] - maxElement;
            if (maxDiff > 0)
                maxElement = result[k];
            else if (maxDiff < maxElement - MAX_LOG10_ERROR_TO_STOP_EARLY) {
                break;
            }


        }
        

        return MathUtils.normalizeFromLog10(Arrays.copyOf(result,maxElementIdx+1),false, true);
    }

    /**
     * From a given variant context, extract a given subset of alleles, and update genotype context accordingly,
     * including updating the PL's, and assign genotypes accordingly
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param allelesToUse                      alleles to subset
     * @param assignGenotypes                   true: assign hard genotypes, false: leave as no-call
     * @param ploidy                            number of chromosomes per sample (pool)
     * @return                                  GenotypesContext with new PLs
     */
    public GenotypesContext subsetAlleles(final VariantContext vc,
                                                      final List<Allele> allelesToUse,
                                                      final boolean assignGenotypes,
                                                      final int ploidy) {
        // the genotypes with PLs
        final GenotypesContext oldGTs = vc.getGenotypes();
        List<Allele> NO_CALL_ALLELES = new ArrayList<Allele>(ploidy);
        
        for (int k=0; k < ploidy; k++)
            NO_CALL_ALLELES.add(Allele.NO_CALL);

        // samples
        final List<String> sampleIndices = oldGTs.getSampleNamesOrderedByName();

        // the new genotypes to create
        final GenotypesContext newGTs = GenotypesContext.create();

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;


        // create the new genotypes
        for ( int k = 0; k < oldGTs.size(); k++ ) {
            final Genotype g = oldGTs.get(sampleIndices.get(k));
            if ( !g.hasLikelihoods() ) {
                newGTs.add(new Genotype(g.getSampleName(), NO_CALL_ALLELES, Genotype.NO_LOG10_PERROR, null, null, false));
                continue;
            }

            // create the new likelihoods array from the alleles we are allowed to use
            final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
            double[] newLikelihoods;
            if ( numOriginalAltAlleles == numNewAltAlleles) {
                newLikelihoods = originalLikelihoods;
            } else {
                newLikelihoods = PoolGenotypeLikelihoods.subsetToAlleles(originalLikelihoods, ploidy, vc.getAlleles(),allelesToUse);

                // might need to re-normalize
                newLikelihoods = MathUtils.normalizeFromLog10(newLikelihoods, false, true);
            }

            // if there is no mass on the (new) likelihoods, then just no-call the sample
            if ( MathUtils.sum(newLikelihoods) > VariantContextUtils.SUM_GL_THRESH_NOCALL ) {
                newGTs.add(new Genotype(g.getSampleName(), NO_CALL_ALLELES, Genotype.NO_LOG10_PERROR, null, null, false));
            }
            else {
                Map<String, Object> attrs = new HashMap<String, Object>(g.getAttributes());
                if ( numNewAltAlleles == 0 )
                    attrs.remove(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY);
                else
                    attrs.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods(newLikelihoods));

                // if we weren't asked to assign a genotype, then just no-call the sample
                if ( !assignGenotypes || MathUtils.sum(newLikelihoods) > VariantContextUtils.SUM_GL_THRESH_NOCALL )
                    newGTs.add(new Genotype(g.getSampleName(), NO_CALL_ALLELES, Genotype.NO_LOG10_PERROR, null, attrs, false));
                else
                    newGTs.add(assignGenotype(g, newLikelihoods, allelesToUse, ploidy, attrs));
            }
        }

        return newGTs;

    }
    /**
     * Assign genotypes (GTs) to the samples in the Variant Context greedily based on the PLs
     *
     * @param originalGT           the original genotype
     * @param newLikelihoods       the PL array
     * @param allelesToUse         the list of alleles to choose from (corresponding to the PLs)
     * @param numChromosomes        Number of chromosomes per pool
     * @param attrs                the annotations to use when creating the genotype
     *
     * @return genotype
     */
    private static Genotype assignGenotype(final Genotype originalGT, final double[] newLikelihoods, final List<Allele> allelesToUse, 
                                           final int numChromosomes, final Map<String, Object> attrs) {
        final int numNewAltAlleles = allelesToUse.size() - 1;
        


        // find the genotype with maximum likelihoods
        final int PLindex = numNewAltAlleles == 0 ? 0 : MathUtils.maxElementIndex(newLikelihoods);

        final int[] mlAlleleCount = PoolGenotypeLikelihoods.getAlleleCountFromPLIndex(allelesToUse.size(), numChromosomes, PLindex);
        final ArrayList<String> alleleFreqs = new ArrayList<String>();
        final ArrayList<Integer> alleleCounts = new ArrayList<Integer>();


        for (int k=1; k < mlAlleleCount.length; k++) {
            alleleCounts.add(mlAlleleCount[k]);
            final String freq = String.format(VariantContextUtils.makePrecisionFormatStringFromDenominatorValue((double)numChromosomes), ((double)mlAlleleCount[k] / (double)numChromosomes));
            alleleFreqs.add(freq);
            
        }

        // per-pool logging of AC and AF
        attrs.put(MAXIMUM_LIKELIHOOD_AC_KEY, alleleCounts.size() == 1 ? alleleCounts.get(0) : alleleCounts);
        attrs.put(MAXIMUM_LIKELIHOOD_AF_KEY, alleleFreqs.size() == 1 ? alleleFreqs.get(0) : alleleFreqs);
        
        if (newLikelihoods.length > MAX_LENGTH_FOR_POOL_PL_LOGGING)
            attrs.remove(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY);
        ArrayList<Allele> myAlleles = new ArrayList<Allele>();

        // add list of called ML genotypes to alleles list
        // TODO - too unwieldy?
        int idx = 0;
        for (int mlind = 0; mlind < mlAlleleCount.length; mlind++) {
            for (int k=0; k < mlAlleleCount[mlind]; k++)
                myAlleles.add(idx++,allelesToUse.get(mlind));
        }

        final double qual = numNewAltAlleles == 0 ? Genotype.NO_LOG10_PERROR : GenotypeLikelihoods.getQualFromLikelihoods(PLindex, newLikelihoods);
        return new Genotype(originalGT.getSampleName(), myAlleles, qual, null, attrs, false);
    }

}
