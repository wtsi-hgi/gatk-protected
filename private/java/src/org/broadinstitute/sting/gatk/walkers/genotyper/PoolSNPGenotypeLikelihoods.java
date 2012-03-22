package org.broadinstitute.sting.gatk.walkers.genotyper;


import net.sf.samtools.SAMUtils;
import org.broadinstitute.sting.gatk.walkers.poolcaller.ErrorModel;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.fragments.FragmentCollection;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

/**
 * Stable, error checking version of the Bayesian genotyper.  Useful for calculating the likelihoods, priors,
 * and posteriors given a pile of bases and quality scores
 *
 * Suppose we have bases b1, b2, ..., bN with qualities scores q1, q2, ..., qN.  This object
 * calculates:
 *
 * P(G | D) = P(G) * P(D | G)
 *
 * where
 *
 * P(D | G) = sum_i log10 P(bi | G)
 *
 * and
 *
 * P(bi | G) = 1 - P(error | q1) if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for homozygous genotypes and for heterozygous genotypes:
 *
 * P(bi | G) = 1 - P(error | q1) / 2 + P(error | q1) / 6 if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for each of the 10 unique diploid genotypes AA, AC, AG, .., TT
 *
 * Everything is stored as arrays indexed by DiploidGenotype.ordinal() values in log10 space.
 *
 * The priors contain the relative probabilities of each genotype, and must be provided at object creation.
 * From then on, you can call any of the add() routines to update the likelihoods and posteriors in the above
 * model.
 */
public class PoolSNPGenotypeLikelihoods implements Cloneable {
    public final static double DEFAULT_PCR_ERROR_RATE = 1e-4;

    protected final static int FIXED_PLOIDY = 2;
    protected final static int MAX_PLOIDY = FIXED_PLOIDY + 1;
    protected final static double ploidyAdjustment = log10(FIXED_PLOIDY);
    protected final static double log10_3 = log10(3.0);
    protected final int TWO_N;

    protected boolean VERBOSE = false;

    //
    // The fundamental data arrays associated with a Genotype Likelhoods object
    //
    protected double[] log10Likelihoods = null;
    protected double[] log10Posteriors = null;
    private final double[] genotypeZeros;

    protected PoolGenotypePriors priors = null;

    // TODO: don't calculate this each time through
    protected double log10_PCR_error_3;
    protected double log10_1_minus_PCR_error;
    protected final int nSamplesPerPool;
    HashMap<String, ErrorModel> perLaneErrorModels;

    protected final boolean treatAllReadsAsSinglePool;


    double[][][] logPLs;
    private double[][] logMismatchProbabilityArray;    
    static private final double qualVec[] = new double[SAMUtils.MAX_PHRED_SCORE+1];
    //RecursivePLHolder plHolder;

        /**
        * Create a new GenotypeLikelhoods object with given priors and PCR error rate for each pool genotype
        *
        * @param priors          priors
        * @param PCR_error_rate  the PCR error rate
        */
    public PoolSNPGenotypeLikelihoods(PoolGenotypePriors priors, double PCR_error_rate, HashMap<String, ErrorModel> perLaneErrorModels, boolean treatAllReadsAsSinglePool) {
        this.priors = priors;
        log10_PCR_error_3 = log10(PCR_error_rate) - log10_3;
        log10_1_minus_PCR_error = log10(1.0 - PCR_error_rate);
        nSamplesPerPool = priors.getNSamplesPerPool();
        TWO_N = 2*nSamplesPerPool;
        genotypeZeros = new double[priors.getPriors().length];
        setToZero();
        this.perLaneErrorModels = perLaneErrorModels;
        this.treatAllReadsAsSinglePool = treatAllReadsAsSinglePool;
        fillCache();

        //createLikelihoodDataStorage();
        //plHolder = new RecursivePLHolder(TWO_N, BaseUtils.BASES.length);
    }

    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        PoolSNPGenotypeLikelihoods c = (PoolSNPGenotypeLikelihoods)super.clone();
        c.priors = priors;
        c.log10Likelihoods = log10Likelihoods.clone();
        c.log10Posteriors = log10Posteriors.clone();
        return c;
    }
        //
    protected void setToZero() {
        log10Likelihoods = genotypeZeros.clone();                 // likelihoods are all zeros
        log10Posteriors = priors.getPriors().clone();     // posteriors are all the priors
        logPLs = new double[1+TWO_N][1+TWO_N][1+TWO_N];
    }

    /**
     * Returns an array of log10 likelihoods for each genotype, indexed by DiploidGenotype.ordinal values()
     * @return likelihoods array
     */
    public double[] getLikelihoods() {
        return log10Likelihoods;
    }

    /**
     * Returns an array of posteriors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return raw log10 (not-normalized posteriors) as a double array
     */
    public double[] getPosteriors() {
        return log10Posteriors;
    }



    /**
     * Returns an array of posteriors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return normalized posterors as a double array
     */
/*    public double[] getNormalizedPosteriors() {
        double[] normalized = new double[log10Posteriors.length];
        double sum = 0.0;

        // for precision purposes, we need to add (or really subtract, since everything is negative)
        // the largest posterior value from all entries so that numbers don't get too small
        double maxValue = log10Posteriors[0];
        for (int i = 1; i < log10Posteriors.length; i++) {
            if ( maxValue < log10Posteriors[i] )
                maxValue = log10Posteriors[i];
        }

        // collect the posteriors
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double posterior = Math.pow(10, getPosterior(g) - maxValue);
            normalized[g.ordinal()] = posterior;
            sum += posterior;
        }

        // normalize
        for (int i = 0; i < normalized.length; i++)
            normalized[i] /= sum;

        return normalized;
    }
  */


    /**
     * Returns an array of priors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return log10 prior as a double array
     */
    public double[] getPriors() {
        return priors.getPriors();
    }

    // -------------------------------------------------------------------------------------
    //
    // add() routines.  These are the workhorse routines for calculating the overall genotype
    // likelihoods given observed bases and reads.  Includes high-level operators all the
    // way down to single base and qual functions.
    //
    // -------------------------------------------------------------------------------------

    /**
     * Updates likelihoods and posteriors to reflect the additional observations contained within the
     * read-based pileup up by calling add(observedBase, qualityScore) for each base / qual in the
     * pileup
     *
     * @param pileup                    read pileup
     * @param ignoreBadBases            should we ignore bad bases?
     * @param capBaseQualsAtMappingQual should we cap a base's quality by its read's mapping quality?
     * @param minBaseQual               the minimum base quality at which to consider a base valid
     * @return the number of good bases found in the pileup
     */
    public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        int n = 0;


        for (String laneID : perLaneErrorModels.keySet() ) {
            // get pileup for this lane
            ReadBackedPileup perLanePileup;
            if (treatAllReadsAsSinglePool)
                perLanePileup = pileup;
            else
                perLanePileup = pileup.getPileupForLane(laneID);

            if (perLanePileup == null || perLanePileup.isEmpty())
                continue;

            ErrorModel errorModel = perLaneErrorModels.get(laneID);
            n += add(perLanePileup, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual, errorModel);
            if (treatAllReadsAsSinglePool)
                break;
            
        }
        log10Likelihoods = convertToVectorIndex(logPLs, TWO_N);
        return n;
    }

    public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual, ErrorModel errorModel) {
        int n=0;
        int[] numSeenBases = new int[BaseUtils.BASES.length];

        // count number of elements in pileup
        for (PileupElement elt : pileup) {
            byte obsBase = elt.getBase();
            byte qual = qualToUse(elt, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);
            if ( qual == 0 )
                continue;
            int idx = 0;
            for (byte base:BaseUtils.BASES)
                numSeenBases[idx++] += (base == obsBase?1:0);


            n++;

        }
        add(numSeenBases[0], numSeenBases[1], numSeenBases[2], numSeenBases[3], errorModel);
        return n;
    }
/*    public int add(List<PileupElement> overlappingPair, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        final PileupElement p1 = overlappingPair.get(0);
        final PileupElement p2 = overlappingPair.get(1);
        String laneID = PoolCallerEngine.getLaneIDFromReadGroupString(p1.getRead().getReadGroup().getReadGroupId(), treatAllReadsAsSinglePool);
        ErrorModel errorModel = perLaneErrorModels.get(laneID);
        // todo - are pairs in a fragment always sequenced in same lane?

        final byte observedBase1 = p1.getBase();
        final byte qualityScore1 = qualToUse(p1, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);
        final byte observedBase2 = p2.getBase();
        final byte qualityScore2 = qualToUse(p2, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);

        if ( qualityScore1 == 0 ) {
            if ( qualityScore2 == 0 ) // abort early if we didn't see any good bases
                return 0;
            else {
                return add(observedBase2, qualityScore2, (byte)0, (byte)0, errorModel);
            }
        } else {
            return add(observedBase1, qualityScore1, observedBase2, qualityScore2, errorModel);
        }
    }     */

    /**
     * Calculates the pool's probability for all possible allele counts for all bases. Calculation is based on the error model
     * generated by the reference sample on the same lane. The probability is given by :
     *
     * Pr(ac=jA,jC,jG,jT| pool, errorModel) = sum_over_all_Qs ( Pr(ac=jA,jC,jG,jT) * Pr(errorModel_q) *
     * [jA*(1-eq)/2n + eq/3*(jc+jg+jt)/2N)^nA *   jC*(1-eq)/2n + eq/3*(ja+jg+jt)/2N)^nC *
     * jG*(1-eq)/2n + eq/3*(jc+ja+jt)/2N)^nG * jT*(1-eq)/2n + eq/3*(jc+jg+ja)/2N)^nT
     *
     *  log Pr(ac=jA,jC,jG,jT| pool, errorModel) = logsum( Pr(ac=jA,jC,jG,jT) * Pr(errorModel_q) *
     * [jA*(1-eq)/2n + eq/3*(jc+jg+jt)/2N)^nA *   jC*(1-eq)/2n + eq/3*(ja+jg+jt)/2N)^nC *
     * jG*(1-eq)/2n + eq/3*(jc+ja+jt)/2N)^nG * jT*(1-eq)/2n + eq/3*(jc+jg+ja)/2N)^nT)
     * = logsum(logPr(ac=jA,jC,jG,jT) + log(Pr(error_Model(q)
     * )) + nA*log(jA/2N(1-eq)+eq/3*(2N-jA)/2N) + nC*log(jC/2N(1-eq)+eq/3*(2N-jC)/2N)
     * + log(jG/2N(1-eq)+eq/3*(2N-jG)/2N) + log(jT/2N(1-eq)+eq/3*(2N-jT)/2N)
     *
     * Let Q(j,k) = log(j/2N*(1-e[k]) + (2N-j)/2N*e[k]/3)
     *
     * Then logPr(ac=jA,jC,jG,jT|D,errorModel) = logPR(ac=Ja,jC,jG,jT) + logsum_k( logPr (errorModel[k],
     * nA*Q(jA,k) +  nC*Q(jC,k) + nG*Q(jG,k) + nT*Q(jT,k))
     *
     * If pileup data comes from several error models (because lanes can have different error models),
     * Pr(Ac=j|D,E1,E2) = sum(Pr(AC1=j1|D,E1,E2) * Pr(AC2=j-j2|D,E1,E2))
     * = sum(Pr(AC1=j1|D,E1)*Pr(AC2=j-j1|D,E2)) from j=0..2N
     *
     * So, for each lane, build error model and combine lanes.
     * To store model, can do
     * for jA=0:2N
     *  for jC = 0:2N-jA
     *   for jG = 0:2N-jA-jC
     *    for jT = 0:2N-jA-jC-jG
     *      Q(jA,jC,jG,jT)
     *      for k = minSiteQual:maxSiteQual
     *        likelihood(jA,jC,jG,jT) = logsum(logPr (errorModel[k],nA*Q(jA,k) +  nC*Q(jC,k) + nG*Q(jG,k) + nT*Q(jT,k))
     *
     *
     *
     * where:
      */
    private int add(int nA, int nC, int nG, int nT, ErrorModel errorModel) {
        int minQ = errorModel.getMinQualityScore();
        int maxQ = errorModel.getMaxQualityScore();
        for (int jA = 0; jA <= TWO_N; jA++) {
            for (int jC = 0; jC <= TWO_N - jA; jC++) {
                for (int jG = 0; jG <= TWO_N - jA - jC; jG++) {
                    int jT = TWO_N - jA - jC - jG;
                    // for observed base X, add Q(jX,k) to likelihood vector for all k in error model
                    //likelihood(jA,jC,jG,jT) = logsum(logPr (errorModel[k],nA*Q(jA,k) +  nC*Q(jC,k) + nG*Q(jG,k) + nT*Q(jT,k))
                    double[] acVec = new double[maxQ - minQ + 1];

                    for (int k=minQ; k<=maxQ; k++)
                        acVec[k-minQ] = nA*logMismatchProbabilityArray[jA][k] +
                                nC*logMismatchProbabilityArray[jC][k] +
                                nG*logMismatchProbabilityArray[jG][k] +
                                nT*logMismatchProbabilityArray[jT][k];

                    logPLs[jA][jC][jG] += MathUtils.logDotProduct(errorModel.getErrorModelVector(), acVec);
                }
            }
        }


        return 1;
    }

    public static double[] convertToVectorIndex(final double[][][] plMatrix, final int TWO_N) {
        // sum_0 to M i = M*(M+1)/2
        // sum_j=0 to j=(M-jA) (M-jA-j) = (M-jA)(M-jA+1)-(M-jA)(M-jA+1)/2 = (M-jA)(M-jA+1)/2
        // todo -closed form solution to this!!
        int len =0;
        for (int jA = 0; jA <= TWO_N; jA++) {
            for (int jC = 0; jC <= TWO_N - jA; jC++) {
                len += 1+TWO_N - jA - jC;
              }
        }

        //int len = (BaseUtils.BASES.length) * (BaseUtils.BASES.length +1)/2;

        double plVec[] = new double[len];

        int idx = len-1;
        for (int jA = 0; jA <= TWO_N; jA++) {
            for (int jC = 0; jC <= TWO_N - jA; jC++) {
                for (int jG = 0; jG <= TWO_N - jA - jC; jG++) {
                    // reverse ordering to match DiploidGenotype ordering for SNP case, i.e.
                    // ordering for N=1 would be AA,AC,AG,AT,CC,CG,CT,GG,GT
                    plVec[idx--] = plMatrix[jA][jC][jG];

                }
            }
        }
        return plVec;
    }

     /**
     * Helper function that returns the phred-scaled base quality score we should use for calculating
     * likelihoods for a pileup element.  May return 0 to indicate that the observation is bad, and may
     * cap the quality score by the mapping quality of the read itself.
     *
     * @param p
     * @param ignoreBadBases
     * @param capBaseQualsAtMappingQual
     * @param minBaseQual
     * @return
     */
    private static byte qualToUse(PileupElement p, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        if ( ignoreBadBases && !BaseUtils.isRegularBase( p.getBase() ) )
            return 0;

        byte qual = p.getQual();

        if ( qual > SAMUtils.MAX_PHRED_SCORE )
            throw new UserException.MalformedBAM(p.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, p.getRead().getReadName()));
        if ( capBaseQualsAtMappingQual )
            qual = (byte)Math.min((int)qual, p.getMappingQual());
        if ( (int)qual < minBaseQual )
            qual = (byte)0;

        return qual;
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // helper routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Return a string representation of this object in a moderately usable form
     *
     * @return string representation
     */
    public String toString() {
        double sum = 0;
        StringBuilder s = new StringBuilder();
        for (int i=0; i < priors.getPriors().length; i++) {
            s.append(String.format("%d %.10f ", i, log10Likelihoods[i]));
            sum += Math.pow(10,log10Likelihoods[i]);
        }
        s.append(String.format(" %f", sum));
        return s.toString();
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // Validation routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    public boolean validate() {
        return validate(true);
    }

    public boolean validate(boolean throwException) {
        try {
            priors.validate(throwException);

            for (int i=0; i < priors.getPriors().length; i++) {
                String bad = null;

                if ( ! MathUtils.wellFormedDouble(log10Likelihoods[i]) || ! MathUtils.isNegativeOrZero(log10Likelihoods[i]) ) {
                    bad = String.format("Likelihood %f is badly formed", log10Likelihoods[i]);
                } else if ( ! MathUtils.wellFormedDouble(log10Posteriors[i]) || ! MathUtils.isNegativeOrZero(log10Posteriors[i]) ) {
                    bad = String.format("Posterior %f is badly formed", log10Posteriors[i]);
                }

                if ( bad != null ) {
                    throw new IllegalStateException(String.format("At %d: %s", i, bad));
                }
            }
        } catch ( IllegalStateException e ) {
            if ( throwException )
                throw new RuntimeException(e);
            else
                return false;
        }

        return true;
    }

    //
    // Constant static data
    //
    private final static double[] baseZeros = new double[BaseUtils.BASES.length];

    static {
        for ( byte base : BaseUtils.BASES ) {
            baseZeros[BaseUtils.simpleBaseToBaseIndex(base)] = 0.0;
        }
    }
       
      /*
    class RecursivePLHolder {
        int size;
        int recursionLevel;
        Object[] pls;
        RecursivePLHolder(int size, int recursionLevel) {
            this.size = size;
            this.recursionLevel = recursionLevel;
            if (recursionLevel> 1) {
                pls = new RecursivePLHolder[size];
                for (int k=0; k < size; k++) {
                    pls[k] = new RecursivePLHolder(size-k,recursionLevel-1);
                }
            }
            else {
                pls = new Double[TWO_N];
            }
        }
        
        Object get(int k) {
            return pls[k];
        }
        void set(int[] indx, Double val[]) {
            if (indx.length != recursionLevel)
                throw new ReviewedStingException("BUG: inconsistent vector size with recursion levels");
 
            if (recursionLevel == 1)
                pls[indx[0]] = val;
            else {
                RecursivePLHolder s = (RecursivePLHolder)pls[indx[0]]; 
                int[] subIdx = new int[indx.length-1];
                for (int k=1; k < indx.length; k++)
                    subIdx[k-1] = indx[k];
                s.set(subIdx, val);
            }
        }
    }   */


    static {
        // cache 10^(-k/10)
        for (int j=0; j <= SAMUtils.MAX_PHRED_SCORE; j++)
            qualVec[j] = Math.pow(10.0,-(double)j/10.0);
    }

    private void fillCache() {
        // cache Q(j,k) = log10(j/2N*(1-ek) + (2N-j)/2N*ek) for j = 0:2N

        logMismatchProbabilityArray = new double[1+TWO_N][1+SAMUtils.MAX_PHRED_SCORE];
        for (int i=0; i <= TWO_N; i++) {
            for (int j=0; j <= SAMUtils.MAX_PHRED_SCORE; j++) {
                double phi = (double)i/TWO_N;
                logMismatchProbabilityArray[i][j] = Math.log10(phi * (1.0-qualVec[j]) + qualVec[j]/3.0 * (1.0-phi));
            }
        }
    }

}
