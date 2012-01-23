package org.broadinstitute.sting.gatk.walkers.poolcaller;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.picard.util.MathUtil;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.ArrayList;

/**
 * The allele count model for the pool calling framework.
 *
 * @author carneiro
 * @since 7/21/11
 */
public class AlleleCountModel extends ProbabilityModel implements Cloneable {
    private int maxAlleleCount;
    private double minCallQuall;
    private ErrorModel errorModel;
    private final double THETA = 0.001; // Human heterozygozity rate
    private byte bestAltBase;
    private byte refBase;
    ArrayList<Allele> altAlleles = new ArrayList<Allele>();
    private double phredScaledConfidence;
    private int idxOfMaxAC;
  //  private double[][] normalizedPosteriors;

    public AlleleCountModel(int maxAlleleCount, ErrorModel errorModel, int matches, int mismatches, double minCallQual) {
        this.maxAlleleCount = maxAlleleCount;
        this.minCallQuall = minCallQual;
        this.errorModel = errorModel;
        model = calculateLog10ProbabilityDistribution(errorModel, matches, mismatches);
    }

    public AlleleCountModel(int maxAlleleCount, ErrorModel errorModel, Integer[] numSeenBases,  double minCallQual, byte refBase) {
        this.maxAlleleCount = maxAlleleCount;
        this.minCallQuall = minCallQual;
        this.errorModel = errorModel;
        this.refBase = refBase;
        model = calculateLog10ProbabilityDistributionForAllBases(errorModel, numSeenBases, refBase);

    }
 /*   public AlleleCountModel(AlleleCountModel acm) {
        maxAlleleCount = acm.getMaxAlleleCount();
        model = new double [maxAlleleCount+1];
        for (int ac=0; ac<=maxAlleleCount; ac++) {
            model[ac] = acm.getLog10ProbabilityGivenAC(ac);
        }
    }
       */
    public int getMaxAlleleCount() {
        return maxAlleleCount;
    }


    public void merge(AlleleCountModel mergeModel) {
        // I need to know which model is the largest so I can iterate on him (the j index) and guarantee that
        // once we run out of terms it is because the smaller one (the k index) has exhausted all it's terms
        AlleleCountModel largerModel = (maxAlleleCount >= mergeModel.getMaxAlleleCount()) ? this : mergeModel;
        AlleleCountModel smallerModel = (maxAlleleCount < mergeModel.getMaxAlleleCount()) ? this : mergeModel;

        double [] result = new double [largerModel.getMaxAlleleCount() + smallerModel.getMaxAlleleCount() + 1];
        for (int i=0; i<result.length; i++) {
            double [] sumTerms = new double[i+1];
            int index = 0;
            for (int j= Math.min(i, largerModel.getMaxAlleleCount()); j>=0; j--) {
                int k = Math.min(i - j, smallerModel.getMaxAlleleCount());
                // For models of different size, some combinations are impossible. Once we reach one, no reason in trying any further
                // it means that smaller model is not big enough to keep generating sumTerms.
                if (j + k < i)
                    break;

                sumTerms[index] = largerModel.getLog10ProbabilityGivenAC(j) + smallerModel.getLog10ProbabilityGivenAC(k);
                index++;
            }
        result[i] = MathUtils.log10sumLog10(sumTerms);
        }
        model = result;
    }

    @Requires({"ac>=0", "ac<=maxAlleleCount"})
    @Ensures({"result <= 0", "! Double.isInfinite(result)", "! Double.isNaN(result)"})
    public double getLog10ProbabilityGivenAC(int ac) {
        return model[ac];
    }

    @Requires({"ac>=0", "ac<=maxAlleleCount", "log10p <= 0", "! Double.isInfinite(log10p)", "! Double.isNaN(log10p)"})
    public void setLog10ProbabilityGivenAC(int ac, double log10p) {
        model[ac] = log10p;
    }


    /**
     * The prior probability of an allele count being observed based solely on the human heterozygozity rate
     * and the number of samples
     *
     * @param ac given allele count
     * @return the prior probability of ac
     */
    public double log10PriorAC (int ac) {
        // prior probability for a given allele count is THETA/AC.
        if (ac > 0)
            return Math.log10(THETA/ac);

        // if allele count is 0, the prior is one minus the sum of all other possibilities.
        double result = 0.0;
        for (int i=1; i<=maxAlleleCount; i++)
            result += THETA/i;
        return Math.log10(1-result);
    }

    public double log10PriorAC(Integer[] acVector, byte refBase) {
        int refBaseIndex = BaseUtils.simpleBaseToBaseIndex(refBase);

        // todo - generalize later for multiallelic priors.
        int numNZElements = 0;
        int ac = 0;
        for (int k=0; k < acVector.length; k++) {
            if (acVector[k] > 0)
                numNZElements++;

            if (k == refBaseIndex)
                continue;
            ac += acVector[k];
        }

        // sanity check: only 2 non-zero elements,
        if (numNZElements > 2)
            return Double.NEGATIVE_INFINITY;

        // ac contains sum of all non-ref elements in acVector
        return  log10PriorAC(ac);
    }

    /**
     * Calculates the pool's probability for all possible allele counts. Calculation is based on the error model
     * generated by the reference sample on the same lane. The probability is given by :
     *
     * Pr(ac=j | pool, errorModel) = sum_over_all_Qs ( Pr(ac=j) * Pr(errorModel_q) * [ (n-j/2n) * (1-e) + (j/n)*e]^m * [(n-j/n)*e + (j/n) * (1-e)]^(1-m)
     *
     * where:
     *  n = number of chromosomes
     *  e = probability of an error at a given Q level (e.g. Q30 = 0.001, Q20 = 0.01, ...)
     *  m = number of mismatches
     *
     * @param errorModel
     * @param matches
     * @param mismatches
     * @return
     */
    public double [] calculateLog10ProbabilityDistribution (ErrorModel errorModel, int matches, int mismatches) {
        double [] p = new double[maxAlleleCount+1];
        for (int ac=0; ac<=maxAlleleCount; ac++) {
            p[ac] = calculateLog10ProbabilityGivenAC(ac, errorModel, matches, mismatches);
        }
        return p;
    }

    /**
     * Calculates the pool's probability for all possible allele counts for all bases. Calculation is based on the error model
     * generated by the reference sample on the same lane. The probability is given by :
     *
     * Pr(ac=jA,jC,jG,jT| pool, errorModel) = sum_over_all_Qs ( Pr(ac=jA,jC,jG,jT) * Pr(errorModel_q) *
     * [jA*(1-eq)/2n + eq/3*(jc+jg+jt)/2N)^nA *   jC*(1-eq)/2n + eq/3*(ja+jg+jt)/2N)^nC *
     * jG*(1-eq)/2n + eq/3*(jc+ja+jt)/2N)^nG * jT*(1-eq)/2n + eq/3*(jc+jg+ja)/2N)^nT
     *
     *
     * where:
     *  n = number of chromosomes
     *  e = probability of an error at a given Q level (e.g. Q30 = 0.001, Q20 = 0.01, ...)
     *  m = number of mismatches
     *
     * @param errorModel
     * @param numSeenBases
     * @param refBase
     * @return
     */
    public double [] calculateLog10ProbabilityDistributionForAllBases(ErrorModel errorModel, Integer[] numSeenBases, byte refBase) {
        double [][] p = new double[BaseUtils.BASES.length][maxAlleleCount+1];
        int refPos = BaseUtils.simpleBaseToBaseIndex(refBase);

        double[] pMax = new double[BaseUtils.BASES.length];
        int[] mlEstimatorIdx = new int[BaseUtils.BASES.length];

        for (int altInd = 0; altInd < BaseUtils.BASES.length; altInd++) {
            if (altInd == refPos) {
                mlEstimatorIdx[altInd] = -1;
                pMax[altInd] = Double.NEGATIVE_INFINITY;
                continue;
            }

            for (int ac=0; ac<=maxAlleleCount; ac++) {
                // for each alt allele, get AC distribution based on pileup
                // todo - assume for now true multiallelic case has prior = 0
                Integer[] acVectorToCompute = new Integer[BaseUtils.BASES.length];
                // todo - apparent redundant initialization but an array of Integer[] is initialized by null references, not zeros.
                for (int k=0; k < BaseUtils.BASES.length; k++)
                      acVectorToCompute[k] = 0;
                acVectorToCompute[refPos] = maxAlleleCount - ac;
                acVectorToCompute[altInd] = ac;

                p[altInd][ac] = calculateLog10ProbabilityGivenACforAllBases(errorModel, acVectorToCompute, numSeenBases, refBase);
            }
            //normalizedPosteriors[altInd] = MathUtils.normalizeFromLog10(p[altInd],true);
            mlEstimatorIdx[altInd] = MathUtils.maxElementIndex(p[altInd]);
            pMax[altInd] = p[altInd][mlEstimatorIdx[altInd]];
        }


        int bestAlleleIdx = MathUtils.maxElementIndex(pMax);
        idxOfMaxAC = mlEstimatorIdx[bestAlleleIdx];

        if (idxOfMaxAC == 0) {
            // best alt AC is zero, so this is a reference site
            bestAltBase = refBase;
        }
        else {
            bestAltBase = BaseUtils.BASES[bestAlleleIdx];
        }

        //altAlleles =
        //return p[bestAlleleIdx];
       // is the most likely frequency conformation AC=0 for all alternate alleles?

        // determine which alternate alleles have AF>0

        // calculate p(f>0):
        final double[] normalizedPosteriors = MathUtils.normalizeFromLog10(p[bestAlleleIdx]);
        final double PofF = 1.0 - normalizedPosteriors[0];
 /*
        if ( !bestGuessIsRef || UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(normalizedPosteriors[0]);
            if ( Double.isInfinite(phredScaledConfidence) )
                phredScaledConfidence = -10.0 * AFresult.log10PosteriorOfAFzero;
        } else {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(PofF);

  */

        if (idxOfMaxAC == 0) {
            // ref call

            phredScaledConfidence = QualityUtils.phredScaleErrorRate(PofF);
            if ( Double.isInfinite(phredScaledConfidence) ) {
                double sum =  0.0;
                for (int i = 1; i < normalizedPosteriors.length; i++) {
                    sum += p[bestAlleleIdx][i];
                }
                phredScaledConfidence = (MathUtils.compareDoubles(sum, 0.0) == 0 ? 0 : -10.0 * sum);
            }
        }
        else {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(normalizedPosteriors[0]);
            if ( Double.isInfinite(phredScaledConfidence) )
                phredScaledConfidence = -10.0 * p[bestAlleleIdx][0];
            // alt call
        }
        return p[bestAlleleIdx];
    }




    @Requires({"ac >= 0", "ac <= maxAlleleCount", "errorModel != null", "matches >= 0", "mismatches >= 0" })
    @Ensures({"result <= 0", "! Double.isInfinite(result)", "! Double.isNaN(result)"})
    private double calculateLog10ProbabilityGivenAC(int ac, ErrorModel errorModel, int matches, int mismatches) {
        double log10PAC = log10PriorAC(ac);

        // for each quality probability in the model, calculate the probability of the allele count = ac
        // we skip Q0 because it's meaningless.
        double [] acc = new double[errorModel.size()]; // we're skipping Q0 so we don't need maxQualityScore + 1 here.
        for (byte qual = errorModel.getMinQualityScore(); qual <= errorModel.getMaxQualityScore(); qual++) {
            final int i = qual - errorModel.getMinQualityScore();
            acc[i] = calculateLog10ProbabilityOfACGivenQual(qual, ac, errorModel, matches, mismatches);
        }
        return log10PAC + MathUtils.log10sumLog10(acc);
    }


    private double calculateLog10ProbabilityGivenACforAllBases(ErrorModel errorModel, Integer[] acVector,  Integer[] numSeenBases, byte refBase) {
        double   log10PAC = log10PriorAC(acVector, refBase);
        // for each quality probability in the model, calculate the probability of the allele count = ac
        // we skip Q0 because it's meaningless.
        double [] acc = new double[errorModel.size()]; // we're skipping Q0 so we don't need maxQualityScore + 1 here.
        for (byte qual = errorModel.getMinQualityScore(); qual <= errorModel.getMaxQualityScore(); qual++) {
             final int i = qual - errorModel.getMinQualityScore();
            acc[i] = calculateLog10ProbabilityOfACGivenQual(qual, acVector, errorModel,  numSeenBases, refBase);
        }
        return log10PAC + MathUtils.log10sumLog10(acc);
    }


    @Requires({"ac >= 0", "ac <= maxAlleleCount", "errorModel != null", "qual >= errorModel.getMinQualityScore()", "qual <= errorModel.getMaxQualityScore()", "matches >= 0", "mismatches >= 0" })
    @Ensures({"result <= 0", "! Double.isInfinite(result)", "! Double.isNaN(result)"})
    private double calculateLog10ProbabilityOfACGivenQual(byte qual, int ac, ErrorModel errorModel, int matches, int mismatches) {
        double p = (double) ac / maxAlleleCount;
        double q = 1 - p;
        double e = MathUtils.phredScaleToProbability(qual);
        double x = Math.log10(q * (1-e) + p * e);
        double y = Math.log10(q * e + p * (1-e));
        return errorModel.getErrorGivenQual(qual) + matches * x + mismatches * y;
    }

    private double calculateLog10ProbabilityOfACGivenQual(byte qual, Integer[] acVector, ErrorModel errorModel, Integer[] numSeenBases, byte refBase) {
        double e = MathUtils.phredScaleToProbability(qual);
        double t1 = (1-e)/maxAlleleCount;
        double t2 = e/(3*maxAlleleCount);

        Integer[] fnov = new Integer[BaseUtils.BASES.length];
        fnov[0] = acVector[1] + acVector[2] + acVector[3];//numC + numG + numT;
        fnov[1] = acVector[0] + acVector[2] + acVector[3];//numA + numG + numT;
        fnov[2] = acVector[0] + acVector[1] + acVector[3];//numA + numC + numT;
        fnov[3] = acVector[0] + acVector[1] + acVector[2];//numA + numC + numG;

        Double[] p1 = MathUtils.vectorLog10(MathUtils.vectorSum(MathUtils.scalarTimesVector(t1, acVector), MathUtils.scalarTimesVector(t2, fnov)));

        return errorModel.getErrorGivenQual(qual) + MathUtils.dotProduct(p1,numSeenBases);
    }

    /**
     * Performs a quick check to see if the maximum likelihood is at AC == 0 and that
     * it is above the calling threshold.
     *
     * @return true if the allele count is zero, false if AC > 0 or it is not confident enough to make any calls.
     */
    public boolean isACZero() {
        return getMaximumLikelihoodIndex() == 0 && isConfidentlyCalled();
    }

    public int getMaximumLikelihoodIndex()  {
        return idxOfMaxAC;
    }

    /**
     * Determines whether or not this model has a maximum likelihood that passes the calling threshold
     *
     * @return true if the most likely AC passes the calling threshold.
     */
    public boolean isConfidentlyCalled() {
        return phredScaledConfidence > minCallQuall;
    }

    /**
     * The absolute value of the log10 likelihood of the call.
     *   For AC==0 this is the likelihood that the call is AC=0.
     *   For AC>0 this is the likelihood that the call is not AC=0.
     *
     * @return The absolute value of the log10 likelihood of the call.
     */
    public double getNegLog10PError() {
 /*       double result;
        if (isACZero()) {
            result = Math.abs(model[0]);
        }
        else {
            double [] normalizedModel = MathUtils.normalizeFromLog10(model);
            result = Math.abs(Math.log10(1 - normalizedModel[0]));
        }            */
        return phredScaledConfidence;
    }

    /**
     * Only call if site quality > Q(1/2N) with 99% confidence (point where cumsum(site) == log10(0.99))
     * @return
     */
    public boolean isErrorModelPowerfulEnough() {
        return errorModel.hasPowerForMaxAC(maxAlleleCount);
    }

    public byte getAltBase() {
        return bestAltBase;
    }

    public boolean isVariant() {
        return (refBase != bestAltBase);
    }

}

