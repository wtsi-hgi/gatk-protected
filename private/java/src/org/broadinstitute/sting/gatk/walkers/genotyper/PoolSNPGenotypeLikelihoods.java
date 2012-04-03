package org.broadinstitute.sting.gatk.walkers.genotyper;


import net.sf.samtools.SAMUtils;
import org.broadinstitute.sting.gatk.walkers.poolcaller.ErrorModel;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;

import java.util.*;

import static java.lang.Math.log;
import static java.lang.Math.log10;

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
public class PoolSNPGenotypeLikelihoods extends PoolGenotypeLikelihoods/* implements Cloneable*/ {

    // TODO: don't calculate this each time through
    //protected double log10_PCR_error_3;
    //protected double log10_1_minus_PCR_error;


    //    double[][][] logPLs;
    private double[][] logMismatchProbabilityArray;
    static private final double qualVec[] = new double[SAMUtils.MAX_PHRED_SCORE+1];

    /**
     * Create a new GenotypeLikelhoods object with given priors and PCR error rate for each pool genotype
     *
     * @param priors          priors
     * @param perLaneErrorModels error model objects for each lane
     * @param ignoreLaneInformation  If true, lane info is ignored
     */
    public PoolSNPGenotypeLikelihoods(List<Allele> alleles, PoolGenotypePriors priors, HashMap<String, ErrorModel> perLaneErrorModels, boolean ignoreLaneInformation) {

        super(alleles,priors, perLaneErrorModels, ignoreLaneInformation);
        fillCache();
    }
    public PoolSNPGenotypeLikelihoods(final List<Allele> alleles, final double[] logLikelihoods, final int ploidy,
                                      final HashMap<String, ErrorModel> perLaneErrorModels, final boolean ignoreLaneInformation) {
        super(alleles, logLikelihoods, ploidy, perLaneErrorModels, ignoreLaneInformation);
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
            if (ignoreLaneInformation)
                perLanePileup = pileup;
            else
                perLanePileup = pileup.getPileupForLane(laneID);

            if (perLanePileup == null || perLanePileup.isEmpty())
                continue;

            ErrorModel errorModel = perLaneErrorModels.get(laneID);
            n += add(perLanePileup, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual, errorModel);
            if (ignoreLaneInformation)
                break;

        }
//        log10Likelihoods = convertToVectorIndex(logPLs, numChromosomes);
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
        if (VERBOSE)
            System.out.format("numSeenBases: %d %d %d %d\n",numSeenBases[0],numSeenBases[1],numSeenBases[2],numSeenBases[3]);

        add(numSeenBases[0], numSeenBases[1], numSeenBases[2], numSeenBases[3], errorModel);
        return n;
    }

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
    private int add(final int nA, final int nC, final int nG, final int nT, final ErrorModel errorModel) {
        final int minQ = errorModel.getMinQualityScore();
        final int maxQ = errorModel.getMaxQualityScore();

        List<Allele> myAlleles = new ArrayList<Allele>(alleles);

        final int[] initialState = getInitialStateVector(nAlleles,numChromosomes);
        
        int[] newInitialState;
        
        if (myAlleles.size() < BaseUtils.BASES.length) {
             newInitialState = new int[BaseUtils.BASES.length];
            // likelihood only defined for subset of possible alleles. Fill then with other alleles to have all possible ones,
            for (byte b : BaseUtils.BASES) {
                // if base is not included in myAlleles, add new allele
                boolean found = false;
                for (Allele a: alleles) {
                    if (a.getBases()[0] == b) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    myAlleles.add(Allele.create(b,false));
                }


            }

            int k=0;
            for (int is: initialState)
                newInitialState[k++] = is;
        }
        else
            newInitialState = initialState;

        SumIterator iterator = new SumIterator(newInitialState,numChromosomes);

        int kk=0;
        // compute permutation vector to figure out mapping from indices to bases
        int idx = 0;
        int[] alleleIndices = new int[myAlleles.size()];
        for (byte b : BaseUtils.BASES) {
            kk=0;
            for (Allele a: myAlleles) {
                if (a.getBases()[0] == b) {
                    alleleIndices[idx++] = kk;
                    break;
                }
                kk++;
            }
        }

        while (iterator.hasNext()) {
            // for observed base X, add Q(jX,k) to likelihood vector for all k in error model
            //likelihood(jA,jC,jG,jT) = logsum(logPr (errorModel[k],nA*Q(jA,k) +  nC*Q(jC,k) + nG*Q(jG,k) + nT*Q(jT,k))
            double[] acVec = new double[maxQ - minQ + 1];
            int[] currentCnt = iterator.getCurrentVector();
            int jA = currentCnt[alleleIndices[0]];
            int jC = currentCnt[alleleIndices[1]];
            int jG = currentCnt[alleleIndices[2]];
            int jT = currentCnt[alleleIndices[3]];
            for (int k=minQ; k<=maxQ; k++)
                acVec[k-minQ] = nA*logMismatchProbabilityArray[jA][k] +
                        nC*logMismatchProbabilityArray[jC][k] +
                        nG*logMismatchProbabilityArray[jG][k] +
                        nT*logMismatchProbabilityArray[jT][k];

            double pl = MathUtils.logDotProduct(errorModel.getErrorModelVector(), acVec);
            if (alleles.size() < BaseUtils.BASES.length)
                currentCnt = Arrays.copyOfRange(currentCnt,0,alleles.size());
            setLogPLs(currentCnt, pl);
            kk++;
            iterator.next();
        }

//System.out.format("Put k=%d\n",kk);

        return 1;
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


    static {
        // cache 10^(-k/10)
        for (int j=0; j <= SAMUtils.MAX_PHRED_SCORE; j++)
            qualVec[j] = Math.pow(10.0,-(double)j/10.0);
    }

    private void fillCache() {
        // cache Q(j,k) = log10(j/2N*(1-ek) + (2N-j)/2N*ek) for j = 0:2N

        logMismatchProbabilityArray = new double[1+numChromosomes][1+SAMUtils.MAX_PHRED_SCORE];
        for (int i=0; i <= numChromosomes; i++) {
            for (int j=0; j <= SAMUtils.MAX_PHRED_SCORE; j++) {
                double phi = (double)i/numChromosomes;
                logMismatchProbabilityArray[i][j] = Math.log10(phi * (1.0-qualVec[j]) + qualVec[j]/3.0 * (1.0-phi));
            }
        }
    }
}
