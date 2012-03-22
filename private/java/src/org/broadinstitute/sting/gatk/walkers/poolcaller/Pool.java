package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 7/21/11
 * Time: 1:29 PM
 *
 * This is a site based implementation of a pool.
 * A pool object will contain all the information pertinent to a pool in one given site
 */

public class Pool {
    private String name;
    //private ReadBackedPileup pileup;
    private AlleleCountModel alleleCountModel;
    //private int matches;
    //private int mismatches;
    //private byte referenceSequenceBase;
    //private Set<Filters> filters;
    private int coverage;
    private Integer calledAC;
    private Allele calledAllele;
    private Allele refAllele;
    private Genotype poolGenotype;
    private boolean isVariant;
   // private double log10LikelihoodCall;
    
    public static final String MAXIMUM_LIKELIHOOD_AC_KEY = "MLAC";
    public static final String NINETY_FIVE_PCT_CONFIDENCE_INTERVAL_KEY = "AC95";

    

    public Pool(String name, ReadBackedPileup pileup, ErrorModel errorModel, int maxAlleleCount, double minCallQual,
                Allele refAllele, boolean doAlleleDiscovery, List<Allele> allelesToTest) {
        this.name = name;
 //       this.pileup = pileup;
        this.refAllele = refAllele;
//        this.maxAlleleCount = maxAlleleCount;
//        this.referenceSequenceBase = referenceSequenceBase;

        byte [] data = new byte[0];
        if (pileup != null)
            data = pileup.getBases();
        //else


        coverage = data.length;

        int idx = 0;
        Integer[] numSeenBases = new Integer[BaseUtils.BASES.length];
        for (byte base:BaseUtils.BASES)
            numSeenBases[idx++] = MathUtils.countOccurrences(base, pileup.getBases());

        // todo - generalize for indels
     //   System.out.format("A:%d C:%d G:%d T:%t\n",numSeenBases[0],numSeenBases[1],numSeenBases[2],numSeenBases[3]);
        alleleCountModel = new AlleleCountModel(maxAlleleCount, errorModel, numSeenBases, minCallQual, refAllele, doAlleleDiscovery, allelesToTest);

       // isConfidentlyCalled = alleleCountModel.isConfidentlyCalled();


       fillGenotypeInformation();


  //      log10LikelihoodCall = alleleCountModel.getMaximumLikelihood();
    }

    /**
     * @return the name of the pool (sample name of the pool in the bam file)
     */
    public String getName() {
        return name;
    }

    public int getCoverage() {
        return coverage;
    }

    public Pool mergePool(Pool other) {
        alleleCountModel.merge(other.getAlleleCountModel());
        fillGenotypeInformation();
        return this;
    }

    private void fillGenotypeInformation() {
        calledAC = alleleCountModel.getMaximumLikelihoodIndex();
        if (calledAC == 0) {
            calledAllele = refAllele;
            isVariant = false;
        }
        else {
            calledAllele = alleleCountModel.getAltAllele();
            isVariant = true;
        }

        List<Allele> alleleList = new LinkedList<Allele>();
        alleleList.add(calledAllele);

        Map<String, Object> attributes = new HashMap<String, Object>();
        attributes.put(VCFConstants.DEPTH_KEY, coverage);
        attributes.put(MAXIMUM_LIKELIHOOD_AC_KEY, calledAC);
        Pair<Integer,Integer> P =  alleleCountModel.get95PctACConfidenceInterval();
        attributes.put(NINETY_FIVE_PCT_CONFIDENCE_INTERVAL_KEY, String.format("%d,%d",P.first,P.second));

        poolGenotype = new Genotype(name, alleleList, -getAlleleCountModel().getNegLog10PError(),null,attributes,false);
    }


    /**
     * @return the allele count model
     */
    public AlleleCountModel getAlleleCountModel() {
        return alleleCountModel;
    }

    /**
     * @return Whether or not the site is filtered
     */
    public boolean isFiltered() {
        return !isCalled();
    }

    /**
     * @return whether or not the site is called
     */
    public boolean isCalled() {
        return calledAC != null;
    }


    public boolean isVariant() {
        return alleleCountModel.isVariant();
    }
    /**
     * Builds the Genotype object for the pool. It takes the most frequent alternate allele if the
     * pool is variant or the ref allele if it isn't.
     *
     * @return the Genotype object of the pool.
     */
    public Genotype getGenotype() {
        return poolGenotype;
    }

}
