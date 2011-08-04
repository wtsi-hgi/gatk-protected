package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

import java.util.*;
import java.util.logging.Filter;

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
    private ReadBackedPileup pileup;
    private int maxAlleleCount;
    private AlleleCountModel alleleCountModel;
    private int matches;
    private int mismatches;
    private byte referenceSequenceBase;
    private Set<Filters> filters;
    private Map<String, Object> attributes;
    private boolean isConfidentlyCalled;
    private Integer calledAC;
    private byte calledAllele;
    private double log10LikelihoodCall;

    public Pool(PoolParameters p) {
        name = p.name;
        pileup = p.pileup;
        maxAlleleCount = p.maxAlleleCount;
        referenceSequenceBase = p.referenceSequenceBase;

        byte [] data = pileup.getBases();
        int coverage = data.length;
        matches = MathUtils.countOccurrences(p.referenceSequenceBase, data);
        mismatches = coverage - matches;

        alleleCountModel = new AlleleCountModel(new AlleleCountModelParameters(maxAlleleCount, p.errorModel, matches, mismatches, p.minCallQual));

        // make the call and apply filters
        filters = new HashSet<Filters>();
        isConfidentlyCalled = alleleCountModel.isConfidentlyCalled();
        if (!alleleCountModel.isConfidentlyCalled())
            filters.add(Filters.LOW_QUAL);
        else if (!alleleCountModel.isErrorModelPowerfulEnough())
            filters.add(Filters.LOW_POWER);
        else {
            calledAC = alleleCountModel.getMaximumLikelihoodIndex();
            calledAllele = (calledAC == 0) ?  referenceSequenceBase : BaseUtils.baseIndexToSimpleBase(BaseUtils.mostFrequentBaseIndexNotRef(pileup.getBaseCounts(), referenceSequenceBase));
        }
        log10LikelihoodCall = alleleCountModel.getMaximumLikelihood();

        calculateAttributes();
    }

    /**
     * @return the name of the pool (sample name of the pool in the bam file)
     */
    public String getName() {
        return name;
    }

    /**
     * @param name The name of the pool (sample name of the pool in the bam file)
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return The pileup for the pool in its location
     */
    public ReadBackedPileup getPileup() {
        return pileup;
    }

    /**
     * @return the number of bases in the pool pileup
     */
    public int size() {
        return pileup.size();
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

    /**
     *
     */
    public boolean isRef() {
        return BaseUtils.baseIndexToSimpleBase(MathUtils.maxElementIndex(pileup.getBaseCounts())) == referenceSequenceBase;
    }

    /**
     * Returns a list of the filters applied to this call. Empty list if nothing was filtered.
     */
    public Set<String> getFilters() {
        Set<String> result = new HashSet<String>(filters.size());
        for (Filters f : filters) {
            result.add(f.toString());    // maybe have a place to get the string definition of each filter in the future?
        }
        return result;
    }

    public Map<String, Object> getAttributes() {
        return attributes;
    }

    private double calculateAlleleFrequency() {
        byte [] bases = pileup.getBases();
        return (double) MathUtils.countOccurrences(calledAllele, bases) / bases.length;
    }

    private int calculateAlleleDepth (byte allele) {
        return MathUtils.countOccurrences(allele, pileup.getBases());
    }

    private List<Integer> calculateAllelicDepths() {
        List<Integer> allelicDepths = new LinkedList<Integer>();
        allelicDepths.add(calculateAlleleDepth(calledAllele));
        allelicDepths.add(calculateAlleleDepth(referenceSequenceBase));
        return allelicDepths;
    }

    private double calculateMappingQualityRMS() {
        return MathUtils.rms(pileup.getMappingQuals());
    }

    private int calculateMappingQualityZero() {
        return pileup.getNumberOfMappingQualityZeroReads();
    }

    private void calculateAttributes() {
        attributes = new HashMap<String, Object>(11);
        attributes.put("AC", calledAC);
        attributes.put("AF", calculateAlleleFrequency());
        attributes.put("DP", pileup.getBases().length);
        attributes.put("AD", calculateAllelicDepths());
        attributes.put("MQ", calculateMappingQualityRMS());
        attributes.put("MQ0", calculateMappingQualityZero());

/*
        headerLines.add(new VCFInfoHeaderLine("AC", 1, VCFHeaderLineType.Integer, "Allele count in the site, number of alternate alleles across all pools"));
        headerLines.add(new VCFInfoHeaderLine("AF", 1, VCFHeaderLineType.Float, "Allele frequency in the site. Proportion of the alternate alleles across all pools"));
        headerLines.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Total depth in the site. Sum of the depth of all pools"));
        headerLines.add(new VCFInfoHeaderLine("MQ", 1, VCFHeaderLineType.Float, "RMS mapping quality of all reads in the site"));
        headerLines.add(new VCFInfoHeaderLine("MQ0", 1, VCFHeaderLineType.Integer, "Total number of mapping quality zero reads in the site"));
        headerLines.add(new VCFFormatHeaderLine("AD", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Allelic depths for the ref and alt alleles in the order listed"));
        headerLines.add(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Read Depth (only filtered reads used for calling)"));
        headerLines.add(new VCFFormatHeaderLine("GQ", 1, VCFHeaderLineType.Float, "Genotype Quality"));
        headerLines.add(new VCFFormatHeaderLine("AL", 3, VCFHeaderLineType.Integer, "Allele count likelihood and the 5% confidence interval"));

*/
    }

    /**
     * Builds the Genotype object for the pool. It takes the most frequent alternate allele if the
     * pool is variant or the ref allele if it isn't.
     *
     * @return the Genotype object of the pool.
     */
    public Map<String, Genotype> getGenotypes() {
        Map<String, Genotype> poolGenotype = new HashMap<String, Genotype>(1);
        byte base;

        // I don't need the reference base to be counted as we are looking for the most common alternate allele.
        if (isRef()) {
            base = referenceSequenceBase;
        }
        else {
            int [] baseCounts = pileup.getBaseCounts();
            baseCounts[BaseUtils.simpleBaseToBaseIndex(referenceSequenceBase)] = -1;
            base = BaseUtils.mostFrequentSimpleBase(baseCounts);
        }

        List<Allele> alleleList = new LinkedList<Allele>();
        alleleList.add(Allele.create(base, isRef()));
        Genotype g = new Genotype(name, alleleList, -log10LikelihoodCall);
        poolGenotype.put(name, g);
        return poolGenotype;
    }

}
