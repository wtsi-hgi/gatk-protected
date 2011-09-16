package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 7/21/11
 * Time: 4:33 PM
 *
 * This is the concept implementation of a site in the replication/validation framework. A site must contain
 * a reference sample with an error model and a list of pools.
 */
public class Lane {
    private String name;
    private ReferenceSample referenceSample;
    private List<Pool> pools;
    private ErrorModel errorModel;
    private AlleleCountModel alleleCountModel;
    private Set<String> filters;
    private Map<String, Object> attributes;

    public Lane(String name, ReadBackedPileup lanePileup, String referenceSampleName, Collection<Byte> trueReferenceBases, byte referenceSequenceBase, byte minQualityScore, byte maxQualityScore, byte phredScaledPrior,int maxAlleleCount, double minCallQual, double minPower) {
        this.name = name;
        this.referenceSample = new ReferenceSample(referenceSampleName, lanePileup.getPileupForSampleName(referenceSampleName), trueReferenceBases);
        this.errorModel = new ErrorModel(minQualityScore, maxQualityScore, phredScaledPrior, referenceSample, minPower);

        Collection<String> poolNames = lanePileup.getSampleNames();
        poolNames.remove(referenceSampleName);
        this.pools = new LinkedList<Pool>();
        for (String poolName : poolNames) {
            pools.add(new Pool(poolName, lanePileup.getPileupForSampleName(poolName),errorModel,referenceSequenceBase,maxAlleleCount,minCallQual));
        }

        this.filters = new TreeSet<String>();
        for (Pool pool : pools) {
            // make the first pool's alleleCountModel our base model and then "recursively" merge it with all the subsequent pools
            if (alleleCountModel == null)
                this.alleleCountModel = new AlleleCountModel(pool.getAlleleCountModel());
            else
                this.alleleCountModel.merge(pool.getAlleleCountModel());

            // add all filters from this pool
            filters.addAll(pool.getFilters());
            attributes.putAll(pool.getAttributes());
        }
    }

    private Lane() {}


    public static Lane debugLane(String name, ReadBackedPileup lanePileup, String referenceSampleName, Collection<Byte> trueReferenceBases, byte referenceSequenceBase, byte minQualityScore, byte maxQualityScore, byte phredScaledPrior, int maxAlleleCount, double minCallQual, double minPower) {
        Lane lane = new Lane();
        lane.name = name;
        lane.referenceSample = new ReferenceSample(referenceSampleName, lanePileup.getPileupForSampleName(referenceSampleName), trueReferenceBases);
        lane.errorModel = new ErrorModel(minQualityScore, maxQualityScore, phredScaledPrior, lane.referenceSample, minPower);

        lane.pools = new LinkedList<Pool>();
        lane.pools.add(new Pool("POOL1", lanePileup, lane.errorModel, referenceSequenceBase, maxAlleleCount, minCallQual));

        for (Pool pool : lane.pools) {
            lane.alleleCountModel = new AlleleCountModel(pool.getAlleleCountModel());
            lane.filters = pool.getFilters();
            lane.attributes = pool.getAttributes();
        }
        return lane;
    }

    public Set<String> getFilters() {
        return filters;
    }

    public Map<String, Object> getAttributes() {
        return attributes;
    }

    public AlleleCountModel getAlleleCountModel() {
        return alleleCountModel;
    }

    public Map<String, Genotype> getGenotypes() {
        Map<String,Genotype> laneMap = new HashMap<String, Genotype>();
        for (Pool pool : pools) {
            laneMap.putAll(pool.getGenotypes());
        }
        return laneMap;
    }
}
