package org.broadinstitute.sting.gatk.walkers.poolcaller;

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

    public Lane(LaneParameters p) {
        name = p.name;
        referenceSample = new ReferenceSample(new ReferenceSampleParameters(p.referenceSampleName, p.lanePileup.getPileupForSampleName(p.referenceSampleName), p.trueReferenceBases));
        errorModel = new ErrorModel(new ErrorModelParameters(p.minQualityScore, p.maxQualityScore, p.phredScaledPrior, referenceSample, p.minPower));

        Collection<String> poolNames = p.lanePileup.getSampleNames();
        poolNames.remove(p.referenceSampleName);
        pools = new LinkedList<Pool>();
        for (String poolName : poolNames) {
            pools.add(new Pool(new PoolParameters(poolName,
                                                  p.lanePileup.getPileupForSampleName(poolName),
                                                  errorModel,
                                                  p.referenceSequenceBase,
                                                  p.maxAlleleCount,
                                                  p.minCallQual)));
        }

        filters = new TreeSet<String>();
        for (Pool pool : pools) {
            // make the first pool's alleleCountModel our base model and then "recursively" merge it with all the subsequent pools
            if (alleleCountModel == null)
                alleleCountModel = new AlleleCountModel(pool.getAlleleCountModel());
            else
                alleleCountModel.merge(pool.getAlleleCountModel());

            // add all filters from this pool
            filters.addAll(pool.getFilters());
            attributes.putAll(pool.getAttributes());
        }
    }

    private Lane() {}


    public static Lane debugLane(LaneParameters p) {
        Lane lane = new Lane();
        lane.name = p.name;
        lane.referenceSample = new ReferenceSample(new ReferenceSampleParameters(p.referenceSampleName, p.lanePileup.getPileupForSampleName(p.referenceSampleName), p.trueReferenceBases));
        lane.errorModel = new ErrorModel(new ErrorModelParameters(p.minQualityScore, p.maxQualityScore, p.phredScaledPrior, lane.referenceSample, p.minPower));

        lane.pools = new LinkedList<Pool>();
        lane.pools.add(new Pool(new PoolParameters("POOL1", p.lanePileup, lane.errorModel, p.referenceSequenceBase, p.maxAlleleCount, p.minCallQual)));

        lane.filters = new TreeSet<String>();
        for (Pool pool : lane.pools) {
            lane.alleleCountModel = new AlleleCountModel(pool.getAlleleCountModel());
            lane.filters.addAll(pool.getFilters());
            lane.attributes.putAll(pool.getAttributes());
        }
        return lane;
    }

    public Set<String> getFilters() {
        return filters;
    }

    public Map<String, ?> getAttributes() {
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
