package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.utils.pileup.PileupElement;
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
    private boolean isVariant;

    public Lane(String name, ReadBackedPileup lanePileup, String referenceSampleName, Collection<Byte> trueReferenceBases, byte referenceSequenceBase,
                byte minQualityScore, byte maxQualityScore, byte phredScaledPrior,int maxAlleleCount, double minCallQual, double minPower,  int minRefDepth, boolean doAlleleDiscovery) {
        this.name = name;
        this.referenceSample = new ReferenceSample(referenceSampleName, lanePileup.getPileupForSample(referenceSampleName), trueReferenceBases);
        this.errorModel = new ErrorModel(minQualityScore, maxQualityScore, phredScaledPrior, referenceSample, minPower);

        Collection<String> poolNames = lanePileup.getSamples();
        poolNames.remove(referenceSampleName);
        this.pools = new LinkedList<Pool>();
        for (String poolName : poolNames) {
            pools.add(new Pool(poolName, lanePileup.getPileupForSample(poolName),errorModel,referenceSequenceBase,maxAlleleCount,minCallQual, minRefDepth, doAlleleDiscovery));
        }

        this.filters = new TreeSet<String>();
        this.isVariant = false;

        for (Pool pool : pools) {
            // make the first pool's alleleCountModel our base model and then "recursively" merge it with all the subsequent pools
            if (alleleCountModel == null)
                this.alleleCountModel = pool.getAlleleCountModel();
            else
                this.alleleCountModel.merge(pool.getAlleleCountModel());

            // add all filters from this pool
            filters.addAll(pool.getFilters());
            attributes.putAll(pool.getAttributes());
            isVariant |= pool.isVariant();
        }
        // add reference sample depth
        attributes.put("RD",errorModel.getReferenceDepth());
    }

    private Lane() {}


    public static Lane debugLane(String name, ReadBackedPileup lanePileup, String referenceSampleName, Collection<Byte> trueReferenceBases,
                                 byte referenceSequenceBase, byte minQualityScore, byte maxQualityScore, byte phredScaledPrior, int maxAlleleCount, double minCallQual, double minPower, int minRefDepth, boolean doAlleleDiscovery) {
        Lane lane = new Lane();
        lane.name = name;
        lane.referenceSample = new ReferenceSample(referenceSampleName, lanePileup.getPileupForSample(referenceSampleName), trueReferenceBases);
        lane.errorModel = new ErrorModel(minQualityScore, maxQualityScore, phredScaledPrior, lane.referenceSample, minPower);

        Collection<String> allSamples = new HashSet<String>();
        for (String sample: lanePileup.getSamples()) {
            if (sample.compareToIgnoreCase(referenceSampleName)!=0)
                allSamples.add(sample);

        }
        ReadBackedPileup samplePileup = lanePileup.getPileupForSamples(allSamples);
        lane.pools = new LinkedList<Pool>();
        lane.pools.add(new Pool("POOL1", samplePileup, lane.errorModel, referenceSequenceBase, maxAlleleCount, minCallQual, minRefDepth, doAlleleDiscovery));

        for (Pool pool : lane.pools) {
            lane.alleleCountModel = pool.getAlleleCountModel();
            lane.filters = pool.getFilters();
            lane.attributes = pool.getAttributes();
            lane.isVariant |= pool.isVariant();
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

    public boolean isVariant() {
        return isVariant;
    }
}
