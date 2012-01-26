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

    private List<Pool> pools;
    private AlleleCountModel alleleCountModel;
    private boolean isVariant;

    public Lane(String name, ReadBackedPileup lanePileup, String referenceSampleName, ErrorModel errorModel, byte referenceSequenceBase,
                int maxAlleleCount, double minCallQual, int minRefDepth,
                boolean doAlleleDiscovery, boolean treatAllSamplesAsSinglePool) {
        this.name = name;


        Collection<String> poolNames = lanePileup.getSamples();
        poolNames.remove(referenceSampleName);
        this.pools = new LinkedList<Pool>();

        if(treatAllSamplesAsSinglePool) {
            pools.add(new Pool("Pool1", lanePileup.getPileupForSamples(poolNames),errorModel,referenceSequenceBase,maxAlleleCount,minCallQual, minRefDepth, doAlleleDiscovery));
        }
        else {
            // regular case: pileup is stratified by pool names.
            // Get then all samples in list except for reference sample
            for (String poolName : poolNames)
                pools.add(new Pool(poolName, lanePileup.getPileupForSample(poolName),errorModel,referenceSequenceBase,maxAlleleCount,minCallQual, minRefDepth, doAlleleDiscovery));
        }


        this.isVariant = false;

        for (Pool pool : pools) {
            // make the first pool's alleleCountModel our base model and then "recursively" merge it with all the subsequent pools
            if (alleleCountModel == null)
                this.alleleCountModel = pool.getAlleleCountModel();
            else
                this.alleleCountModel.merge(pool.getAlleleCountModel());

             isVariant |= pool.isVariant();
        }
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
