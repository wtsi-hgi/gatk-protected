package org.broadinstitute.sting.gatk.walkers.replication_validation;

import com.google.java.contract.Requires;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.AlleleCount;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

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


    public Lane(String name, ReadBackedPileup lanePileup, String referenceSampleName, Collection<Byte> trueReferenceBases, byte referenceSequenceBase, byte minQualityScore, byte maxQualityScore, byte phredScaledPrior, int maxAlleleCount) {
        this.name = name;
        referenceSample = new ReferenceSample(referenceSampleName, lanePileup.getPileupForSampleName(referenceSampleName), trueReferenceBases);
        errorModel = new ErrorModel(minQualityScore, maxQualityScore, phredScaledPrior, referenceSample);

        Collection<String> poolNames = lanePileup.getSampleNames();
        poolNames.remove(referenceSampleName);
        pools = new LinkedList<Pool>();
        for (String poolName : poolNames) {
            pools.add(new Pool(poolName, lanePileup.getPileupForSampleName(poolName), errorModel, referenceSequenceBase, maxAlleleCount));
        }

        for (Pool pool : pools) {
            // make the first pool's alleleCountModel our base model and then "recursively" merge it with all the subsequent pools
            if (alleleCountModel == null)
                alleleCountModel = new AlleleCountModel(pool.getACModel());
            else
                alleleCountModel.merge(pool.getACModel());
        }
    }

    public List<Pool> getPools() {
        return pools;
    }

    public ErrorModel getErrorModel() {
        return errorModel;
    }

    public AlleleCountModel getAlleleCountModel() {
        return alleleCountModel;
    }
}
