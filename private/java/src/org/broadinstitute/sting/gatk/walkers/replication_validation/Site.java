package org.broadinstitute.sting.gatk.walkers.replication_validation;

import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 7/21/11
 * Time: 6:13 PM
 *
 * This is the concept implementation of a site in the replication/validation framework.
 */
public class Site {
    private Set<Lane> lanes;
    private Set<String> laneIDs;
    private AlleleCountModel alleleCountModel;

    public Site(ReadBackedPileup sitePileup, String referenceSampleName, Collection<Byte> trueReferenceBases, byte referenceSequenceBase, byte minQualityScore, byte maxQualityScore, byte phredScaledPrior, int maxAlelleCount) {
        laneIDs = parseLaneIDs(sitePileup.getReadGroups());
        lanes = new HashSet<Lane>(laneIDs.size());
        for (String laneID : laneIDs) {
            lanes.add(new Lane(laneID, sitePileup.getPileupForLane(laneID), referenceSampleName, trueReferenceBases, referenceSequenceBase, minQualityScore, maxQualityScore, phredScaledPrior, maxAlelleCount));
        }

        for (Lane lane : lanes) {
            // make the first pool's alleleCountModel our base model and then "recursively" merge it with all the subsequent pools
            if (alleleCountModel == null)
                alleleCountModel = new AlleleCountModel(lane.getAlleleCountModel());
            else
                alleleCountModel.merge(lane.getAlleleCountModel());
        }
    }

    public Set<Lane> getLanes() {
        return lanes;
    }

    public Set<String> getLaneIDs() {
        return laneIDs;
    }

    public AlleleCountModel getAlleleCountModel() {
        return alleleCountModel;
    }

    /**
     * GATK Engine creates readgroups of the form XXX.Y.Z
     * XXX.Y is the unique lane identifier.
     *     Z is the id of the sample to make the read group id unique
     * This function returns the list of lane identifiers.
     *
     * @param readGroups readGroups A collection of read group strings (obtained from the alignment context pileup)
     * @return a collection of lane ids.
     */
    private Set<String> parseLaneIDs(Collection<String> readGroups) {
        HashSet<String> result = new HashSet<String>();
        for (String readGroup : readGroups) {
            result.add(getLaneIDFromReadGroupString(readGroup));
        }
        return result;
    }

    /**
     * GATK Engine creates readgroups of the form XXX.Y.Z
     * XXX.Y is the unique lane identifier.
     *     Z is the id of the sample to make the read group id unique
     *
     * @param readGroupID the read group id string
     * @return just the lane id (the XXX.Y string)
     */
    private String getLaneIDFromReadGroupString(String readGroupID) {
        String [] parsedID = readGroupID.split("\\.");
        return parsedID[0] + "." + parsedID[1];
    }
}
