package org.broadinstitute.sting.gatk.walkers.poolcaller;


import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

import java.util.*;

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
    private double minCallQual;
    private byte referenceSequenceBase;
    private ReadBackedPileup pileup;
    private Set<String> filters;
    private Map<String, Object> attributes;



    /**
     * General constructor for a Site object. It creates a site with all the lanes represented in it's pileup,
     * the error model based on the reference sample pileup and calculates the allele count model for the entire
     * lane.
     */

    public Site(SiteParameters p) {
        this.referenceSequenceBase = p.referenceSequenceBase;
        this.minCallQual = p.minCallQual;
        this.pileup = p.sitePileup;

        laneIDs = parseLaneIDs(p.sitePileup.getReadGroups());
        lanes = new HashSet<Lane>(laneIDs.size());
        for (String laneID : laneIDs) {
            lanes.add(new Lane(new LaneParameters(laneID,
                                                  p.sitePileup.getPileupForLane(laneID),
                                                  p.referenceSampleName,
                                                  p.trueReferenceBases,
                                                  p.referenceSequenceBase,
                                                  p.minQualityScore,
                                                  p.maxQualityScore,
                                                  p.phredScaledPrior,
                                                  p.maxAlelleCount,
                                                  p.minCallQual,
                                                  p.minPower)));
        }

        for (Lane lane : lanes) {
            // make the first pool's alleleCountModel our base model and then "recursively" merge it with all the subsequent pools
            if (alleleCountModel == null)
                alleleCountModel = new AlleleCountModel(lane.getAlleleCountModel());
            else
                alleleCountModel.merge(lane.getAlleleCountModel());

            // Add all the filters from this lane
            filters.addAll(lane.getFilters());
            attributes.putAll(lane.getAttributes());
        }
    }


    /**
     * Special case constructor when site consists of only one pool *use debugSite instead*
     * @param laneID the lane id
     * @param lane the lane object
     */
    private Site(String laneID, Lane lane, double minCallQual, byte referenceSequenceBase, ReadBackedPileup sitePileup) {
        this.lanes = new HashSet<Lane>(1);
        this.laneIDs = new HashSet<String>(1);
        this.lanes.add(lane);
        this.laneIDs.add(laneID);
        this.pileup = sitePileup;
        this.minCallQual = minCallQual;
        this.referenceSequenceBase = referenceSequenceBase;
        alleleCountModel = new AlleleCountModel(lane.getAlleleCountModel());
        filters = lane.getFilters();
        attributes = lane.getAttributes();
    }

    /**
     * Same as Site constructor but will create only one lane containing only one pool with all samples in the SAM record.
     * Only for debug purposes.
     */
    public static Site debugSite(SiteParameters p) {
        String laneID = "LANE1";
        Lane lane = Lane.debugLane(new LaneParameters(laneID,
                                                p.sitePileup,
                                                p.referenceSampleName,
                                                p.trueReferenceBases,
                                                p.referenceSequenceBase,
                                                p.minQualityScore,
                                                p.maxQualityScore,
                                                p.phredScaledPrior,
                                                p.maxAlelleCount,
                                                p.minCallQual,
                                                p.minPower));

        return new Site(laneID, lane, p.minCallQual, p.referenceSequenceBase, p.sitePileup);
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

    public double getNegLog10PError() {
        return alleleCountModel.getNegLog10PError();
    }

    public Collection<Allele> getAlleles() {
        Set<Allele> alternateAlleles = new HashSet<Allele>();
        alternateAlleles.add(Allele.create(referenceSequenceBase, true));
        if (alleleCountModel.isConfidentlyCalled()) {
            int [] baseCounts = pileup.getBaseCounts();
            // if reference base is the most common, make it negative so we get the second most common.
            baseCounts[BaseUtils.simpleBaseToBaseIndex(referenceSequenceBase)] = -1;

            // add the most common alternate allele
            alternateAlleles.add(Allele.create(BaseUtils.baseIndexToSimpleBase(MathUtils.maxElementIndex(baseCounts)), false));
        }

        else {
            // todo -- not sure what to do with not confident calls
            alternateAlleles.add(Allele.create(referenceSequenceBase, true)); // adding reference allele because I don't know how to handle this case.
        }
        return alternateAlleles;
    }

    public Map<String, Genotype> getGenotypes() {
        Map<String, Genotype> siteGenotypeMap = new HashMap<String, Genotype>();
        for (Lane lane : lanes) {
            Map<String, Genotype> laneGenotypeMap = lane.getGenotypes();
            siteGenotypeMap.putAll(laneGenotypeMap);
        }
        return siteGenotypeMap;
    }

    public Set<String> getFilters() {
        return filters;
    }

    public Map<String, Object> getAttributes () {
        return attributes;
    }



}
