package org.broadinstitute.sting.gatk.walkers.poolcaller;


import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.*;

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
    private int minRefDepth;
    private ErrorModel errorModel;
    private ReadBackedPileup pileup;
    private Set<String> filters;
    private Map<String, Object> attributes;
    private final boolean filteredRefSampleCall;
    private ReferenceContext refContext;
    private UnifiedArgumentCollection UAC;
    private VariantContext siteVC;


    /**
     * General constructor for a Site object. It creates a site with all the lanes represented in it's pileup,
     * the error model based on the reference sample pileup and calculates the allele count model for the entire
     * lane.
     */

    public Site(ReferenceContext refContext, AlignmentContext context, String referenceSampleName,Collection<Byte> trueReferenceBases,
                byte minQualityScore, byte maxQualityScore,
                byte phredScaledPrior, int maxAlleleCount, UnifiedArgumentCollection UAC, double minPower, int minRefDepth, boolean filteredRefSampleCall,
                boolean treatAllSamplesAsSinglePool) {

        //this.referenceSequenceBase = referenceSequenceBase;
        this.refContext = refContext;
        this.UAC = UAC;
        double minCallQual = UAC.STANDARD_CONFIDENCE_FOR_CALLING;
        this.minRefDepth = minRefDepth;
        this.filteredRefSampleCall = filteredRefSampleCall;

        //this.pileup = sitePileup;
        //this.doDiscoveryMode = doDiscoveryMode;
        ReferenceSample referenceSample = new ReferenceSample(referenceSampleName, context.getBasePileup().getPileupForSample(referenceSampleName), trueReferenceBases);
        this.errorModel = new ErrorModel(minQualityScore, maxQualityScore, phredScaledPrior, referenceSample, minPower);

        // so start: generate an empty VC in case we need to emit at all sites
        Set<Allele> alleles = new HashSet<Allele>();
        alleles.add(Allele.create(refContext.getBase(), true));
        siteVC = new VariantContextBuilder("PoolCaller", refContext.getLocus().getContig(), refContext.getLocus().getStart(), refContext.getLocus().getStart(), alleles).make();

        if (context.hasBasePileup()) {
            filters = new HashSet<String>();
            pileup = context.getBasePileup();
            boolean doDiscoveryMode = (UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY);

            if (treatAllSamplesAsSinglePool) {
                laneIDs = new HashSet<String>(1);
                laneIDs.add("Lane1");
            }
            else
                laneIDs = parseLaneIDs(pileup.getReadGroups());

            lanes = new HashSet<Lane>(laneIDs.size());
            if (treatAllSamplesAsSinglePool) {

                lanes.add(new Lane("Lane1",pileup, referenceSampleName, errorModel, refContext.getBase(),
                        maxAlleleCount,minCallQual, minRefDepth, doDiscoveryMode, treatAllSamplesAsSinglePool));

            }
            else {
                for (String laneID : laneIDs) {
                    lanes.add(new Lane(laneID,pileup.getPileupForLane(laneID),referenceSampleName,errorModel, refContext.getBase(),
                            maxAlleleCount,minCallQual, minRefDepth,
                            doDiscoveryMode, treatAllSamplesAsSinglePool));
                }
                }

            for (Lane lane : lanes) {
                // make the first pool's alleleCountModel our base model and then "recursively" merge it with all the subsequent pools
                if (alleleCountModel == null)
                    alleleCountModel = lane.getAlleleCountModel();
                else
                    alleleCountModel.merge(lane.getAlleleCountModel());

                // Add all the filters from this lane
                // todo - this is wrong, we need to merge filters and attributes from different lanes!
                calculateFilters();
                calculateAttributes();
            }

            siteVC = new VariantContextBuilder("PoolCaller",
                    refContext.getLocus().getContig(),
                    refContext.getLocus().getStart(),
                    refContext.getLocus().getStop(),
                    getAlleles())
                    .genotypes(GenotypesContext.copy(getGenotypes().values()))
                    .log10PError(-1 * getNegLog10PError())
                    .filters(filters)
                    .attributes(attributes).make();
        }
    }

    private void calculateFilters() {
        // make the call and apply filters
//        filters = new HashSet<Filters>();
        if (!alleleCountModel.isConfidentlyCalled())
            filters.add(Filters.LOW_QUAL.toString());
        if (!alleleCountModel.isErrorModelPowerfulEnough())
            filters.add(Filters.LOW_POWER.toString());

        if (!errorModel.hasData())
            filters.add(Filters.NO_BASES_IN_REFERENCE_SAMPLE.toString());
        else if (errorModel.getReferenceDepth() < minRefDepth)
            filters.add(Filters.LOW_REFERENCE_SAMPLE_DEPTH.toString());

        if (filteredRefSampleCall)
            filters.add(Filters.FILTERED_REFERENCE_SAMPLE_CALL.toString());

    }
    private int calculateAlleleDepth (byte allele) {
        return MathUtils.countOccurrences(allele, pileup.getBases());
    }

    private List<Integer> calculateAllelicDepths(Collection<Allele> alleles) {
        List<Integer> allelicDepths = new LinkedList<Integer>();
        
        for (Allele a: alleles)
            allelicDepths.add(calculateAlleleDepth(a.getBases()[0]));
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
        attributes.put("AC", alleleCountModel.getMaximumLikelihoodIndex());
        attributes.put("AF", (double)alleleCountModel.getMaximumLikelihoodIndex()/alleleCountModel.getMaxAlleleCount());
        attributes.put("DP", pileup.getBases().length);
        attributes.put("AD", calculateAllelicDepths(getAlleles()));
        attributes.put("MQ", calculateMappingQualityRMS());
        attributes.put("MQ0", calculateMappingQualityZero());
        attributes.put("RD", errorModel.getReferenceDepth());
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
        alternateAlleles.add(Allele.create(refContext.getBase(), true));
        if (alleleCountModel.isVariant() ) {
            // add the most common alternate allele
            alternateAlleles.add(Allele.create(alleleCountModel.getAltBase(), false));
        }

        else {
            // todo -- not sure what to do with not confident calls
            //alternateAlleles.add(Allele.create(referenceSequenceBase, true)); // adding reference allele because I don't know how to handle this case.
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

    public boolean isVariant() {
        return (alleleCountModel.isVariant());
    }

    public boolean needToEmitCall() {
        return (isVariant() || (UAC.OutputMode == UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES));
    }

    public VariantContext getCallFromSite() {
        return siteVC;
    }

}
