package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 3/12/12
 * Time: 12:19 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class PoolGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {

    public static final String DUMMY_LANE = "Lane1";
    public static final String DUMMY_POOL = "Pool1";
    //protected Set<String> laneIDs;
    public enum Model {
        SNP,
        INDEL,
        POOLSNP,
        POOLINDEL,
        BOTH
    }

    final protected PoolCallerUnifiedArgumentCollection UAC;
    protected PoolGenotypePriors priors = null;

    protected PoolGenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        super(UAC,logger);
        if (UAC instanceof PoolCallerUnifiedArgumentCollection)
            this.UAC = (PoolCallerUnifiedArgumentCollection)UAC;
        else
            this.UAC = new PoolCallerUnifiedArgumentCollection(); // dummy copy

 /*       if (!(UAC instanceof PoolCallerUnifiedArgumentCollection)) 
            throw new ReviewedStingException("BUG: incorrect construction of PoolGenotypeLikelihoodsCalculationModel");

        PoolCallerUnifiedArgumentCollection pUAC = (PoolCallerUnifiedArgumentCollection)UAC;
        this.treatAllReadsAsSinglePool = pUAC.TREAT_ALL_READS_AS_SINGLE_POOL;
        this.referenceSampleName = pUAC.referenceSampleName;
        //this.laneIDs = laneIDs;
        this.minQualityScore = pUAC.minQualityScore;
        this.maxQualityScore = pUAC.maxQualityScore;
        this.sitePhredScalePrior = pUAC.phredScaledPrior;
        this.samplesPerPool = pUAC.nSamplesPerPool;
        this.minPower = pUAC.minPower;
        this.referenceSampleRod = pUAC.referenceSampleRod;
   */     
    }


    protected Collection<Allele> getTrueAlleles(final RefMetaDataTracker tracker,
                                              final ReferenceContext ref,
                                              Map<String,AlignmentContext> contexts) {
        // Get reference base from VCF or Reference
        AlignmentContext context = contexts.get(UAC.referenceSampleName);

        VariantContext referenceSampleVC = tracker.getFirstValue(UAC.referenceSampleRod, context.getLocation());

        ArrayList<Allele> trueReferenceAlleles = new ArrayList<Allele>();

        // Site is not a variant, take from the reference
        if (referenceSampleVC == null) {
            trueReferenceAlleles.add(Allele.create(ref.getBase(),true));
        }
        // Site has a VCF entry -- is variant
        else {
            Genotype referenceGenotype = referenceSampleVC.getGenotype(UAC.referenceSampleName);
            List<Allele> referenceAlleles = referenceGenotype.getAlleles();
            for (Allele allele : referenceAlleles) {
                if (!trueReferenceAlleles.contains(allele))
                    trueReferenceAlleles.add(allele);
            }
        }
        return trueReferenceAlleles;
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
    public static Set<String> parseLaneIDs(Collection<String> readGroups) {
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
    public static String getLaneIDFromReadGroupString(String readGroupID) {
        String [] parsedID = readGroupID.split("\\.");
        return parsedID[0] + "." + parsedID[1];
    }

}
