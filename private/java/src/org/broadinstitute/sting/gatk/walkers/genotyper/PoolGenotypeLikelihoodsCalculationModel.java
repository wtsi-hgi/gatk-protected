package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.collections.Pair;
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

    protected PoolGenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        super(UAC,logger);
        if (UAC instanceof PoolCallerUnifiedArgumentCollection)
            this.UAC = (PoolCallerUnifiedArgumentCollection)UAC;
        else
            this.UAC = new PoolCallerUnifiedArgumentCollection(); // dummy copy

    }


    protected Collection<Allele> getTrueAlleles(final RefMetaDataTracker tracker,
                                              final ReferenceContext ref,
                                              Map<String,AlignmentContext> contexts) {
        // Get reference base from VCF or Reference
        AlignmentContext context = contexts.get(UAC.referenceSampleName);
        ArrayList<Allele> trueReferenceAlleles = new ArrayList<Allele>();

        VariantContext referenceSampleVC = null;

        if (tracker != null)
            referenceSampleVC = tracker.getFirstValue(UAC.referenceSampleRod, context.getLocation());

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
       // System.out.println(readGroupID);
        String [] parsedID = readGroupID.split("\\.");
        if (parsedID.length > 1)
            return parsedID[0] + "." + parsedID[1];
        else
            return parsedID[0] + ".0";
    }


    /** Wrapper class that encapsulates likelihood object and sample name
     *
     */
    protected static class PoolGenotypeData {

        public final String name;
        public final PoolSNPGenotypeLikelihoods GL;
        public final int depth;

        public PoolGenotypeData(final String name, final PoolSNPGenotypeLikelihoods GL, final int depth) {
            this.name = name;
            this.GL = GL;
            this.depth = depth;
        }
    }

    // determines the alleles to use
    protected List<Allele> determineAlternateAlleles(final List<PoolGenotypeData> sampleDataList) {


        final int REFERENCE_IDX = 0;
        final List<Allele> allAlleles = sampleDataList.get(0).GL.getAlleles();
        double[] likelihoodSums = new double[allAlleles.size()];

        // based on the GLs, find the alternate alleles with enough probability
        for ( PoolGenotypeData sampleData : sampleDataList ) {
            final Pair<int[],Double> mlACPair = sampleData.GL.getMostLikelyACCount();
            final double topLogGL = mlACPair.second;

            if (sampleData.GL.getAlleles().size() != allAlleles.size())
                throw new ReviewedStingException("BUG: inconsistent size of alleles!");

            // ref allele is always first in array list
            if (sampleData.GL.alleles.get(0).isNonReference())
                throw new ReviewedStingException("BUG: first allele in list is not reference!");

            double refGL = sampleData.GL.getLikelihoods()[REFERENCE_IDX];

            // check if maximum likelihood AC is all-ref for current pool. If so, skip
            if (mlACPair.first[REFERENCE_IDX] == sampleData.GL.numChromosomes)
                continue;

            // most likely AC is not all-ref: for all non-ref alleles, add difference of max likelihood and all-ref likelihood
            for (int i=0; i < mlACPair.first.length; i++) {
                if (i==REFERENCE_IDX) continue;

                if (mlACPair.first[i] > 0)
                    likelihoodSums[i] += topLogGL - refGL;

            }
        }

        final List<Allele> allelesToUse = new ArrayList<Allele>();
        for ( int i = 0; i < likelihoodSums.length; i++ ) {
            if ( likelihoodSums[i] > 0.0 )
                allelesToUse.add(allAlleles.get(i));
        }

        return allelesToUse;
    }

}
