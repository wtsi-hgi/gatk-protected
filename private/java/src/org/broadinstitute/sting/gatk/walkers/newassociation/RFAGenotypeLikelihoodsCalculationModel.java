package org.broadinstitute.sting.gatk.walkers.newassociation;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.MultiallelicGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.newassociation.features.old.BinaryFeatureAggregator;
import org.broadinstitute.sting.gatk.walkers.newassociation.features.old.ReadFeatureAggregator;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.Sample;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A genotype likelihood calculation model for use with read feature association. In order to maintain
 * 'genotype-free'-ness, ploidy is empirically estimated via a goodness-of-fit test [though perhaps
 * it is something that may be marginalized over], and genotype likelihoods are likelihoods of allele counts
 * under the calculated empirical ploidy. These likelihoods are marginalized over by the ExactAFCalculation to generate
 * P[AC=Y|D], which can then be marginalized over to give a skew. Thus one could forsee marginalizing additionally over
 * the ploidy: Pr[Skew] = Sum_P Pr[Skew|P]Pr[P|D] (and Pr[P|D] = Pr[D|P]/Sum_P{Pr[P|D]} -- etc.
 */
public class RFAGenotypeLikelihoodsCalculationModel  {
    private static double FLAT_ERROR = 1e-4;

    public RFAGenotypeLikelihoodsCalculationModel() {

    }

    public Map<String,GenotypeLikelihoods> getLikelihoods(RefMetaDataTracker tracker, ReferenceContext ref, Map<String,BinaryFeatureAggregator> aggregators, int ploidy) {
        Map<String,GenotypeLikelihoods> likelihoods = new HashMap<String,GenotypeLikelihoods>(aggregators.size());
        for ( Map.Entry<String,BinaryFeatureAggregator> aggregatorEntry : aggregators.entrySet() ) {
            double[] lHood = new double[1+ploidy];
            for ( int ac = 0; ac < lHood.length; ac ++ ) {
                lHood[ac] = calculateLikelihoodFromAggregator(ac,ploidy,aggregatorEntry.getValue());
            }
            likelihoods.put(aggregatorEntry.getKey(),GenotypeLikelihoods.fromLog10Likelihoods(MathUtils.normalizeFromLog10(lHood)));
        }

        return likelihoods;
    }

    public double calculateLikelihoodFromAggregator(int alleleCount, int ploidy, BinaryFeatureAggregator aggregator) {
        // todo -- aggregator should have some granularized metric for quality, and there should be an
        //      -- inner loop a la DiploidSNPGenotypeLikelihoods.java to map this to a non-flat error
        // genotypes are assumed to be multiploid, but bi-allelic, and
        // Pr[Read not showing variant|chrom is variant] = FLAT_ERROR, and the reverse
        // math: because of normalization at the very end, we can ignore typical normalizing coefficients - binomials, etc

        double p = aggregator.getnAberrant()*( (1.0-FLAT_ERROR)*(0.0+alleleCount)/ploidy + FLAT_ERROR * (ploidy-alleleCount-0.0)/ploidy);
        p += aggregator.getnNotAberrant()*( (FLAT_ERROR)*(0.0+alleleCount)/ploidy + (1.0-FLAT_ERROR) * (ploidy-alleleCount-0.0)/ploidy);

        return Math.log10(p);
    }

    public Map<String,Genotype> getLikelihoods(RefMetaDataTracker tracker, ReferenceContext ref, Map<String,BinaryFeatureAggregator> aggregators) {
        // use a likelihood ratio test to determine the ploidy ( <= 10) that "best" represents the data
        int bestPloidy = 2;
        int maxPloidy = 2; // todo -- raise me after generalizing exact AF calculation
        double bestLikelihood = Double.NEGATIVE_INFINITY;
        Map<String,GenotypeLikelihoods> bestLikelihoods = null;
        for ( int ploidy = 2; ploidy <= maxPloidy; ploidy++ ) {
            Map<String,GenotypeLikelihoods> likelihoods = getLikelihoods(tracker,ref,aggregators,ploidy);
            // now just sum
            double lhood = 0.0;
            for ( Map.Entry<String,GenotypeLikelihoods> entry : likelihoods.entrySet() ) {
                for ( double d : entry.getValue().getAsVector() ) {
                    lhood += d;
                }
            }

            if ( 2*(lhood - bestLikelihood) > 3.8414 ) { // significant
                bestPloidy = ploidy;
                bestLikelihood = lhood;
                bestLikelihoods = likelihoods;
            }
        }

        // normalize likelihoods and phred-scale
        for ( Map.Entry<String,GenotypeLikelihoods> entry : bestLikelihoods.entrySet() ) {
            double[] lHoods = entry.getValue().getAsVector();
        }

        // cast the likelihoods into a genotype object
        Map<String,Genotype> bestGenos = new HashMap<String,Genotype>(bestLikelihoods.size());
        ArrayList<Allele> alleles = new ArrayList<Allele>(2);
        alleles.add(Allele.create(ref.getBase()));
        alleles.add(Allele.NO_CALL);
        for ( Map.Entry<String,GenotypeLikelihoods> entry : bestLikelihoods.entrySet() ) {
            HashMap<String,Object> attributes = new HashMap<String,Object>(1);
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY,entry.getValue());
            Genotype g = new Genotype(entry.getKey(),alleles,Genotype.NO_NEG_LOG_10PERROR,null,attributes,false);
            bestGenos.put(entry.getKey(),g);
        }

        return bestGenos;
    }
}
