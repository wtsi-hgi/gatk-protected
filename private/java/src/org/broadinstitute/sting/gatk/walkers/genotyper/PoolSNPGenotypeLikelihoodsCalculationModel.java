package org.broadinstitute.sting.gatk.walkers.genotyper;
/*
 * Copyright (c) 2010.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

public class PoolSNPGenotypeLikelihoodsCalculationModel extends PoolGenotypeLikelihoodsCalculationModel {

    //private final boolean useAlleleFromVCF;

//    private final double[] likelihoodSums = new double[4];

    protected PoolSNPGenotypeLikelihoodsCalculationModel( UnifiedArgumentCollection UAC,  Logger logger) {
        super(UAC, logger);

    }

    public VariantContext getLikelihoods(final RefMetaDataTracker tracker,
                                         final ReferenceContext ref,
                                         Map<String, AlignmentContext> contexts,
                                         final AlignmentContextUtils.ReadOrientation contextType,
                                         final List<Allele> alternateAllelesToUse,
                                         final boolean useBAQedPileup,
                                         final GenomeLocParser locParser) {


        final byte refBase = ref.getBase();
        final int indexOfRefBase = BaseUtils.simpleBaseToBaseIndex(refBase);
        
        Collection<Allele> trueReferenceAlleles = getTrueAlleles(tracker, ref, contexts);

        // start making the VariantContext
        final GenomeLoc loc = ref.getLocus();
        final List<Allele> alleles = new ArrayList<Allele>();
        alleles.add(Allele.create(refBase, true));
        final VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), loc.getStop(), alleles);

        // Build error model for site based on reference sample, and keep stratified for each lane.
        HashMap<String, ErrorModel> perLaneErrorModels = new HashMap<String, ErrorModel>();
        AlignmentContext refContext = contexts.get(UAC.referenceSampleName);

        ReadBackedPileup refPileup = null;
        if (refContext != null && refContext.hasBasePileup()) {
            refPileup = refContext.getBasePileup();

            Set<String> laneIDs = new TreeSet<String>();
            if (UAC.TREAT_ALL_READS_AS_SINGLE_POOL || UAC.IGNORE_LANE_INFO)
                laneIDs.add(PoolGenotypeLikelihoodsCalculationModel.DUMMY_LANE);
            else
                laneIDs = parseLaneIDs(refPileup.getReadGroups());
            // build per-lane error model for all lanes present in ref sample
            for (String laneID : laneIDs) {
                // get reference pileup for this lane
                ReadBackedPileup refLanePileup = refPileup;
                // subset for this lane
                if (refPileup != null && !(UAC.TREAT_ALL_READS_AS_SINGLE_POOL || UAC.IGNORE_LANE_INFO))
                    refLanePileup = refPileup.getPileupForLane(laneID);

                //ReferenceSample referenceSample = new ReferenceSample(UAC.referenceSampleName, refLanePileup, trueReferenceAlleles);
                perLaneErrorModels.put(laneID, new ErrorModel(UAC.minQualityScore, UAC.maxQualityScore, UAC.phredScaledPrior,  refLanePileup, trueReferenceAlleles, UAC.minPower));
            }


        }
        else
            return null;

        // calculate the pool likelihoods
        if (UAC.TREAT_ALL_READS_AS_SINGLE_POOL) {
            AlignmentContext mergedContext = AlignmentContextUtils.joinContexts(contexts.values());
            Map<String,AlignmentContext> newContext = new HashMap<String,AlignmentContext>();
            newContext.put(DUMMY_POOL,mergedContext);
            contexts = newContext;
        }
        ArrayList<PoolGenotypeData> GLs = new ArrayList<PoolGenotypeData>(contexts.size());
        final List<Allele> allAlleles = new ArrayList<Allele>(alleles);  // this contains only ref Allele up to now

        if (alternateAllelesToUse == null) {
            // add all possible alleles 
            for (byte b: BaseUtils.BASES) {
                if (refBase != b)
                    allAlleles.add(Allele.create(b));
            }
        }  else {
            allAlleles.addAll(alternateAllelesToUse);
        }
        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            // skip reference sample
            if (sample.getKey().equals(UAC.referenceSampleName))
                continue;

            ReadBackedPileup pileup = AlignmentContextUtils.stratify(sample.getValue(), contextType).getBasePileup();
            if ( useBAQedPileup )
                pileup = createBAQedPileup( pileup );

            // create the GenotypeLikelihoods object
            final PoolSNPGenotypeLikelihoods GL = new PoolSNPGenotypeLikelihoods(allAlleles, null, UAC.nSamplesPerPool*2, perLaneErrorModels, UAC.IGNORE_LANE_INFO);
            final int nGoodBases = GL.add(pileup, true, true, UAC.MIN_BASE_QUALTY_SCORE);
            if ( nGoodBases > 0 )
                GLs.add(new PoolGenotypeData(sample.getKey(), GL, getFilteredDepth(pileup)));
        }

        // find the alternate allele(s) that we should be using
        if ( alternateAllelesToUse != null ) {
            alleles.addAll(alternateAllelesToUse);
        } else if ( UAC.GenotypingMode == GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            final VariantContext vc = UnifiedGenotyperEngine.getVCFromAllelesRod(tracker, ref, ref.getLocus(), true, logger, UAC.alleles);

            // ignore places where we don't have a SNP
            if ( vc == null || !vc.isSNP() )
                return null;

            alleles.addAll(vc.getAlternateAlleles());
        } else {

            alleles.addAll(determineAlternateAlleles( GLs));

            // if there are no non-ref alleles...
            if ( alleles.size() == 1 ) {
                // if we only want variants, then we don't need to calculate genotype likelihoods
                if ( UAC.OutputMode == UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                    return builder.make();

                // otherwise, choose any alternate allele (it doesn't really matter)
                alleles.add(Allele.create(BaseUtils.baseIndexToSimpleBase(indexOfRefBase == 0 ? 1 : 0)));
            }
        }

        // create the alternate alleles and the allele ordering (the ordering is crucial for the GLs)
        final int numAlleles = alleles.size();
        final int numAltAlleles = numAlleles - 1;



        builder.alleles(alleles);

        // create the genotypes; no-call everyone for now
        final GenotypesContext genotypes = GenotypesContext.create();
        final List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        for ( PoolGenotypeData sampleData : GLs ) {
            // extract from multidimensional array
            final double[] myLikelihoods = PoolGenotypeLikelihoods.subsetToAlleles(sampleData.GL.getLikelihoods(),sampleData.GL.numChromosomes,
            allAlleles, alleles);

            // normalize in log space so that max element is zero.
            final GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(MathUtils.normalizeFromLog10(myLikelihoods, false, true));

            final HashMap<String, Object> attributes = new HashMap<String, Object>();
            attributes.put(VCFConstants.DEPTH_KEY, sampleData.depth);
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, likelihoods);
            genotypes.add(new Genotype(sampleData.name, noCall, Genotype.NO_LOG10_PERROR, null, attributes, false));
        }

        return builder.genotypes(genotypes).make();
    }


    public ReadBackedPileup createBAQedPileup( final ReadBackedPileup pileup ) {
        final List<PileupElement> BAQedElements = new ArrayList<PileupElement>();
        for( final PileupElement PE : pileup ) {
            final PileupElement newPE = new BAQedPileupElement( PE );
            BAQedElements.add( newPE );
        }
        return new ReadBackedPileupImpl( pileup.getLocation(), BAQedElements );
    }

    public class BAQedPileupElement extends PileupElement {
        public BAQedPileupElement( final PileupElement PE ) {
            super(PE.getRead(), PE.getOffset(), PE.isDeletion(), PE.isBeforeDeletedBase(), PE.isAfterDeletedBase(), PE.isBeforeInsertion(), PE.isAfterInsertion(), PE.isNextToSoftClip());
        }

        @Override
        public byte getQual( final int offset ) { return BAQ.calcBAQFromTag(getRead(), offset, true); }
    }

}
