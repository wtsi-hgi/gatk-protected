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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.indels.PairHMMIndelErrorModel;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

public class PoolIndelGenotypeLikelihoodsCalculationModel extends PoolGenotypeLikelihoodsCalculationModel {
    private final int HAPLOTYPE_SIZE;

    private final int minIndelCountForGenotyping;
    private final boolean getAlleleListFromVCF;

    private boolean DEBUG = false;
    private final boolean doMultiAllelicCalls = true;
    private boolean ignoreSNPAllelesWhenGenotypingIndels = false;
    private PairHMMIndelErrorModel pairModel;

    private static ThreadLocal<HashMap<PileupElement, LinkedHashMap<Allele, Double>>> indelLikelihoodMap =
            new ThreadLocal<HashMap<PileupElement, LinkedHashMap<Allele, Double>>>() {
                protected synchronized HashMap<PileupElement, LinkedHashMap<Allele, Double>> initialValue() {
                    return new HashMap<PileupElement, LinkedHashMap<Allele, Double>>();
                }
            };

    private LinkedHashMap<Allele, Haplotype> haplotypeMap;

    // gdebug removeme
    // todo -cleanup
    private GenomeLoc lastSiteVisited;
    private ArrayList<Allele> alleleList;

    static {
        indelLikelihoodMap.set(new HashMap<PileupElement, LinkedHashMap<Allele, Double>>());
    }


    protected PoolIndelGenotypeLikelihoodsCalculationModel(final UnifiedArgumentCollection UAC, final Logger logger) {
        super(UAC, logger);


        pairModel = new PairHMMIndelErrorModel(UAC.INDEL_GAP_OPEN_PENALTY, UAC.INDEL_GAP_CONTINUATION_PENALTY,
                UAC.OUTPUT_DEBUG_INDEL_INFO, !UAC.DONT_DO_BANDED_INDEL_COMPUTATION);
        alleleList = new ArrayList<Allele>();
        getAlleleListFromVCF = UAC.GenotypingMode == GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
        minIndelCountForGenotyping = UAC.MIN_INDEL_COUNT_FOR_GENOTYPING;
        HAPLOTYPE_SIZE = UAC.INDEL_HAPLOTYPE_SIZE;
        DEBUG = UAC.OUTPUT_DEBUG_INDEL_INFO;

        haplotypeMap = new LinkedHashMap<Allele, Haplotype>();
        ignoreSNPAllelesWhenGenotypingIndels = UAC.IGNORE_SNP_ALLELES;
    }



    private final static EnumSet<VariantContext.Type> allowableTypes = EnumSet.of(VariantContext.Type.INDEL, VariantContext.Type.MIXED);

    public VariantContext getLikelihoods(final RefMetaDataTracker tracker,
                                         final ReferenceContext ref,
                                         final Map<String, AlignmentContext> contexts,
                                         final AlignmentContextUtils.ReadOrientation contextType,
                                         final List<Allele> alternateAllelesToUse,
                                         final boolean useBAQedPileup,
                                         final GenomeLocParser locParser) {

        if (tracker == null)
            return null;

        //dummy, no functionality yet
        return null;
    }


}