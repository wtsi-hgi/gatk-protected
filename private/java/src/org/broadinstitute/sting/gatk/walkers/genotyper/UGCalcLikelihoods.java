/*
 * Copyright (c) 2010 The Broad Institute
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

import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.downsampling.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;

import java.util.*;


/**
 * Uses the UG engine to determine per-sample genotype likelihoods and emits them as a VCF (using PLs).
 * Absolutely not supported or recommended for public use.
 * Run this as you would the UnifiedGenotyper, except that you must additionally pass in a VCF bound to
 * the name 'allele' so we know which alternate allele to use at each site.
 */
@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = ReadTransformer.ApplicationTime.ON_INPUT)
@Reference(window=@Window(start=-200,stop=200))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.BY_SAMPLE, toCoverage=250)
public class UGCalcLikelihoods extends LocusWalker<List<VariantCallContext>, Integer> implements TreeReducible<Integer> {

    @ArgumentCollection private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // control the output
    @Output(doc="File to which variants should be written",required=true)
    protected VariantContextWriter writer = null;

    // the calculation arguments
    private UnifiedGenotyperEngine UG_engine = null;

    // enable deletions in the pileup
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    public void initialize() {
        // get all of the unique sample names
        Set<String> samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());

        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, samples, GATKVariantContextUtils.DEFAULT_PLOIDY);

        // initialize the header
        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        writer.writeHeader(new VCFHeader(headerInfo, samples)) ;
    }

    public List<VariantCallContext> map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        List<VariantCallContext> vccList = new LinkedList<VariantCallContext>();

        List<RefMetaDataTracker> useTrackers = new LinkedList<RefMetaDataTracker>();
        // Allow for multiple records in UAC.alleles, even at same locus:
        if ( UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            for (VariantContext vc : tracker.getValues(UAC.alleles, rawContext.getLocation()))
                useTrackers.add(new MatchFirstLocRefAltRefMetaDataTracker(tracker, vc));
        }
        else
            useTrackers.add(tracker);

        for (RefMetaDataTracker t : useTrackers) {
            final Map<String, org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = new HashMap<String, org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap>();
            VariantContext call = UG_engine.calculateLikelihoods(t, refContext, rawContext, perReadAlleleLikelihoodMap);
            if (call != null)
                vccList.add(new VariantCallContext(call, true));
        }

        return vccList;
    }

    public Integer reduceInit() { return 0; }

    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    public Integer reduce(List<VariantCallContext> value, Integer sum) {
        if ( value == null )
            return sum;

        try {
            for (VariantCallContext vcc : value)
                writer.add(vcc);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(e.getMessage());
        }

        return sum + value.size();
    }

    public void onTraversalDone(Integer sum) {
        logger.info(String.format("Visited variants: %d", sum));
    }
}