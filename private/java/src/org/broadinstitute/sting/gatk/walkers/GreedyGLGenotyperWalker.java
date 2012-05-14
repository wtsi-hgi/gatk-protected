package org.broadinstitute.sting.gatk.walkers;


import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;


public class GreedyGLGenotyperWalker extends RodWalker<Integer, Integer>  implements TreeReducible<Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VariantContextWriter vcfWriter = null;

     @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public List<RodBinding<VariantContext>> variants;

    @Argument(fullName="excludeMultiallelics", shortName="excludeMultiallelics", doc="exclude multi allelic records", required=false)
    private Boolean excludeMultiallelics = false;


    private UnifiedGenotyperEngine engine;
    public void initialize() {
        // Initialize VCF header
        engine = new UnifiedGenotyperEngine(this.getToolkit(), new UnifiedArgumentCollection()); // dummy
        Set<String> rodNames = SampleUtils.getRodNamesWithVCFHeader(getToolkit(), null);

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);

        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));

    }
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> vcs = tracker.getValues(variants, context.getLocation());

        if ( vcs == null || vcs.size() == 0) {
            return 0;
        }

        int written = 0;
        for (VariantContext vc : vcs) {
            if (excludeMultiallelics && !vc.isBiallelic())
                continue;

            GenotypesContext genotypes = VariantContextUtils.assignDiploidGenotypes(vc);

            // finish constructing the resulting VC
            GenomeLoc loc = getToolkit().getGenomeLocParser().createGenomeLoc(vc);
            VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), loc.getStop(), vc.getAlleles());
            builder.log10PError(vc.getPhredScaledQual()/-10.0);
            builder.filters(vc.getFilters());
            if ( vc.hasReferenceBaseForIndel() ) {
                builder.referenceBaseForIndel(vc.getReferenceBaseForIndel());
            } else {
                builder.referenceBaseForIndel(ref.getBase());
            }
            builder.genotypes(genotypes);
            builder.attributes(vc.getAttributes());
            VariantContext vcCall = builder.make();
            builder = new VariantContextBuilder(vcCall);
            // re-compute chromosome counts
            VariantContextUtils.calculateChromosomeCounts(builder, false);

            vcfWriter.add(builder.make());
            written++;
        }
        return written;
    }
    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(Integer value, Integer sum) { return value + sum; }

    //@Override
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }


    public void onTraversalDone(Integer result) {
        logger.info(result + " records processed.");
    }

}