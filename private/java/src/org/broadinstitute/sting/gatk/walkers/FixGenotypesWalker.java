package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 7/21/11
 * Time: 10:31 AM
 * To change this template use File | Settings | File Templates.
 */
public class FixGenotypesWalker extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    public void initialize() {
        // Initialize VCF header
        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(variantCollection.variants.getName());

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);

         vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));

    }
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
         if ( tracker == null )
             return 0;

         Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());

         if ( vcs == null || vcs.size() == 0) {
             return 0;
         }

         for (VariantContext vc : vcs) {
             if (vc.isIndel())
                 vc = modifyGLs(vc, ref);

            if (vc.isPolymorphic())
                vcfWriter.add(vc);
         }
        return 1;
    }
    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(Integer value, Integer sum) { return value + sum; }

    public void onTraversalDone(Integer result) {
        logger.info(result + " records processed.");
    }

    private VariantContext modifyGLs(VariantContext vc, ReferenceContext ref) {
        Map<String,Genotype> genotypes = new HashMap<String,Genotype> (vc.getGenotypes());
        for (String sample: genotypes.keySet()) {
            Genotype g = genotypes.get(sample), newg = null;
            if (g.isCalled()) {


                GenotypeLikelihoods gls =  g.getLikelihoods();
                double[] glvec = gls.getAsVector();

                boolean isAllZeros = true;
                for (int k=0; k < glvec.length; k++) {
                    if (glvec[k] != 0) {
                        isAllZeros = false;
                        break;
                    }
                }


                if (isAllZeros) {
                    ArrayList<Allele> a = new ArrayList<Allele>();
                    a.add(Allele.NO_CALL);
                    a.add(Allele.NO_CALL);
                    newg = new Genotype(sample,a);
                }
                else
                    newg = g;
            }
            else
                newg = g;
            genotypes.put(sample,newg);

        }

        VariantContext sub =  new VariantContext( vc.getSource(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), genotypes, vc.getNegLog10PError(), vc.getFilters(), vc.getAttributes(), ref.getBase());
        HashMap<String, Object> attributes = new HashMap<String, Object>(sub.getAttributes());

        VariantContextUtils.calculateChromosomeCounts(sub, attributes, false);

        return VariantContext.modifyAttributes(sub, attributes);

    }
}