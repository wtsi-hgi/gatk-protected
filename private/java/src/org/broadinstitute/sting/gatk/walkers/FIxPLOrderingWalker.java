package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 7/21/11
 * Time: 10:31 AM
 * To change this template use File | Settings | File Templates.
 */
public class FIxPLOrderingWalker extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;
    private final String variantRodName = "variant";

    public void initialize() {
        // Initialize VCF header
        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(variantRodName);

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);

         vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));

    }
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
         if ( tracker == null )
             return 0;

         Collection<VariantContext> vcs = tracker.getValues(VariantContext.class, variantRodName, context.getLocation());

         if ( vcs == null || vcs.size() == 0) {
             return 0;
         }

         for (VariantContext vc : vcs) {
             if (vc.isIndel() && !vc.isBiallelic())
                 vc = modifyGLs(vc);
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

    private VariantContext modifyGLs(VariantContext vc) {
        int numAlleles = vc.getNAlleles();
        Map<String,Genotype> genotypes = new HashMap<String,Genotype> (vc.getGenotypes());
        for (String sample: genotypes.keySet()) {
            Genotype g = genotypes.get(sample);
            if (!g.hasLikelihoods())
                continue;
            GenotypeLikelihoods gls =  g.getLikelihoods();
            double[] glvec = gls.getAsVector();
            double[][] glmatrix = new double[numAlleles][numAlleles];

            // read first in column-wide ordering
            int k=0;
            for (int i=0; i < numAlleles; i++) {
                for (int j=i; j < numAlleles; j++) {
                    glmatrix[i][j] = glvec[k++];
                }
            }

            k=0;
            // now write in correct ordering
            for (int j=0; j < numAlleles; j++) {
                for (int i=0; i <= j; i++){
                    glvec[k++] = glmatrix[i][j];
                }
            }
            HashMap<String,Object> attributes = new HashMap<String,Object>(g.getAttributes());
            //GenotypeLikelihoods likelihoods = new GenotypeLikelihoods(GL.getLikelihoods());
            GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(glvec);
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, likelihoods);

            Genotype newg = new Genotype(sample, g.getAlleles(),g.getNegLog10PError(),g.getFilters(),attributes,g.isPhased());
            genotypes.put(sample,newg);

        }

        return new VariantContext( vc.getSource(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), genotypes, vc.getNegLog10PError(), vc.getFilters(), vc.getAttributes());
    }
}