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

package org.broadinstitute.sting.gatk.walkers.callset_assessment;

import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Map;

/**
 * Assesses GLs at truth sites.
 * Use --variant and --truth
 */
public class AssessGLConcordance extends RodWalker<Integer, Integer>   {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Input(fullName="truth", shortName = "truth", doc="Input VCF truth file", required=true)
    public RodBinding<VariantContext> truthTrack;


    private double totalError = 0.0;
    private int totalNumSites = 0;
    private int nSamples = 0;
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        VariantContext variant = tracker.getFirstValue(variantCollection.variants, context.getLocation());
        if ( variant == null )
            return 0;

        VariantContext truth = tracker.getFirstValue(truthTrack, context.getLocation());
        if ( truth == null )
            return 0;

        int nAlleles = variant.getNAlleles();

        if (truth.getNAlleles() != nAlleles)
            return 0;

        double[] emptyGL = new double[nAlleles*(nAlleles+1)/2];

        double siteError = 0.0;
        nSamples = variant.getNSamples();

        for ( final String sample: variant.getSampleNames() ) {
            Genotype evalGenotype = variant.getGenotype(sample);

            if ( !truth.hasGenotype(sample) )
                 continue;

            double[] evalGL, truthGL;

            if ( evalGenotype.isNoCall() )  {
                evalGL = emptyGL;
            }
            else
                evalGL = evalGenotype.getLikelihoods().getAsVector();


            Genotype truthGenotype = truth.getGenotype(sample);
            if ( truthGenotype.isNoCall() )
                truthGL = emptyGL;

            else
                truthGL = truthGenotype.getLikelihoods().getAsVector();


            double[] normalizedEvalGLs = MathUtils.normalizeFromLog10(evalGL);
            double[] normalizedTruthGLs = MathUtils.normalizeFromLog10(truthGL);
   /*
            System.out.print("\n"+sample+" Eval:");
            for (int k=0; k < normalizedEvalGLs.length; k++)
                System.out.format("%4.3f ",normalizedEvalGLs[k]);
            System.out.print("\nTruth:");
            for (int k=0; k < normalizedTruthGLs.length; k++)
                System.out.format("%4.3f ",normalizedTruthGLs[k]);
      */
            // compute MSE between eval and truth GLs
            double err = 0.0, norm = 0.0;

            for (int k=0; k < normalizedEvalGLs.length; k++) {
                err += (normalizedEvalGLs[k]-normalizedTruthGLs[k])*(normalizedEvalGLs[k]-normalizedTruthGLs[k]);
              //  norm +=  normalizedTruthGLs[k] * normalizedTruthGLs[k];
            }

            siteError += err; ///norm;
        }

        totalNumSites++;
        totalError += siteError; // /nSamples;
        return 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {
        System.out.format("Number of Sites:%d Number of Samples:%d Total Error:%4.1f\n", totalNumSites, nSamples, totalError);
    }
}
