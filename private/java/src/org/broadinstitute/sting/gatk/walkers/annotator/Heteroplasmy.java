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
package org.broadinstitute.sting.gatk.walkers.annotator;

import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.PoolAFCalculationModel;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineCount;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

public class Heteroplasmy extends InfoFieldAnnotation {
    public static final String NUM_VARIANT_SAMPLES_KEY = "NV";
    public static final String HETEROPLASMY_KEY = "HP";

    private String[] keyNames = { NUM_VARIANT_SAMPLES_KEY, HETEROPLASMY_KEY };
    private VCFInfoHeaderLine[] descriptions = { new VCFInfoHeaderLine(NUM_VARIANT_SAMPLES_KEY, 1, VCFHeaderLineType.Integer, "Humber of samples with a non-reference genotype"),
            new VCFInfoHeaderLine(HETEROPLASMY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Mean heteroplasmy (dosage fraction of variant allele) across non-variant samples")};

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        final GenotypesContext genotypes = vc.getGenotypes();

        // todo- avoid complications from
        if (!vc.isBiallelic())
            return null;

        int numVariantSamples = 0;
        double heteroplasmySum = 0.0;
        
        for (Genotype g: genotypes) {
            if (g.hasAttribute(PoolAFCalculationModel.MAXIMUM_LIKELIHOOD_AC_KEY)) {
                String s = g.getAttributeAsString(PoolAFCalculationModel.MAXIMUM_LIKELIHOOD_AC_KEY,"");
                int numAlts = Integer.valueOf(s);
                if (numAlts>0) numVariantSamples++;

                // AF will be per-pool heteroplasmy
                s = g.getAttributeAsString(PoolAFCalculationModel.MAXIMUM_LIKELIHOOD_AF_KEY,"");
                double af = Double.valueOf(s);
                if (numAlts>0) heteroplasmySum += af;
         }
       }

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(NUM_VARIANT_SAMPLES_KEY, String.format("%d", numVariantSamples));
        if (numVariantSamples > 0)
            map.put(HETEROPLASMY_KEY, String.format("%4.3f",heteroplasmySum/(double)numVariantSamples));
        return map;
    }

    public List<String> getKeyNames() {
        return Arrays.asList(keyNames);
    }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(descriptions); }
     /*
    private static ArrayList<Integer> parseString(String s) {
        String[] pieces = s.split(",");
        ArrayList<Integer> vals = new ArrayList<Integer>();
        for (String sp : pieces) {
            vals.add(Integer.valueOf(sp));

        }
        return vals;
    }    */
}
