/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Emits a PED and INFO field file from a VCF file
 *
 * <p>
 *     Takes a VCF file and writes out a PED/INFO file suitable for loading into Haploview.
 * </p>
 *
 * <h2>Input</h2>
 * <p>A VCF file with genotypes</p>
 *
 * <h2>Output</h2>
 * <p>A PED and INFO files.
 * See http://www.broadinstitute.org/science/programs/medical-and-population-genetics/haploview/input-file-formats-0 for
 * information on these files.  This info file is a 'Marker Information File'
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *     -T $WalkerName -V my.vcf -pedOut my.ped -infoOut my.info
 * </pre>
 *
 * @author Mark DePristo
 * @since 2011
 */
public class VariantsToPED extends RodWalker<VariantContext, Collection<VariantContext>> {
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(shortName = "pedOut", fullName = "pedOut", doc="File to which results should be written",required=true)
    protected PrintStream pedOutput;

    @Output(shortName = "infoOut", fullName = "infoOut", doc="File to which results should be written",required=true)
    protected PrintStream infoOutput;

    @Override
    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return null;

        final VariantContext vc = tracker.getFirstValue(variantCollection.variants, ref.getLocus());
        if ( vc != null && vc.isSNP() && vc.isBiallelic() )
            return vc;
        else
            return null;
    }


    @Override
    public Collection<VariantContext> reduceInit() {
        return new ArrayList<VariantContext>();
    }

    @Override
    public Collection<VariantContext> reduce(VariantContext singleton, Collection<VariantContext> sum) {
        if ( singleton != null )
            sum.add(singleton);
        return sum;
    }

    @Override
    public void onTraversalDone(Collection<VariantContext> sum) {
        // write the per sample genotypes
        final List<String> inputNames = Arrays.asList(variantCollection.variants.getName());
        final Collection<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), inputNames);

        for ( final String sample : samples ) {
            List<String> genotypes = new ArrayList<String>(sum.size());
            for ( final VariantContext vc : sum ) {
                final Genotype g = vc.getGenotype(sample);
                if ( g.isNoCall() )
                    genotypes.add("0 0");
                else
                    genotypes.add(String.format("%s %s", g.getAllele(0).getBaseString(), g.getAllele(1).getBaseString()));
            }

            pedOutput.printf("%s\t%s\t0\t0\t1\t0\t%s%n", sample, sample, Utils.join("\t", genotypes));
        }

        // write the info field
        for ( final VariantContext vc : sum ) {
            infoOutput.printf("%s\t%d%n", vc.hasID() ? vc.getID() : "SNP-" + vc.getStart(), vc.getStart());
        }
    }
}
