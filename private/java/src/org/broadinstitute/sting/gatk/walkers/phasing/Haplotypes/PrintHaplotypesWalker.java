/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.phasing.Haplotypes;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Allows;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.phasing.AnnotateTrioPhasingInheritanceNoRecombinationWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static org.broadinstitute.sting.utils.codecs.vcf.VCFUtils.getVCFHeadersFromRods;

/**
 * Walks along all variant ROD loci and uses the phase information to divide up the genome into phased segments.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})

public class PrintHaplotypesWalker extends RodWalker<Integer, Integer> {
    @Output(doc = "File to which results should be written", required = true)
    protected PrintStream out;

    @Argument(doc = "sample to emit", required = false)
    protected String sample = null;

    private boolean prevHadInheritanceInfo;
    private GenomeLoc prevLoc;

    public void initialize() {
        this.prevHadInheritanceInfo = true;
        this.prevLoc = null;
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public Integer reduceInit() {
        return 0;
    }

    /**
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        GenomeLoc curLocus = ref.getLocus();
        if (prevLoc == null || prevLoc.compareContigs(curLocus) != 0)
            prevHadInheritanceInfo = true;

        Collection<VariantContext> vcs = tracker.getAllVariantContexts(ref, curLocus);
        for (VariantContext vc : vcs) {
            if (vc.isFiltered())
                continue;

            Genotype gt = vc.getGenotype(sample);
            if (gt == null || gt.isFiltered() || gt.isNoCall() || !gt.isHet())
                continue;

            String inheritance = gt.getAttributeAsStringNoException(AnnotateTrioPhasingInheritanceNoRecombinationWalker.INHERITANCE_KEY);
            if (!gt.isPhased() && (inheritance == null || !prevHadInheritanceInfo))
                out.println("PHASE_BREAK");

            String[] sources;
            if (inheritance != null) {
                sources = inheritance.split("\\|");
            }
            else {
                sources = new String[gt.getPloidy()];
                for (int i = 0; i < sources.length; i++) {
                    sources[i] = new Integer(i+1).toString();
                }
            }

            out.print(curLocus);
            int ind = 0;
            for (Allele all : gt.getAlleles()) {
                out.print("\t" + sources[ind] + ": " + (all.isReference() ? 0 : 1));
                ind++;
            }
            out.println();

            prevHadInheritanceInfo = (inheritance != null);
            prevLoc = curLocus;
        }

        return 1;
    }

    public Integer reduce(Integer addIn, Integer runningCount) {
        if (addIn == null)
            addIn = 0;

        return runningCount + addIn;
    }

    /**
     * @param result the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(Integer result) {
        System.out.println("map was called " + result + " times.");
    }
}