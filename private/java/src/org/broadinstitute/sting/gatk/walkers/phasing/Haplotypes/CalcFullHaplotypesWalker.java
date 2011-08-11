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
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Allows;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static org.broadinstitute.sting.utils.codecs.vcf.VCFUtils.getVCFHeadersFromRods;

/**
 * Walks along all variant ROD loci and uses the phase information to divide up the genome into phased segments.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})

public class CalcFullHaplotypesWalker extends RodWalker<Integer, Integer> {
    @Output(doc = "File to which results should be written", required = true)
    protected PrintStream out;

    @Argument(doc = "sample to emit", required = false)
    protected String sample = null;

    @Argument(doc = "only include physically-phased results", required = false)
    protected boolean requirePQ = false;

    private HaplotypeTracker hapTracker = null;

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    public void initialize() {
        Set<String> samples;
        if (sample == null) {
            samples = getAllSamplesFromVCFHeaders(getToolkit());
        }
        else {
            samples = new HashSet<String>();
            samples.add(sample);
        }
        this.hapTracker = new HaplotypeTracker(samples, requirePQ, getToolkit(), out);
    }

    public static Set<String> getAllSamplesFromVCFHeaders(GenomeAnalysisEngine toolkit) {
        Set<String> samples = new HashSet<String>();
        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(toolkit, null);
        for (VCFHeader header : rodNameToHeader.values()) {
            for (String sample : header.getGenotypeSamples())
                samples.add(sample);
        }

        return samples;
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

        hapTracker.trackSite(tracker.getValues(variantCollection.variants, context.getLocation()), ref);

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
        hapTracker.finalizeAllHaplotypes();
        System.out.println("map was called " + result + " times.");
    }
}