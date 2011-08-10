package org.broadinstitute.sting.gatk.walkers.phasing.Haplotypes;

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

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.phasing.ReadBackedPhasingWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class HaplotypeTracker {
    private Map<String, Haplotype> waitingHaplotypes = null;
    private boolean requirePQ;
    private GenomeAnalysisEngine toolKit = null;

    public HaplotypeTracker(Collection<String> samples, boolean requirePQ, GenomeAnalysisEngine toolKit, PrintStream out) {
        this.waitingHaplotypes = new HashMap<String, Haplotype>();
        for (String sample : samples) {
            waitingHaplotypes.put(sample, null);
        }

        this.requirePQ = requirePQ;
        this.toolKit = toolKit;
        Haplotype.out = out;
    }

    public void trackSite(RefMetaDataTracker tracker, ReferenceContext ref) {
        GenomeLoc curLocus = ref.getLocus();
        finalizeDoneHaplotypes(curLocus);

        int curPosition = curLocus.getStop();
        int prevPosition = curPosition - 1;

        // Extend the haplotypes to include up to this position (BUT EXCLUSIVE OF THIS POSITION):
        for (Map.Entry<String, Haplotype> sampleHapEntry : waitingHaplotypes.entrySet()) {
            Haplotype waitingHaplotype = sampleHapEntry.getValue();

            if (waitingHaplotype == null) {// changed to a new contig:
                // Set the new haplotype to extend from [1, prevPosition]
                if (prevPosition >= 1) {
                    GenomeLoc startInterval = toolKit.getGenomeLocParser().createGenomeLoc(curLocus.getContig(), 1, prevPosition);
                    waitingHaplotype = new Haplotype(startInterval, sampleHapEntry.getKey(), toolKit.getGenomeLocParser());
                    sampleHapEntry.setValue(waitingHaplotype);
                }
            }
            else
                waitingHaplotype.extend(prevPosition);
        }

        Collection<VariantContext> vcs = tracker.getValues(VariantContext.class);
        for (VariantContext vc : vcs) {
            if (vc.isFiltered())
                continue;

            for (Map.Entry<String, Genotype> sampleGtEntry : vc.getGenotypes().entrySet()) {
                String sample = sampleGtEntry.getKey();
                if (waitingHaplotypes.get(sample) == null) // an irrelevant sample
                    continue;
                Genotype gt = sampleGtEntry.getValue();

                if (gt.isHet()) {
                    Haplotype sampleHap = waitingHaplotypes.get(sample);
                    if (sampleHap == null)
                        throw new ReviewedStingException("EVERY sample should have a haplotype [by code above and getToolkit().getSamples()]");

                    // Terminate the haplotype before here:
                    if (!gt.isPhased() || (requirePQ && !gt.hasAttribute(ReadBackedPhasingWalker.PQ_KEY))) {
                        sampleHap.finalizeHaplotype();

                        // Start a new haplotype from the current position:
                        sampleHap = new Haplotype(curLocus, sample, toolKit.getGenomeLocParser());
                        waitingHaplotypes.put(sample, sampleHap);
                    }
                    else {
                        sampleHap.extend(curPosition);
                    }

                    sampleHap.addPhasedHet(gt);
                }
            }
        }
    }

    private void finalizeDoneHaplotypes(GenomeLoc curLocus) {
        for (Map.Entry<String, Haplotype> sampleHapEntry : waitingHaplotypes.entrySet()) {
            Haplotype waitingHaplotype = sampleHapEntry.getValue();

            if (waitingHaplotype != null) {
                if (curLocus == null || !waitingHaplotype.interval.onSameContig(curLocus)) {
                    sampleHapEntry.setValue(null);

                    // Set the output haplotype to terminate at the end of its contig:
                    int contigLength = getContigLength(waitingHaplotype.interval.getContig());
                    waitingHaplotype.extend(contigLength);
                    waitingHaplotype.finalizeHaplotype();
                }
            }
        }
    }

    private int getContigLength(String contig) {
        return toolKit.getGenomeLocParser().getContigInfo(contig).getSequenceLength();
    }

    public void finalizeAllHaplotypes() {
        finalizeDoneHaplotypes(null);
    }
}