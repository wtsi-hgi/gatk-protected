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

package org.broadinstitute.sting.gatk.walkers.diagnostics;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

/**
 * Walks along reference and calculates the GC content for each interval.
 * <p/>
 * <p/>
 * <h2>Input</h2>
 * <p>
 * <ul>
 * <li>A reference file</li>
 * <li>A bam file OR multiple bams</li>
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * Tab-deliminated text file showing average coverage per read group per interval.
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CoverageByRG \
 *   -I input.bam \
 *   -o output.txt \
 *   -L input.intervals
 * </pre>
 * <p/>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CoverageByRG \
 *   -I input1.bam \
 *   -I input2.bam \
 *   -o output.txt \
 *   -L input.intervals
 * </pre>
 */
public class CoverageByRG extends LocusWalker<LinkedHashMap<String, Long>, LinkedHashMap<String, Long>> {

    @Output
    PrintStream out;

    GATKReportTable reportTable;
    List<String> readGroupIds;


    public void initialize() {
        //Sets up our report table columns (by Read Groups + GCcontent)

        readGroupIds = new LinkedList<String>();

        reportTable = new GATKReportTable("CoverageByRG", "A table with the coverage per interval for each read group", true);
        reportTable.addPrimaryKey("Interval", true);
        for (SAMReadGroupRecord RG : getToolkit().getSAMFileHeader().getReadGroups()) {
            readGroupIds.add(RG.getReadGroupId());
            reportTable.addColumn(RG.getReadGroupId(), 0, true);
        }
        reportTable.addColumn("GCContent", 0, true);
    }

    public boolean isReduceByInterval() {
        // onTraversalDone is called after every interval
        return true;
    }

    @Override
    public LinkedHashMap<String, Long> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        ReadBackedPileup pileup = context.getBasePileup();

        LinkedHashMap<String, Long> output = new LinkedHashMap<String, Long>();

        for (String RG : readGroupIds) {
            // This pileup is null when empty so a check must be added
            ReadBackedPileup RGpileup = pileup.getPileupForReadGroup(RG);
            output.put(RG, (RGpileup == null) ? 0L : (long) RGpileup.depthOfCoverage());
        }

        byte base = ref.getBase();
        output.put("GCContent", (base == BaseUtils.G || base == BaseUtils.C) ? 1L : 0L);

        return output;
    }

    @Override
    public LinkedHashMap<String, Long> reduceInit() {
        return new LinkedHashMap<String, Long>();

    }

    @Override
    public LinkedHashMap<String, Long> reduce(LinkedHashMap<String, Long> value, LinkedHashMap<String, Long> sum) {
        if (sum.isEmpty()) {
            return value;
        }

        for (String RG : value.keySet()) {
            // Find the sum with the same RG
            sum.put(RG, value.get(RG) + sum.get(RG));
        }

        return sum;
    }

    @Override
    public void onTraversalDone(List<Pair<GenomeLoc, LinkedHashMap<String, Long>>> results) {
        for (Pair<GenomeLoc, LinkedHashMap<String, Long>> intervalPair : results) {
            GenomeLoc interval = intervalPair.getFirst();
            LinkedHashMap<String, Long> counts = intervalPair.getSecond();

            // Get coverage by taking total counts and diving by interval length
            for (String key : counts.keySet()) {
                reportTable.set(interval.toString(), key, (double) counts.get(key) / (double) interval.size());
            }

        }

        reportTable.write(out);

    }
}
