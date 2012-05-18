/*
 * Copyright (c) 2012, The Broad Institute
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
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Gather;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportGatherer;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.*;

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
 *   -g "ReadGroupID1 ReadGroupID2"
 *   -L input.intervals
 * </pre>
 */
@PartitionBy(PartitionType.INTERVAL)
public class CoverageByRG extends LocusWalker<LinkedHashMap<String, Long>, LinkedHashMap<String, Long>> /*implements TreeReducible<LinkedHashMap<String, Long>> */ {


    @Output
    @Gather(GATKReportGatherer.class)
    PrintStream out;

    @Argument(fullName = "groupRGs", shortName = "g", doc = "This parameter will take the Read Groups provided (separated by a space) and sum them up in a new column", required = false)
    private List<String> groups = new LinkedList<String>();

    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();


    GATKReportTable reportTable;

    HashMap<String, String> rgGroups = new HashMap<String, String>();
    List<HashSet<String>> readGroupIds = new LinkedList<HashSet<String>>();

    final String columnInterval = "Interval";
    final String columnGC = "GCContent";
    final String columnIntervalSize = "IntervalSize";
    final String columnVariants = "Variants";


    /*
    @Override
    public LinkedHashMap<String, Long> treeReduce(LinkedHashMap<String, Long> left, LinkedHashMap<String, Long> right) {
        // Add everything in right to left
        for ( String key : right.keySet() ) {
            left.put(key, left.get(key) + right.get(key) );
        }
        return left;
    }
    */


    public void initialize() {
        int groupIndex = 1;
        for (String groupString : groups) {
            String groupID = "G" + groupIndex;
            String[] rgs = groupString.split(" ");             // Decode the read groups in the grouping argument

            for (String rg : rgs)
                rgGroups.put(rg, groupID);                      // Update the hash with all RGs that correspond to this group

            HashSet<String> groupSet = new HashSet<String>();
            groupSet.addAll(Arrays.asList(rgs));
            readGroupIds.add(groupSet);                         // Add this RG group to the list of RGs

            groupIndex++;
        }

        for (SAMReadGroupRecord RG : getToolkit().getSAMFileHeader().getReadGroups()) {
            String readGroupID = RG.getReadGroupId();
            if (!rgGroups.containsKey(readGroupID)) {
                HashSet<String> rgSet = new HashSet<String>();
                rgSet.add(readGroupID);
                readGroupIds.add(rgSet);                       // Add this RG group to the list of RGs

                rgGroups.put(readGroupID, readGroupID);        // Update the hash with all RGs that correspond to this group
            }
        }

        reportTable = new GATKReportTable("CoverageByRG", "A table with the coverage per interval for each read group", 4 + rgGroups.size(), true);        //Sets up our report table columns (by Read Groups + GCcontent)
        reportTable.addColumn(columnInterval);
        reportTable.addColumn(columnGC);
        reportTable.addColumn(columnIntervalSize);
        reportTable.addColumn(columnVariants);
        for ( String rg : rgGroups.values() )
            reportTable.addColumn(rg);
    }

    public boolean isReduceByInterval() {
        // onTraversalDone is called after every interval
        return true;
    }

    @Override
    public LinkedHashMap<String, Long> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        LinkedHashMap<String, Long> output = new LinkedHashMap<String, Long>();

        VariantContext variantContext = tracker.getFirstValue(dbsnp.dbsnp, ref.getLocus());
        output.put(columnVariants, (variantContext == null) ? 0L : 1L);

        byte base = ref.getBase();                                             // Update site GC content for interval
        output.put(columnGC, (base == BaseUtils.G || base == BaseUtils.C) ? 1L : 0L);

        ReadBackedPileup pileup = context.getBasePileup();                     // Update site pileup for all groups of read groups 
        for (HashSet<String> rgSet : readGroupIds) {
            ReadBackedPileup rgPileup = pileup.getPileupForReadGroups(rgSet);  // This pileup is null when empty so a check must be added
            String rg1 = (String) rgSet.toArray()[0];                          // only need the first rg to determine which column to add            
            output.put(rgGroups.get(rg1), (rgPileup == null) ? 0L : (long) rgPileup.depthOfCoverage());
        }

        return output;
    }

    @Override
    public LinkedHashMap<String, Long> reduceInit() {
        return new LinkedHashMap<String, Long>();

    }

    @Override
    public LinkedHashMap<String, Long> reduce(LinkedHashMap<String, Long> value, LinkedHashMap<String, Long> sum) {
        for (String key : value.keySet()) {
            if (sum.containsKey(key))
                sum.put(key, value.get(key) + sum.get(key));  // sum all keys (if they exist in the sum object).
            else
                sum.put(key, value.get(key));                 // if they don't, just put the first one.

        }

        return sum;
    }

    @Override
    public void onTraversalDone(List<Pair<GenomeLoc, LinkedHashMap<String, Long>>> results) {
        for (Pair<GenomeLoc, LinkedHashMap<String, Long>> intervalPair : results) {
            GenomeLoc interval = intervalPair.getFirst();
            LinkedHashMap<String, Long> counts = intervalPair.getSecond();

            // Get coverage by taking total counts and dividing by interval length
            for (String key : counts.keySet()) {
                if (!key.equals(columnVariants))
                    reportTable.set(interval.toString(), key, (double) counts.get(key) / (double) interval.size());
            }
            reportTable.set(interval.toString(), columnIntervalSize, interval.size());

            if (counts.containsKey(columnVariants))
                reportTable.set(interval.toString(), columnVariants, counts.get(columnVariants));
        }
        GATKReport report = new GATKReport(reportTable);

        report.print(out);
    }
}
