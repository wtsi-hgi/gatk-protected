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

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

/**
 * Reports the distribution of quality scores of bases in the BAM file(s).
 * <p>
 * Runs through every read counting the number of times a Q score has occurred.
 * </p>
 * <h2>Input</h2>
 * <p>
 * One or more BAM files.
 * </p>
 * <h2>Output</h2>
 * <p>
 * A table with the counts per PHRED scaled quality score from 0 to MAX_QUAL (defined in picard tools - currently 90)
 * </p>
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T QualityScoreDistribution
 *      -R reference.fasta
 *      -I input.bam
 *      -o output.tbl
 *  </pre>
 *
 * @author Mauricio Carneiro
 * @since 11/22/11
 */
public class QualityScoreDistribution extends ReadWalker<HashMap<Byte, Long>, HashMap<Byte, Long>> {
    @Output
    PrintStream out;


    public HashMap<Byte, Long> reduceInit() {
        HashMap<Byte, Long> qualsMap = new HashMap<Byte, Long>(QualityUtils.MAX_QUAL_SCORE);
        return initQualsMap(qualsMap);
    }

    @Override
    public HashMap<Byte, Long> map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {
        HashMap<Byte, Long> qualsMap = new HashMap<Byte, Long>(QualityUtils.MAX_QUAL_SCORE);
        byte [] quals = read.getBaseQualities();
        initQualsMap(qualsMap);
        for (int i = 0; i < quals.length; i++) {
            byte q = quals[i];
            if (!qualsMap.containsKey(q))
                throw new ReviewedStingException(String.format("Invalid quality score %d on read %s", q, read));

            int count = read.isReducedRead() ? read.getReducedCount(i) : 1;
            qualsMap.put(q, qualsMap.get(q) + count);
        }

        return qualsMap;
    }

    @Override
    public HashMap<Byte, Long> reduce(HashMap<Byte, Long> value, HashMap<Byte, Long> sum) {
        return combineHashes(value, sum);
    }

    @Override
    public void onTraversalDone(HashMap<Byte, Long> result) {
        super.onTraversalDone(result);

        GATKReportTable table = new GATKReportTable("QualityScoreDistribution", "Distribution of all quality scores found in the input BAM(s)");
        table.addPrimaryKey("Qual", true);
        table.addColumn("Count", true);

        for (Map.Entry<Byte, Long> entry : result.entrySet())
            table.set(entry.getKey(), "Count", entry.getValue());

        GATKReport report = new GATKReport(table);
        report.print(out);
    }

    /**
     * Initializes the HashMap with all Quality Scores from 0 to MAX_QUAL and sets all counts to zero
     *
     * @param map
     * @return it returns the initialized map as a convenience
     */
    private static HashMap<Byte, Long> initQualsMap (HashMap<Byte, Long> map) {
        for (byte i=0; i<=QualityUtils.MAX_QUAL_SCORE; i++)
            map.put(i, 0L);
        return map;
    }

    /**
     * Adds up the two hashes.
     *
     * it assumes that both hashes have exactly the same set of keys (which should be okay because they are all initialized
     * the same way)
     *
     * @param map1
     * @param map2
     * @return the sume of the two hashes
     */
    private static HashMap<Byte, Long> combineHashes (HashMap<Byte, Long> map1, HashMap<Byte, Long> map2) {
        for (Map.Entry<Byte, Long> entry : map1.entrySet()) {
            Byte qual = entry.getKey();
            Long entryCount = entry.getValue();
            Long sumCount = map2.get(qual);

            if (sumCount == null)
                sumCount = 0L;

            map2.put(qual, sumCount + entryCount);
        }
        return map2;
    }



}