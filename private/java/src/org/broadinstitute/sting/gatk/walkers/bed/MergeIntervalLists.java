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

package org.broadinstitute.sting.gatk.walkers.bed;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalSetRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Merges multiple interval lists.
 */
public class MergeIntervalLists extends RodWalker<Integer, Integer> {

    @Input(fullName="targetIntervals", shortName="targetIntervals", doc="interval files to merge", required=true)
    protected List<IntervalBinding<Feature>> intervals;

    @Output(doc="File to which intervals should be written", required=true)
    protected PrintStream writer = null;

    @Argument(fullName = "intervalSetRule", shortName = "intervalSetRule", doc = "Indicates the set merging approach the interval parser should use to combine the various intervals", required = false)
    public IntervalSetRule intervalSetRule = IntervalSetRule.UNION;

    public void initialize() {

        List<GenomeLoc> allIntervals = new ArrayList<GenomeLoc>();
        for ( IntervalBinding intervalBinding : intervals) {
            List<GenomeLoc> intervals = intervalBinding.getIntervals(getToolkit());
            allIntervals = IntervalUtils.mergeListsBySetOperator(intervals, allIntervals, intervalSetRule);
        }

        final GenomeLocSortedSet sorted = IntervalUtils.sortAndMergeIntervals(getToolkit().getGenomeLocParser(), allIntervals, IntervalMergingRule.ALL);
        for ( GenomeLoc loc : sorted )
            writer.println(loc);
    }

    @Override
    public boolean isDone() {
        return true;
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { return 0; }

    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(Integer counter, Integer sum) { return 0; }

    @Override
    public void onTraversalDone(Integer sum) {}
}
