package org.broadinstitute.sting.gatk.walkers.activeregionqc;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.HashSet;

/**
 * Debugging walker that reimplements CountReads as a LocusWalker.
 * Not meant to be used by anyone.
 *
 * User: rpoplin
 * Date: 1/19/12
 */

public class CountReadsByLocus extends LocusWalker<Integer, Integer>{

    private final HashSet<GATKSAMRecord> myReads = new HashSet<GATKSAMRecord>();
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        int count = 0;
        for(PileupElement p : context.getBasePileup()) {
            GATKSAMRecord read = p.getRead();
            if(!myReads.contains(read)) {
                count++;
                myReads.add(read);
            }
        }
        return count;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
