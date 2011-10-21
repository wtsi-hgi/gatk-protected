package org.broadinstitute.sting.gatk.walkers.newassociation.features.old;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/4/11
 * Time: 1:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class MateUnmapped extends BinaryFeatureAggregator {

    public boolean extractFeature(SAMRecord record) {
        return record.getMateUnmappedFlag();
    }

    public boolean featureDefined(SAMRecord record) {
        return record.getReadPairedFlag();
    }

    public MateUnmapped(RFAArgumentCollection col) {
        super(col);
    }
}
