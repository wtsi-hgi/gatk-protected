package org.broadinstitute.sting.gatk.walkers.newassociation.features;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.gatk.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/4/11
 * Time: 1:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class MateOtherContig extends BinaryFeatureAggregator {

    public MateOtherContig(RFAArgumentCollection col) {
        super(col);
    }

    public boolean extractFeature(GATKSAMRecord record) {
        return ! record.getReferenceName().equals(record.getMateReferenceName());
    }

    public boolean featureDefined(GATKSAMRecord read) {
        return read.getReadPairedFlag() && ! read.getMateUnmappedFlag();
    }
}
