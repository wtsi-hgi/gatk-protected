package org.broadinstitute.sting.gatk.walkers.newassociation.features.old;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.gatk.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/4/11
 * Time: 1:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class MateSameStrand extends BinaryFeatureAggregator {

    public boolean extractFeature(GATKSAMRecord record) {
        return record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag();
    }

    public boolean featureDefined(GATKSAMRecord record) {
        return record.getReadPairedFlag() && ! record.getMateUnmappedFlag() && record.getReferenceIndex().equals(record.getMateReferenceIndex());
    }


    public MateSameStrand(RFAArgumentCollection col) {
        super(col);
    }
}
