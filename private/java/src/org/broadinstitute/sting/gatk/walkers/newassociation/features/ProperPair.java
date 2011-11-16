package org.broadinstitute.sting.gatk.walkers.newassociation.features;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.gatk.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 6/8/11
 * Time: 11:59 AM
 * To change this template use File | Settings | File Templates.
 */
public class ProperPair extends BinaryFeatureAggregator {

    public ProperPair(RFAArgumentCollection collection) {
        super(collection);
    }

    public boolean extractFeature(GATKSAMRecord record) {
        return record.getProperPairFlag();
    }

    public boolean featureDefined(GATKSAMRecord record) {
        return record.getReadPairedFlag();
    }
}
