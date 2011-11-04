package org.broadinstitute.sting.gatk.walkers.newassociation.features;

import org.broadinstitute.sting.gatk.walkers.newassociation.RFAArgumentCollection;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/4/11
 * Time: 1:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class AberrantInsertSize extends BinaryFeatureAggregator {

    private int min;
    private int max;

    public AberrantInsertSize(RFAArgumentCollection col) {
        super(col);
        min = col.lowInsertSize;
        max = col.highInsertSize;
    }

    public boolean extractFeature(GATKSAMRecord rec) {
        return rec.getAttribute("AI") != null ? (rec.getAttribute("AI").equals(1)) : (Math.abs(rec.getInferredInsertSize()) > max || Math.abs(rec.getInferredInsertSize()) < min);
    }

    public boolean featureDefined(GATKSAMRecord rec) {
        return rec.getReadPairedFlag() && rec.getProperPairFlag();
    }
}
