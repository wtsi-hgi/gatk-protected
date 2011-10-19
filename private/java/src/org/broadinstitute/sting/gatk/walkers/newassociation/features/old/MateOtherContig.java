package org.broadinstitute.sting.gatk.walkers.newassociation.features.old;

import net.sf.samtools.SAMRecord;
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

    public boolean extractFeature(SAMRecord record) {
        return ! record.getReferenceName().equals(record.getMateReferenceName());
    }

    public boolean featureDefined(SAMRecord read) {
        return read.getReadPairedFlag() && ! read.getMateUnmappedFlag();
    }
}
