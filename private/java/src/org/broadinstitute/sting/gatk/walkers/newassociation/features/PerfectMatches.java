package org.broadinstitute.sting.gatk.walkers.newassociation.features;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 6/27/11
 * Time: 10:48 AM
 * To change this template use File | Settings | File Templates.
 */
public class PerfectMatches extends BinaryFeatureAggregator {
    private short editDistanceMax = 1;

    public boolean extractFeature(SAMRecord read) {
        return read.getAttribute("NM").equals(0);
    }

    public boolean featureDefined(SAMRecord read) {
        return ((Integer) read.getAttribute("NM") < editDistanceMax);
    }

    public PerfectMatches(RFAArgumentCollection col) {
        super(col);
    }
}
