package org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/4/11
 * Time: 12:52 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class BinaryFeatureAggregator extends ReadFeatureAggregator<Boolean> {

    public BinaryFeatureAggregator(RFAArgumentCollection col) {
        super(col);
    }

    public void aggregate(Boolean hasFeature) {
	// now robustified
        mean = ( (hasFeature ? 1 : 0)+nReads*mean)/(++nReads);
        var = mean*(1-mean) + Math.pow(2,1-nReads);
        // todo -- this is a massive bias to the variance. Really we want to estimate the prior for a mixture distribution
        // in order to find the proper variance...
        //var = mean*(1-mean) + Math.pow(nReads,Math.pow(mean,4)-1);
    }
}
