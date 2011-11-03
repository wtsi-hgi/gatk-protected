package org.broadinstitute.sting.gatk.walkers.newassociation.features;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.gatk.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 9/28/11
 * Time: 11:54 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReadFeature {

    public ReadFeature(RFAArgumentCollection collection) { }

    public abstract String getName();

    public abstract String getKey();

    public abstract String getDescription();

    public abstract Object getFeature(GATKSAMRecord read);

    public abstract boolean isDefinedFor(GATKSAMRecord read);
}
