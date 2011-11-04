package org.broadinstitute.sting.gatk.walkers.newassociation.features.old;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.gatk.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 6/27/11
 * Time: 10:35 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReadIndels extends BinaryFeatureAggregator {

    public boolean extractFeature(GATKSAMRecord record) {
        for (CigarElement elem : record.getCigar().getCigarElements() ) {
            if ( elem.getOperator().equals(CigarOperator.INSERTION) || elem.getOperator().equals(CigarOperator.DELETION) ) {
                return true;
            }
        }

        return false;
    }


    public boolean featureDefined(GATKSAMRecord record) {
        return true;
    }

    public ReadIndels(RFAArgumentCollection collection) {
        super(collection);
    }
}
