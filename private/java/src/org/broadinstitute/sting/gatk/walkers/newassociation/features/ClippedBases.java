package org.broadinstitute.sting.gatk.walkers.newassociation.features;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.gatk.walkers.newassociation.RFAArgumentCollection;
import org.broadinstitute.sting.utils.sam.ReadUtils;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 9/28/11
 * Time: 12:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class ClippedBases extends ReadFeature {

    public ClippedBases(RFAArgumentCollection col) {
        super(col);
    }

    public String getName() { return "ClippedBases"; }

    public String getKey() { return "clippedBases"; }

    public String getDescription() { return "the number of clipped bases with Q >= 14"; }

    public Object getFeature(GATKSAMRecord read) {
         int firstClippedToAliStart = read.getUnclippedStart()-read.getAlignmentStart();
        int lastUnclippedToReadEnd = read.getUnclippedEnd()-read.getAlignmentEnd();

        byte[] quals = read.getBaseQualities();
        int nClipped = 0;
        for ( int offset = 0; offset < firstClippedToAliStart; offset++ ) {
            if ( quals[offset] >= 14 ) {
                nClipped++;
            }
        }

        for ( int offset = quals.length - lastUnclippedToReadEnd; offset < quals.length ; offset++ ) {
            if ( quals[offset] >= 14 ) {
                nClipped ++;
            }
        }

        return nClipped;
    }

    public boolean isDefinedFor(GATKSAMRecord read) {
        return ! read.getReadUnmappedFlag();
    }
}
