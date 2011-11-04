package org.broadinstitute.sting.gatk.walkers.reducereads;

import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 8/3/11
 * Time: 2:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class SlidingRead {
    protected GATKSAMRecord read;

    public SlidingRead(GATKSAMRecord read) {
        this.read = read;
    }

    public GATKSAMRecord getRead() {
        return read;
    }

    public GATKSAMRecord trimToVariableRegion(int refStart, int refStop) {
        int start = read.getAlignmentStart();
        int stop = read.getAlignmentEnd();
        GATKSAMRecord clippedRead = read;

        ReadClipper clipper = new ReadClipper(read);

        // check to see if read is contained in region
        if ( start <= refStop && stop >= refStart) {
            if ( start < refStart && stop > refStop )
                clippedRead = clipper.hardClipBothEndsByReferenceCoordinates(refStart-1, refStop+1);
            else if ( start < refStart )
                clippedRead = clipper.hardClipByReferenceCoordinatesLeftTail(refStart-1);
            else if ( stop > refStop )
                clippedRead = clipper.hardClipByReferenceCoordinatesRightTail(refStop+1);
            return clippedRead;
        }
        else
            return new GATKSAMRecord(read.getHeader());
    }

    public SlidingRead clipStart(int refNewStart) {
        if (refNewStart > read.getAlignmentEnd())
            return null;
        if (refNewStart <= read.getAlignmentStart())
            return this;
        ReadClipper readClipper = new ReadClipper(read);
        return new SlidingRead(readClipper.hardClipByReferenceCoordinatesLeftTail(refNewStart-1));
    }


}

