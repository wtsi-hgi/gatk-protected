package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.*;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 8/3/11
 * Time: 2:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class SlidingRead {
    private SAMRecord read;

    public SlidingRead(SAMRecord read) {
        this.read = read;
    }

    public SAMRecord getRead() {
        return read;
    }

    public SAMRecord trimToVariableRegion(int refStart, int refStop) {
        int start = read.getAlignmentStart();
        int stop = read.getAlignmentEnd();
        SAMRecord clippedRead = read;

        ReadClipper clipper = new ReadClipper(read);

        // check to see if read is contained in region
        if ( start < refStop && stop > refStart) {
            if ( start < refStart && stop > refStop )
                return clipper.hardClipBothEndsByReferenceCoordinates(refStart-1, refStop+1);
            if ( start < refStart )
                clippedRead = clipper.hardClipByReferenceCoordinates(-1, refStart-1);
            if ( stop > refStop )
                clippedRead = clipper.hardClipByReferenceCoordinates(refStop+1, -1);
            return clippedRead;
        }
        else
            return new SAMRecord(read.getHeader());
    }

    public SlidingRead clipStart(int refStop) {
        if (refStop >= read.getAlignmentEnd())
            return null;
        ReadClipper readClipper = new ReadClipper(read);
        return new SlidingRead(readClipper.hardClipByReferenceCoordinates(-1, refStop));
    }


}

