package org.broadinstitute.sting.gatk.walkers.misc;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 9/7/11
 * Time: 1:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class OriginalQualScoreHistogramWalker extends ReadWalker<Boolean,Long> {

    @Output
    PrintStream out;

    private Map<Integer,long[]> qualsByOffset = new HashMap<Integer,long[]>(256);

    public void initialize() {
        for ( int o = 0; o < 256; o++ ) {
            qualsByOffset.put(o,new long[41]);
        }
    }

    public Long reduceInit() {
        return 0L;
    }

    public boolean filter(SAMRecord read) {
        return read.getAttribute("OQ") != null;
    }

    public Boolean map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {
        byte[] quals = read.getOriginalBaseQualities();
        for ( int offset = 0; offset < quals.length; offset++ ) {
            qualsByOffset.get(offset)[ (int) quals[offset]]++;
        }

        return null;
    }

    public Long reduce(Boolean m, Long r) {
        return r+1;
    }

    public void onTraversalDone(Long total) {
        logger.info("Finished with "+total.toString()+"reads");
        for ( int o = 0; o < 256; o++) {
            StringBuffer stringBuffer = new StringBuffer();
            long[] quals = qualsByOffset.get(o);
            long obs = 0;
            for ( long c : quals ) {
                obs += c;
            }
            if ( obs == 0 ) { break; }
            double prob = ((double)quals[0])/obs;
            stringBuffer.append(String.format("%.3e",prob));
            for ( int qOffset = 1; qOffset < 41; qOffset ++ ) {
                prob = ((double)quals[qOffset])/obs;
                stringBuffer.append(",");
                stringBuffer.append(String.format("%.3e",prob));
            }
            out.printf("%s%n",stringBuffer);
        }
    }
}
