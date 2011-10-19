package org.broadinstitute.sting.gatk.walkers.newassociation;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.newassociation.regiontraversal.TriggeringReadStash;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 10/11/11
 * Time: 6:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadFeatureGenotyper extends ReadWalker<SAMRecord,TriggeringReadStash> {

    @Output
    StingSAMFileWriter out;

    public TriggeringReadStash reduceInit() {
        return new TriggeringReadStash(getToolkit().getSAMFileHeader());
    }

    public SAMRecord map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker tracker) {
        if ( Math.abs(read.getInferredInsertSize()) > 1000 ) {
            read.setAttribute("LI",1);
        }
        return read;
    }

    public TriggeringReadStash reduce(SAMRecord read, TriggeringReadStash readBin) {
        if ( read == null )
            return readBin;
        Iterable<SAMRecord> triggeredContext = readBin.compress(read);
        for ( SAMRecord r : triggeredContext ) {
            logger.debug(String.format("%s: %d",r.getReadName(), r.getAlignmentStart()));
            out.addAlignment(r);
        }
        return readBin;
    }

    public void onTraversalDone(TriggeringReadStash sum) {
        Iterable<SAMRecord> finalContext = sum.close();
        for ( SAMRecord r : finalContext ) {
            out.addAlignment(r);
        }
    }
}
