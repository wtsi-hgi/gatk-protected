package org.broadinstitute.sting.gatk.walkers.newassociation;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.newassociation.regiontraversal.TriggeringReadStash;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 10/11/11
 * Time: 6:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadFeatureGenotyper extends ReadWalker<GATKSAMRecord,TriggeringReadStash> {

    @Output
    StingSAMFileWriter out;

    public TriggeringReadStash reduceInit() {
        return new TriggeringReadStash(getToolkit().getSAMFileHeader());
    }

    public GATKSAMRecord map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker tracker) {
        if ( Math.abs(read.getInferredInsertSize()) > 1000 ) {
            read.setAttribute("LI",1);
        }
        return read;
    }

    public TriggeringReadStash reduce(GATKSAMRecord read, TriggeringReadStash readBin) {
        if ( read == null )
            return readBin;
        Iterable<GATKSAMRecord> triggeredContext = readBin.compress(read);
        for ( GATKSAMRecord r : triggeredContext ) {
            logger.debug(String.format("%s: %d",r.getReadName(), r.getAlignmentStart()));
            out.addAlignment(r);
        }
        return readBin;
    }

    public void onTraversalDone(TriggeringReadStash sum) {
        Iterable<GATKSAMRecord> finalContext = sum.close();
        for ( GATKSAMRecord r : finalContext ) {
            out.addAlignment(r);
        }
    }
}
