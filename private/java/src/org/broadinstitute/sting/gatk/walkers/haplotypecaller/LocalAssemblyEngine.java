package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 14, 2011
 */
public abstract class LocalAssemblyEngine {

    public enum ASSEMBLER {
        SIMPLE_DE_BRUIJN
    }

    private PrintStream out;
    private IndexedFastaSequenceFile referenceReader;

    protected LocalAssemblyEngine(PrintStream out, IndexedFastaSequenceFile referenceReader) {
        this.out = out;
        this.referenceReader = referenceReader;
    }

    protected PrintStream getOutputStream() { return out; }

    protected IndexedFastaSequenceFile getReferenceReader() { return referenceReader; }

    public abstract ArrayList<Haplotype> runLocalAssembly(ArrayList<GATKSAMRecord> reads, Haplotype refHaplotype);

}
