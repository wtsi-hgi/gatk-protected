package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 14, 2011
 */
public abstract class LocalAssemblyEngine {

    public enum ASSEMBLER {
        SIMPLE_DE_BRUIJN
    }

    protected LocalAssemblyEngine() {
    }

    public abstract ArrayList<Haplotype> runLocalAssembly(ArrayList<GATKSAMRecord> reads, Haplotype refHaplotype, int PRUNE_FACTOR);
}
