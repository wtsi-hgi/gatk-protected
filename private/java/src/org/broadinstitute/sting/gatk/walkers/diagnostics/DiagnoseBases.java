package org.broadinstitute.sting.gatk.walkers.diagnostics;

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

/**
 * Checks the error rate of the bases in a BAM file
 *
 * <p>
 *  Traverse all the reads in the bam file keeping track of the mismatches, insertions and deletions and outputs the error rate of the 
 *  sequences. It skips sites present in a gold standard callset (provided as a rod). Future versions should check for the correct 
 *  genotype on these sites and provide a separate output.
 * </p>
 *
 * <h2>Input</h2>
 * <p>
 *  A BAM file and an optional gold standard callset (dbSNP)
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 *  The mismatch, insertion and deletion error rate in this file.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T DiagnoseBases
 *      -R reference.fasta
 *      -I myFile.bam 
 *      -G goldstandard.vcf
 *      -o table
 *  </pre>
 *
 * @author Mauricio Carneiro
 * @since 2/27/12
 */
public class DiagnoseBases extends LocusWalker<DiagnoseBases.ErrorCounts, DiagnoseBases.ErrorCounts> {    
    @Output
    PrintStream out;

    @Input(shortName = "G", fullName = "gold_standard", doc = "gold standard callset of sites to skip", required = false)
    public RodBinding<VariantContext> goldStandardCallset = null;

    @Override
    public boolean includeReadsWithDeletionAtLoci() {
        return true;
    }

    @Override
    public ErrorCounts map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        ErrorCounts result = new ErrorCounts();

        VariantContext vc = (goldStandardCallset == null) ? null : tracker.getFirstValue(goldStandardCallset);

        if (vc == null || vc.isFiltered()) {
            byte refBase = ref.getBase();
            for (PileupElement p : context.getBasePileup()) {
                if (p.isDeletion())
                    result.incDeletions(p.getOffset());

                if (p.getBase() != refBase)
                    result.incMismatches(p.getOffset());

                if (p.isBeforeInsertion())
                    result.incInsertions(p.getOffset());

                result.incTotalBases(p.getOffset());
            }
        }
        return result;
    }

    @Override
    public ErrorCounts reduceInit() {
        return new ErrorCounts();
    }

    @Override
    public ErrorCounts reduce(ErrorCounts value, ErrorCounts sum) {
        return sum.add(value);
    }

    @Override
    public void onTraversalDone(ErrorCounts result) {
        result.report();
    }

    public class ErrorCounts {
        public long mismatches;
        public long insertions;
        public long deletions;
        public long totalBases;

        public HashMap<Integer, Long> mismatchesByPosition = new HashMap<Integer, Long>(250);
        public HashMap<Integer, Long> insertionsByPosition = new HashMap<Integer, Long>(250);
        public HashMap<Integer, Long> deletionsByPosition  = new HashMap<Integer, Long>(250);
        public HashMap<Integer, Long> totalBasesByPosition = new HashMap<Integer, Long>(250);

        public void incMismatches(int offset) {
            this.mismatches++;
            incMap(this.mismatchesByPosition, offset);
        }

        public void incInsertions(int offset) {
            this.insertions++;
            incMap(this.insertionsByPosition, offset);
        }

        public void incDeletions(int offset) {
            this.deletions++;
            incMap(this.deletionsByPosition, offset);
        }

        public void incTotalBases(int offset) {
            this.totalBases++;
            incMap(this.totalBasesByPosition, offset);
        }
        
        public ErrorCounts add (ErrorCounts other) {
            this.mismatches += other.mismatches;
            this.deletions += other.deletions;
            this.insertions += other.insertions;
            this.totalBases += other.totalBases;
            copyMap(this.mismatchesByPosition, other.mismatchesByPosition);
            copyMap(this.insertionsByPosition, other.insertionsByPosition);
            copyMap(this.deletionsByPosition,  other.deletionsByPosition);
            copyMap(this.totalBasesByPosition, other.totalBasesByPosition);
            return this;
        }
        
        public void report() {
            System.out.println("Error rates:\n");
            System.out.println(String.format("\tMismatch rate: %.2f%%\n\tInsertion rate: %.2f%%\n\tDeletion rate: %.2f%%\n", (double) 100*mismatches/totalBases, (double) 100*insertions/totalBases, (double) 100*deletions/totalBases));
            out.println("ReadPosition\tMismatches\tInsertions\tDeletions\tTotalBases");
            for (int position : totalBasesByPosition.keySet()) {
                out.println(String.format("%d\t%d\t%d\t%d\t%d",
                        position,
                        mismatchesByPosition.containsKey(position) ? mismatchesByPosition.get(position): 0,
                        insertionsByPosition.containsKey(position) ? insertionsByPosition.get(position): 0,
                        deletionsByPosition.containsKey(position) ?  deletionsByPosition.get(position) : 0,
                        totalBasesByPosition.get(position)));
            }
        }
        
    }

    private static void incMap(HashMap<Integer, Long> map, int location) {
        map.put(location, map.containsKey(location) ? map.get(location) + 1 : 1);
    }

    private static void copyMap(HashMap<Integer, Long> dest, HashMap<Integer, Long> source) {
        for (Map.Entry<Integer, Long> v : source.entrySet()) {
            int key = v.getKey();
            long value = v.getValue();
            long total = dest.containsKey(key) ? dest.get(key) + value : value;
            dest.put(key, total);
        }
    }
}
