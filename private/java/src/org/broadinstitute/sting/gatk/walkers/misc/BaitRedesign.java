package org.broadinstitute.sting.gatk.walkers.misc;

import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.alignment.bwa.BWAConfiguration;
import org.broadinstitute.sting.alignment.bwa.BWTFiles;
import org.broadinstitute.sting.alignment.bwa.c.BWACAligner;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.table.TableFeature;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/3/12
 * Time: 8:40 PM
 * To change this template use File | Settings | File Templates.
 */
@Reference(window=@Window(start=0, stop=1200))
public class BaitRedesign extends RodWalker<Integer,Integer> {

    @Input(shortName="tbl",fullName="table",required=true,doc="The input bait list table")
    RodBinding<TableFeature> baitBinding = null;

    @Output(shortName="seq",fullName="baitSequence",required=true,doc="The file to write output bait sequences to (in table format)")
    PrintStream baitSeq = null;

    @Output(shortName="mtd",fullName="baitMetaData",required=true,doc="The file to write output bait metadata to (in table format)." +
            " Metadata includes the number of increments of %1 the bait could handle before the optimal, and the increments of %dist" +
            "from optimal. For downstream tabulation of the overall table that was sampled.")
    PrintStream baitMeta = null;

    BWACAligner aligner = null;
    private SAMFileHeader header = null;

    public void initialize() {
        // intialize BWA bindings for checking specificity of enabled variation
        File targetReferenceFile = getToolkit().getArguments().referenceFile;
        BWTFiles bwtFiles = new BWTFiles(targetReferenceFile.getAbsolutePath());
        BWAConfiguration configuration = new BWAConfiguration();
        configuration.maximumEditDistance = 0.30f;
        configuration.mismatchPenalty = 1;
        aligner = new BWACAligner(bwtFiles,configuration);
        header = new SAMFileHeader();
        SAMSequenceDictionary referenceDictionary =
                ReferenceSequenceFileFactory.getReferenceSequenceFile(targetReferenceFile).getSequenceDictionary();
        header.setSequenceDictionary(referenceDictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        List<TableFeature> baitFeatures = tracker.getValues(baitBinding);
        if ( baitFeatures.size() == 0 ) {
            return 0; // nothing to do here
        } else {
            // exactly one bait, as per contract
            TableFeature bait = null;
            for ( int i = 0; i < baitFeatures.size(); i++ ) {
                bait = baitFeatures.get(i);
                if ( ref.getLocus().getStart() != bait.getLocation().getStart() ) {
                    bait = null;
                }
            }
            if ( bait == null ) {
                return 0;
            }
            if ( Integer.parseInt(bait.get("Variants")) > 0 ) {
                return 0;
            }
            byte[] baitBases = Arrays.copyOfRange(ref.getBases(), 0, bait.getEnd() - bait.getStart());
            byte[] optimalBases = BaitRedesignUtils.getOptimalBases(aligner,baitBases,bait.getLocation());
            sampleAndPrint(baitBases, optimalBases,bait.getLocation());
        }
        return 1;
    }

    private void sampleAndPrint(byte[] initial, byte[] optimal,GenomeLoc pos) {
        // todo -- implement sampling along the path from intial -> optimal and tabulating metadata for output
        baitSeq.printf("%s:%d-%d\t%.2f\t%.2f\t%.2f\t%s%n",pos.getContig(),pos.getStart(),pos.getStop(),
                BaitRedesignUtils.calculateGC(initial),
                BaitRedesignUtils.calculateGC(optimal),
                BaitRedesignUtils.editDist(initial,optimal),new String(optimal));
    }

    public Integer reduce(Integer map, Integer prevReduce) {
        return map + prevReduce;
    }

}
