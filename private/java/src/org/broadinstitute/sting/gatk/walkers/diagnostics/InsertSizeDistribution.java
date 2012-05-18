package org.broadinstitute.sting.gatk.walkers.diagnostics;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.HashSet;

/**
 * This tool computes the insert size distributions (from the ISIZE field of the GATKSAMRecord) for each sample and read group in a BAM
 *
 * <h2>Input</h2>
 * <p>
 *     Any number of BAM files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 *     Emits a standard <a href="https://www.broadinstitute.org/gsa/wiki/index.php/GATKReport">GATKReport</a>
 *     with two tables:
 *
 *     <dl>
 *         <dt>InsertSizeDistributionBySample</dt>
 *         <dd>Table of read counts with each insert size, for each sample</dd>

 *         <dt>InsertSizeDistributionByReadGroup</dt>
 *         <dd>Table of read counts with each insert size, for each read group</dd>
 *     </dl>
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *     -T InsertSizeDistribution -I my.bam -o insert_sizes.gatkreport.txt
 * </pre>

 * @author Kiran Garimella and Mark DePristo
 * @since 2010-2011
 */
public class InsertSizeDistribution extends ReadWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    private GATKReport report;

    public void initialize() {
        HashSet<String> samplesSeen = new HashSet<String>();
        HashSet<String> rgsSeen = new HashSet<String>();

        for (SAMReadGroupRecord rg : this.getToolkit().getSAMFileHeader().getReadGroups()) {
            samplesSeen.add(rg.getSample());
            rgsSeen.add(rg.getReadGroupId());
        }

        report = new GATKReport();
        report.addTable("InsertSizeDistributionBySample", "Table of insert size distributions", 1 + samplesSeen.size());
        report.addTable("InsertSizeDistributionByReadGroup", "Table of insert size distributions", 1 + rgsSeen.size());

        final GATKReportTable sampleTable = report.getTable("InsertSizeDistributionBySample");
        final GATKReportTable rgTable = report.getTable("InsertSizeDistributionByReadGroup");

        sampleTable.addColumn("insertSize");
        rgTable.addColumn("insertSize");

        for ( String sample : samplesSeen )
            sampleTable.addColumn(sample);

        for ( String rg : rgsSeen )
            rgTable.addColumn(rg);
    }

    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
        return (read.getReadPairedFlag() && read.getFirstOfPairFlag());
    }

    @Override
    public Integer map(ReferenceContext referenceContext, GATKSAMRecord samRecord, ReadMetaDataTracker readMetaDataTracker) {
        final GATKReportTable sampleTable = report.getTable("InsertSizeDistributionBySample");
        final GATKReportTable rgTable = report.getTable("InsertSizeDistributionByReadGroup");

        final int insert = Math.abs(samRecord.getInferredInsertSize());
        final String rg = samRecord.getReadGroup().getReadGroupId();
        final String sample = samRecord.getReadGroup().getSample();

        sampleTable.increment(insert, sample);
        rgTable.increment(insert, rg);

        return null;
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer integer, Integer integer1) {
        return null;
    }

    public void onTraversalDone(Integer sum) {
        // Write report.
        report.print(out);
    }
}
