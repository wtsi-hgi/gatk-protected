package org.broadinstitute.sting.gatk.walkers.diagnostics;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;

import java.io.PrintStream;

public class ReadLengthDistribution extends ReadWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    private GATKReport report;

    public void initialize() {
        report = new GATKReport();
        report.addTable("ReadLengthDistribution", "Table of read length distributions");
        GATKReportTable table = report.getTable("ReadLengthDistribution");

        table.addPrimaryKey("readLength");

        for (SAMReadGroupRecord rg : this.getToolkit().getSAMFileHeader().getReadGroups()) {
            table.addColumn(rg.getSample(), 0);
        }
    }

    public boolean filter(ReferenceContext ref, SAMRecord read) {
        return (read.getReadPairedFlag() && read.getFirstOfPairFlag());
    }

    @Override
    public Integer map(ReferenceContext referenceContext, SAMRecord samRecord, ReadMetaDataTracker readMetaDataTracker) {
        GATKReportTable table = report.getTable("ReadLengthDistribution");

        int length = Math.abs(samRecord.getReadLength());
        String sample = samRecord.getReadGroup().getSample();

        table.increment(length, sample);

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
        report.print(out);
    }
}
