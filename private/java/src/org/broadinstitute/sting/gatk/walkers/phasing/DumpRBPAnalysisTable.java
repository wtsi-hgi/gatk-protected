package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.Collection;

public class DumpRBPAnalysisTable extends RodWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @Argument(fullName="sample", shortName="sn", doc="Sample name to extract", required=true)
    public String SAMPLE;

    private String tableName;
    private GATKReport report;

    public void initialize() {
        tableName = "RBPResults." + SAMPLE;

        report = new GATKReport();
        report.addTable(tableName, "RBP results for " + SAMPLE);

        GATKReportTable rbpTable = report.getTable(tableName);
        rbpTable.addPrimaryKey("pk", false);
        rbpTable.addColumn("sample", SAMPLE);
        rbpTable.addColumn("chr", "unknown");
        rbpTable.addColumn("start", 0);
        rbpTable.addColumn("id", "unknown");
        rbpTable.addColumn("ref", "unknown");
        rbpTable.addColumn("alt", "unknown");

        rbpTable.addColumn("truth.GT", "unknown");
        rbpTable.addColumn("truth.AC", 0.0);
        rbpTable.addColumn("truth.AN", 0.0);
        rbpTable.addColumn("truth.GQ", 0.0);
        rbpTable.addColumn("truth.DP", 0.0);
        rbpTable.addColumn("truth.TP", 0.0);

        rbpTable.addColumn("rbp00.GT", "unknown");
        rbpTable.addColumn("rbp00.AC", 0.0);
        rbpTable.addColumn("rbp00.AN", 0.0);
        rbpTable.addColumn("rbp00.GQ", 0.0);
        rbpTable.addColumn("rbp00.DP", 0.0);
        rbpTable.addColumn("rbp00.PQ", 0.0);

        rbpTable.addColumn("rbp01.GT", "unknown");
        rbpTable.addColumn("rbp01.AC", 0.0);
        rbpTable.addColumn("rbp01.AN", 0.0);
        rbpTable.addColumn("rbp01.AF", 0.0);
        rbpTable.addColumn("rbp01.GQ", 0.0);
        rbpTable.addColumn("rbp01.DP", 0.0);
        rbpTable.addColumn("rbp01.PQ", 0.0);

        rbpTable.addColumn("rbp10.GT", "unknown");
        rbpTable.addColumn("rbp10.AC", 0.0);
        rbpTable.addColumn("rbp10.AN", 0.0);
        rbpTable.addColumn("rbp10.AF", 0.0);
        rbpTable.addColumn("rbp10.GQ", 0.0);
        rbpTable.addColumn("rbp10.DP", 0.0);
        rbpTable.addColumn("rbp10.PQ", 0.0);

        rbpTable.addColumn("rbp11.GT", "unknown");
        rbpTable.addColumn("rbp11.AC", 0.0);
        rbpTable.addColumn("rbp11.AN", 0.0);
        rbpTable.addColumn("rbp11.AF", 0.0);
        rbpTable.addColumn("rbp11.GQ", 0.0);
        rbpTable.addColumn("rbp11.DP", 0.0);
        rbpTable.addColumn("rbp11.PQ", 0.0);

        rbpTable.addColumn("beagle00.GT", "unknown");
        rbpTable.addColumn("beagle00.AF", 0.0);
        rbpTable.addColumn("beagle00.GA", 0.0);
        rbpTable.addColumn("beagle00.AR2", 0.0);
        rbpTable.addColumn("beagle00.DR2", 0.0);

        rbpTable.addColumn("beagle01.GT", "unknown");
        rbpTable.addColumn("beagle01.AF", 0.0);
        rbpTable.addColumn("beagle01.GA", 0.0);
        rbpTable.addColumn("beagle01.AR2", 0.0);
        rbpTable.addColumn("beagle01.DR2", 0.0);

        rbpTable.addColumn("beagle10.GT", "unknown");
        rbpTable.addColumn("beagle10.AF", 0.0);
        rbpTable.addColumn("beagle10.GA", 0.0);
        rbpTable.addColumn("beagle10.AR2", 0.0);
        rbpTable.addColumn("beagle10.DR2", 0.0);

        rbpTable.addColumn("beagle11.GT", "unknown");
        rbpTable.addColumn("beagle11.AF", 0.0);
        rbpTable.addColumn("beagle11.GA", 0.0);
        rbpTable.addColumn("beagle11.AR2", 0.0);
        rbpTable.addColumn("beagle11.DR2", 0.0);
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Collection<VariantContext> truths = tracker.getVariantContexts(ref, "truth", null, ref.getLocus(), true, true);
        VariantContext truth = truths.iterator().hasNext() ? truths.iterator().next() : null;

        Collection<VariantContext> rbp00s = tracker.getVariantContexts(ref, "rbp00s", null, ref.getLocus(), true, true);
        VariantContext rbp00 = rbp00s.iterator().hasNext() ? rbp00s.iterator().next() : null;

        Collection<VariantContext> rbp01s = tracker.getVariantContexts(ref, "rbp01s", null, ref.getLocus(), true, true);
        VariantContext rbp01 = rbp01s.iterator().hasNext() ? rbp01s.iterator().next() : null;

        Collection<VariantContext> rbp10s = tracker.getVariantContexts(ref, "rbp10s", null, ref.getLocus(), true, true);
        VariantContext rbp10 = rbp10s.iterator().hasNext() ? rbp10s.iterator().next() : null;

        Collection<VariantContext> rbp11s = tracker.getVariantContexts(ref, "rbp11s", null, ref.getLocus(), true, true);
        VariantContext rbp11 = rbp11s.iterator().hasNext() ? rbp11s.iterator().next() : null;
        
        Collection<VariantContext> beagle00s = tracker.getVariantContexts(ref, "beagle00s", null, ref.getLocus(), true, true);
        VariantContext beagle00 = beagle00s.iterator().hasNext() ? beagle00s.iterator().next() : null;

        Collection<VariantContext> beagle01s = tracker.getVariantContexts(ref, "beagle01s", null, ref.getLocus(), true, true);
        VariantContext beagle01 = beagle01s.iterator().hasNext() ? beagle01s.iterator().next() : null;

        Collection<VariantContext> beagle10s = tracker.getVariantContexts(ref, "beagle10s", null, ref.getLocus(), true, true);
        VariantContext beagle10 = beagle10s.iterator().hasNext() ? beagle10s.iterator().next() : null;

        Collection<VariantContext> beagle11s = tracker.getVariantContexts(ref, "beagle11s", null, ref.getLocus(), true, true);
        VariantContext beagle11 = beagle11s.iterator().hasNext() ? beagle11s.iterator().next() : null;

        if (truth != null && rbp00 != null && rbp01 != null && rbp10 != null && rbp11 != null && beagle00 != null && beagle01 != null && beagle10 != null && beagle11 != null) {
            GATKReportTable rbpTable = report.getTable(tableName);
            GenomeLoc pk = ref.getLocus();
            rbpTable.set(pk, "chr", ref.getLocus().getContig());
            rbpTable.set(pk, "start", ref.getLocus().getContig());
            rbpTable.set(pk, "id", truth.getID());
            rbpTable.set(pk, "ref", truth.getReference());
            rbpTable.set(pk, "alt", truth.getAltAlleleWithHighestAlleleCount());

            Genotype truthG = truth.getGenotype(SAMPLE);
            rbpTable.addColumn("truth.GT", truthG.getGenotypeString(true));
            rbpTable.addColumn("truth.AC", truth.getAttributeAsInt("AC", 0));
            rbpTable.addColumn("truth.AN", truth.getAttributeAsInt("AN", 0));
            rbpTable.addColumn("truth.AF", truth.getAttributeAsInt("AF", 0));
            rbpTable.addColumn("truth.GQ", truthG.getAttributeAsDouble("GQ", 0.0));
            rbpTable.addColumn("truth.DP", truthG.getAttributeAsInt("DP", 0));
            rbpTable.addColumn("truth.TP", truthG.getAttributeAsDouble("TP", 0.0));

            Genotype rbp00G = rbp00.getGenotype(SAMPLE);
            rbpTable.addColumn("rbp00.GT", rbp00G.getGenotypeString(true));
            rbpTable.addColumn("rbp00.AC", rbp00.getAttributeAsInt("AC", 0));
            rbpTable.addColumn("rbp00.AN", rbp00.getAttributeAsInt("AN", 0));
            rbpTable.addColumn("rbp00.AF", rbp00.getAttributeAsInt("AF", 0));
            rbpTable.addColumn("rbp00.GQ", rbp00G.getAttributeAsDouble("GQ", 0.0));
            rbpTable.addColumn("rbp00.DP", rbp00G.getAttributeAsInt("DP", 0));
            rbpTable.addColumn("rbp00.PQ", rbp00G.getAttributeAsDouble("PQ", 0.0));

            Genotype rbp01G = rbp01.getGenotype(SAMPLE);
            rbpTable.addColumn("rbp01.GT", rbp01G.getGenotypeString(true));
            rbpTable.addColumn("rbp01.AC", rbp01.getAttributeAsInt("AC", 0));
            rbpTable.addColumn("rbp01.AN", rbp01.getAttributeAsInt("AN", 0));
            rbpTable.addColumn("rbp01.AF", rbp01.getAttributeAsInt("AF", 0));
            rbpTable.addColumn("rbp01.GQ", rbp01G.getAttributeAsDouble("GQ", 0.0));
            rbpTable.addColumn("rbp01.DP", rbp01G.getAttributeAsInt("DP", 0));
            rbpTable.addColumn("rbp01.PQ", rbp01G.getAttributeAsDouble("PQ", 0.0));

            Genotype rbp10G = rbp10.getGenotype(SAMPLE);
            rbpTable.addColumn("rbp10.GT", rbp10G.getGenotypeString(true));
            rbpTable.addColumn("rbp10.AC", rbp10.getAttributeAsInt("AC", 0));
            rbpTable.addColumn("rbp10.AN", rbp10.getAttributeAsInt("AN", 0));
            rbpTable.addColumn("rbp10.AF", rbp10.getAttributeAsInt("AF", 0));
            rbpTable.addColumn("rbp10.GQ", rbp10G.getAttributeAsDouble("GQ", 0.0));
            rbpTable.addColumn("rbp10.DP", rbp10G.getAttributeAsInt("DP", 0));
            rbpTable.addColumn("rbp10.PQ", rbp10G.getAttributeAsDouble("PQ", 0.0));

            Genotype rbp11G = rbp11.getGenotype(SAMPLE);
            rbpTable.addColumn("rbp00.GT", rbp11G.getGenotypeString(true));
            rbpTable.addColumn("rbp00.AC", rbp11.getAttributeAsInt("AC", 0));
            rbpTable.addColumn("rbp00.AN", rbp11.getAttributeAsInt("AN", 0));
            rbpTable.addColumn("rbp00.AF", rbp11.getAttributeAsInt("AF", 0));
            rbpTable.addColumn("rbp00.GQ", rbp11G.getAttributeAsDouble("GQ", 0.0));
            rbpTable.addColumn("rbp00.DP", rbp11G.getAttributeAsInt("DP", 0));
            rbpTable.addColumn("rbp00.PQ", rbp11G.getAttributeAsDouble("PQ", 0.0));
            
            Genotype beagle00G = beagle00.getGenotype(SAMPLE);
            rbpTable.addColumn("beagle00.GT", beagle00G.getGenotypeString(true));
            rbpTable.addColumn("beagle00.AF", beagle00.getAttributeAsInt("AF", 0));
            rbpTable.addColumn("beagle00.GA", beagle00G.getAttributeAsDouble("GA", 0.0));
            rbpTable.addColumn("beagle00.AR2", beagle00G.getAttributeAsInt("AR2", 0));
            rbpTable.addColumn("beagle00.DR2", beagle00G.getAttributeAsDouble("DR2", 0.0));

            Genotype beagle01G = beagle01.getGenotype(SAMPLE);
            rbpTable.addColumn("beagle01.GT", beagle01G.getGenotypeString(true));
            rbpTable.addColumn("beagle01.AF", beagle01.getAttributeAsInt("AF", 0));
            rbpTable.addColumn("beagle01.GA", beagle01G.getAttributeAsDouble("GA", 0.0));
            rbpTable.addColumn("beagle01.AR2", beagle01G.getAttributeAsInt("AR2", 0));
            rbpTable.addColumn("beagle01.DR2", beagle01G.getAttributeAsDouble("DR2", 0.0));

            Genotype beagle10G = beagle10.getGenotype(SAMPLE);
            rbpTable.addColumn("beagle10.GT", beagle10G.getGenotypeString(true));
            rbpTable.addColumn("beagle10.AF", beagle10.getAttributeAsInt("AF", 0));
            rbpTable.addColumn("beagle10.GA", beagle10G.getAttributeAsDouble("GA", 0.0));
            rbpTable.addColumn("beagle10.AR2", beagle10G.getAttributeAsInt("AR2", 0));
            rbpTable.addColumn("beagle10.DR2", beagle10G.getAttributeAsDouble("DR2", 0.0));

            Genotype beagle11G = beagle11.getGenotype(SAMPLE);
            rbpTable.addColumn("beagle00.GT", beagle11G.getGenotypeString(true));
            rbpTable.addColumn("beagle00.AF", beagle11.getAttributeAsInt("AF", 0));
            rbpTable.addColumn("beagle00.GQ", beagle11G.getAttributeAsDouble("GQ", 0.0));
            rbpTable.addColumn("beagle00.AR2", beagle11G.getAttributeAsInt("AR2", 0));
            rbpTable.addColumn("beagle00.DR2", beagle11G.getAttributeAsDouble("DR2", 0.0));
        }

        return null;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public Integer reduceInit() {
        return null;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    public void onTraversalDone(Integer sum) {
        report.print(out);
    }
}
