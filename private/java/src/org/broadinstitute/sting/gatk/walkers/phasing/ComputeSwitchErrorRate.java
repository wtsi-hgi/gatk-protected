package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.Collection;

public class ComputeSwitchErrorRate extends RodWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    private int markersSeen = 0;
    private boolean switchState = false;
    private int numSwitches = 0;
    private int chunkSize = 1;

    private GATKReport report;

    private boolean isSwitched(Genotype a, Genotype b) {
        return !a.sameGenotype(b, false);
    }

    private void toggleSwitchState() {
        switchState = !switchState;
    }

    public void initialize() {
        report = new GATKReport();

        report.addTable("SwitchSites", "Specifies the genomic locations of the haplotype switches");
        report.addTable("SwitchMetrics", "Specifies metrics regarding the switches");
        report.addTable("GenotypeMatches", "Specifies genotypes that match orientation with the truth table");

        GATKReportTable switchSites = report.getTable("SwitchSites");
        switchSites.addPrimaryKey("locus", false);
        switchSites.addColumn("chrom", "unknown");
        switchSites.addColumn("start", "unknown");
        switchSites.addColumn("chunkSize", "unknown");

        GATKReportTable switchMetrics = report.getTable("SwitchMetrics");
        switchMetrics.addPrimaryKey("pk", false);
        switchMetrics.addColumn("markersSeen", "unknown");
        switchMetrics.addColumn("numSwitches", "unknown");
        switchMetrics.addColumn("switchErrorRate", "unknown");

        GATKReportTable genotypePhaseMatches = report.getTable("GenotypeMatches");
        genotypePhaseMatches.addPrimaryKey("locus", false);
        genotypePhaseMatches.addColumn("chrom", "unknown");
        genotypePhaseMatches.addColumn("start", "unknown");
        genotypePhaseMatches.addColumn("evalGenotype", "unknown");
        genotypePhaseMatches.addColumn("compGenotype", "unknown");
        genotypePhaseMatches.addColumn("match", "unknown");
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            Collection<VariantContext> evals = tracker.getVariantContexts(ref, "eval", null, ref.getLocus(), true, true);
            Collection<VariantContext> comps = tracker.getVariantContexts(ref, "comp", null, ref.getLocus(), true, true);

            VariantContext eval = evals.iterator().hasNext() ? evals.iterator().next() : null;
            VariantContext comp = comps.iterator().hasNext() ? comps.iterator().next() : null;

            if (eval != null && comp != null) {
                Genotype evalG = eval.getGenotype("NA12878");
                Genotype compG = comp.getGenotype("NA12878");

                if (!eval.isFiltered() && evalG.isHet() && evalG.isPhased() && !comp.isFiltered() && compG.isHet() && compG.isPhased()) {
                    markersSeen++;

                    GATKReportTable genotypePhaseMatches = report.getTable("GenotypeMatches");
                    genotypePhaseMatches.set(ref.getLocus(), "chrom", ref.getLocus().getContig());
                    genotypePhaseMatches.set(ref.getLocus(), "start", ref.getLocus().getStart());
                    genotypePhaseMatches.set(ref.getLocus(), "evalGenotype", evalG.getGenotypeString());
                    genotypePhaseMatches.set(ref.getLocus(), "compGenotype", compG.getGenotypeString());
                    genotypePhaseMatches.set(ref.getLocus(), "match", evalG.sameGenotype(compG, false));

                    if (markersSeen == 1) {
                        switchState = isSwitched(evalG, compG);
                    } else {
                        if (switchState == isSwitched(evalG, compG)) {
                            chunkSize++;
                        } else {
                            GATKReportTable switchSites = report.getTable("SwitchSites");

                            switchSites.set(ref.getLocus(), "chrom", ref.getLocus().getContig());
                            switchSites.set(ref.getLocus(), "start", ref.getLocus().getStart());
                            switchSites.set(ref.getLocus(), "chunkSize", chunkSize);

                            toggleSwitchState();
                            chunkSize = 1;
                            numSwitches++;
                        }
                    }
                }
            }
        }

        return null;
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    @Override
    public void onTraversalDone(Integer sum) {
        GATKReportTable switchMetrics = report.getTable("SwitchMetrics");

        switchMetrics.set("pk", "markersSeen", markersSeen);
        switchMetrics.set("pk", "numSwitches", numSwitches);
        switchMetrics.set("pk", "switchErrorRate", (double) numSwitches / (double) markersSeen);

        report.print(out);
    }
}
