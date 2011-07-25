package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class ComputeSwitchErrorRate extends RodWalker<Integer, Integer> {
    @Argument(shortName="f", fullName="familySpec", required=true, doc="Patterns for the family structure (usage: mom+dad=child).  Specify several trios by supplying this argument many times and/or a file containing many patterns.")
    public ArrayList<String> familySpecs = null;

    @Output
    public PrintStream out;

    private int chunkSize = 1;

    private GATKReport report;

    private class Trio {
        private String mother;
        private String father;
        private String child;

        public Trio(String mother, String father, String child) {
            this.mother = mother;
            this.father = father;
            this.child = child;
        }

        public Trio(String familySpec) {
            String[] pieces = familySpec.split("[\\+\\=]");

            this.mother = pieces[0];
            this.father = pieces[1];
            this.child = pieces[2];
        }

        public String getMother() { return mother; }
        public String getFather() { return father; }
        public String getChild() { return child; }
    }

    private ArrayList<Trio> trios = new ArrayList<Trio>();

    public ArrayList<Trio> getFamilySpecsFromCommandLineInput(ArrayList<String> familySpecs) {
        if (familySpecs != null) {
            // Let's first go through the list and see if we were given any files.  We'll add every entry in the file to our
            // spec list set, and treat the entries as if they had been specified on the command line.
            ArrayList<Trio> specs = new ArrayList<Trio>();
            for (String familySpec : familySpecs) {
                File specFile = new File(familySpec);

                try {
                    XReadLines reader = new XReadLines(specFile);

                    List<String> lines = reader.readLines();
                    for (String line : lines) {
                        specs.add(new Trio(line));
                    }
                } catch (FileNotFoundException e) {
                    specs.add(new Trio(familySpec)); // not a file, so must be a family spec
                }
            }

            return specs;
        }

        return new ArrayList<Trio>();
    }

    public void initialize() {
        trios = getFamilySpecsFromCommandLineInput(familySpecs);

        report = new GATKReport();

//        report.addTable("SwitchSites", "Specifies the genomic locations of the haplotype switches");
//        GATKReportTable switchSites = report.getTable("SwitchSites");
//        switchSites.addPrimaryKey("locus", false);
//        switchSites.addColumn("chrom", "unknown");
//        switchSites.addColumn("start", "unknown");
//        switchSites.addColumn("chunkSize", "unknown");

        report.addTable("SwitchMetrics", "Specifies metrics regarding the switches");
        GATKReportTable switchMetrics = report.getTable("SwitchMetrics");
        switchMetrics.addPrimaryKey("sample");
        switchMetrics.addColumn("markersSeen", "unknown");
        switchMetrics.addColumn("numSwitches", "unknown");
        switchMetrics.addColumn("switchErrorRate", "unknown");
        switchMetrics.addColumn("switchState", false, false);

//        report.addTable("GenotypeMatches", "Specifies genotypes that match orientation with the truth table");
//        GATKReportTable genotypePhaseMatches = report.getTable("GenotypeMatches");
//        genotypePhaseMatches.addPrimaryKey("locus", false);
//        genotypePhaseMatches.addColumn("chrom", "unknown");
//        genotypePhaseMatches.addColumn("start", "unknown");
//        genotypePhaseMatches.addColumn("evalGenotype", "unknown");
//        genotypePhaseMatches.addColumn("compGenotype", "unknown");
//        genotypePhaseMatches.addColumn("match", "unknown");
    }

    private boolean isSwitched(Genotype a, Genotype b) {
        return !a.sameGenotype(b, false);
    }

//    private void toggleSwitchState() {
//        switchState = !switchState;
//    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            Collection<VariantContext> evals = tracker.getVariantContexts(ref, "eval", null, ref.getLocus(), true, true);
            Collection<VariantContext> comps = tracker.getVariantContexts(ref, "comp", null, ref.getLocus(), true, true);

            VariantContext eval = evals.iterator().hasNext() ? evals.iterator().next() : null;
            VariantContext comp = comps.iterator().hasNext() ? comps.iterator().next() : null;

            if (eval != null && comp != null) {
                for (Trio trio : trios) {
                    String child = trio.getChild();
                    Genotype evalG = eval.getGenotype(child);
                    Genotype compG = comp.getGenotype(child);

                    if (!eval.isFiltered() && evalG.isHet() && evalG.isPhased() && !comp.isFiltered() && compG.isHet() && compG.isPhased()) {
                        GATKReportTable switchMetrics = report.getTable("SwitchMetrics");

//                        markersSeen++;
                        switchMetrics.increment(child, "markersSeen");

//                        GATKReportTable genotypePhaseMatches = report.getTable("GenotypeMatches");
//                        genotypePhaseMatches.set(ref.getLocus(), "chrom", ref.getLocus().getContig());
//                        genotypePhaseMatches.set(ref.getLocus(), "start", ref.getLocus().getStart());
//                        genotypePhaseMatches.set(ref.getLocus(), "evalGenotype", evalG.getGenotypeString());
//                        genotypePhaseMatches.set(ref.getLocus(), "compGenotype", compG.getGenotypeString());
//                        genotypePhaseMatches.set(ref.getLocus(), "match", evalG.sameGenotype(compG, false));

                        if ((Integer) switchMetrics.get(child, "markersSeen") == 1) {
//                            switchState = isSwitched(evalG, compG);
                            switchMetrics.set(child, "switchState", isSwitched(evalG, compG));
                        } else {
//                            if (switchState == isSwitched(evalG, compG)) {
                            if ((Boolean) switchMetrics.get(child, "switchState") == isSwitched(evalG, compG)) {
//                                chunkSize++;
//                                switchMetrics.increment(trio.getChild(), "chunkSize");
                            } else {
//                                GATKReportTable switchSites = report.getTable("SwitchSites");
//                                switchSites.set(ref.getLocus(), "chrom", ref.getLocus().getContig());
//                                switchSites.set(ref.getLocus(), "start", ref.getLocus().getStart());
//                                switchSites.set(ref.getLocus(), "chunkSize", chunkSize);

//                                toggleSwitchState();
//                                chunkSize = 1;
//                                numSwitches++;

                                switchMetrics.increment(child, "numSwitches");
                            }
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

        //switchMetrics.set("pk", "markersSeen", markersSeen);
        //switchMetrics.set("pk", "numSwitches", numSwitches);
        //switchMetrics.set("pk", "switchErrorRate", (double) numSwitches / (double) markersSeen);

//        for (Trio trio : trios) {
//            switchMetrics.
//        }

        switchMetrics.divideColumns("switchErrorRate", "numSwitches", "markersSeen");

        report.print(out);
    }
}
