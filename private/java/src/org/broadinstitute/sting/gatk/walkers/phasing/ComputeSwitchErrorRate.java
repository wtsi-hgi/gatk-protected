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
        report.addTable("SwitchMetrics", "Specifies metrics regarding the switches");
        GATKReportTable switchMetrics = report.getTable("SwitchMetrics");
        switchMetrics.addPrimaryKey("sample");
        switchMetrics.addColumn("markersSeen", 0);
        switchMetrics.addColumn("numSwitches", 0);
        switchMetrics.addColumn("switchErrorRate", 0);
        switchMetrics.addColumn("switchState", false, false);
    }

    private boolean isSwitched(Genotype a, Genotype b) {
        return !a.sameGenotype(b, false);
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            Collection<VariantContext> evals = tracker.getVariantContexts(ref, "eval", null, ref.getLocus(), true, true);
            Collection<VariantContext> truths = tracker.getVariantContexts(ref, "truth", null, ref.getLocus(), true, true);

            VariantContext eval = evals.iterator().hasNext() ? evals.iterator().next() : null;
            VariantContext truth = truths.iterator().hasNext() ? truths.iterator().next() : null;

            if (eval != null && truth != null) {
                for (Trio trio : trios) {
                    String child = trio.getChild();
                    Genotype evalG = eval.getGenotype(child);
                    Genotype truthG = truth.getGenotype(child);

                    if (!eval.isFiltered() && evalG.isHet() && evalG.isPhased() && !truth.isFiltered() && truthG.isHet() && truthG.isPhased()) {
                        GATKReportTable switchMetrics = report.getTable("SwitchMetrics");

                        switchMetrics.increment(child, "markersSeen");
                        if ((Integer) switchMetrics.get(child, "markersSeen") == 1) {
                            switchMetrics.set(child, "switchState", isSwitched(evalG, truthG));
                        } else {
                            if ((Boolean) switchMetrics.get(child, "switchState") != isSwitched(evalG, truthG)) {
                                switchMetrics.increment(child, "numSwitches");
                                switchMetrics.set(child, "switchState", ! (Boolean) switchMetrics.get(child, "switchState"));
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

        switchMetrics.divideColumns("switchErrorRate", "numSwitches", "markersSeen");

        report.print(out);
    }
}
