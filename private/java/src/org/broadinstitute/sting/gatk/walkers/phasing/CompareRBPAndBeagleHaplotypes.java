package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;

public class CompareRBPAndBeagleHaplotypes extends RodWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @Argument(fullName="sample", shortName="sn", doc="Sample to compare", required=false)
    public String sample;

    private ArrayList<VariantContext> rbpHaplotype = new ArrayList<VariantContext>();
    private ArrayList<VariantContext> beagleHaplotype = new ArrayList<VariantContext>();

    public void printHaplotypes(ArrayList<VariantContext> rbpHaplotype, ArrayList<VariantContext> beagleHaplotype, String sample) {
        out.print("%n   rbp:");
        for (int i = 0; i < rbpHaplotype.size(); i++) {
            Genotype rbpg = rbpHaplotype.get(i).getGenotype(sample);
            Allele allele1 = rbpg.getAllele(0);
            out.print(" " + allele1 + " ");
        }
        out.println();

        for (int i = 0; i < rbpHaplotype.size(); i++) {
            Genotype rbpg = rbpHaplotype.get(i).getGenotype(sample);
            Allele allele2 = rbpg.getAlleles().size() == 1 ? rbpg.getAllele(0) : rbpg.getAllele(1);
            out.print(" " + allele2 + " ");
        }
        out.println();

        out.print("%nbeagle:");
        for (int i = 0; i < beagleHaplotype.size(); i++) {
            Genotype beagleg = beagleHaplotype.get(i).getGenotype(sample);
            Allele allele1 = beagleg.getAllele(0);
            out.print(" " + allele1 + " ");
        }
        out.println();

        for (int i = 0; i < beagleHaplotype.size(); i++) {
            Genotype beagleg = beagleHaplotype.get(i).getGenotype(sample);
            Allele allele2 = beagleg.getAlleles().size() == 1 ? beagleg.getAllele(0) : beagleg.getAllele(1);
            out.print(" " + allele2 + " ");
        }
        out.println();
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            Collection<VariantContext> rbps = tracker.getVariantContexts(ref, "rbp", null, ref.getLocus(), true, true);
            Collection<VariantContext> beagles = tracker.getVariantContexts(ref, "beagle", null, ref.getLocus(), true, true);

            VariantContext rbp = rbps.iterator().hasNext() ? rbps.iterator().next() : null;
            VariantContext beagle = beagles.iterator().hasNext() ? beagles.iterator().next() : null;

            if (rbp != null && beagle != null) {
                Genotype rbpg = rbp.getGenotype(sample);

                if (!rbpg.isPhased()) {
                    printHaplotypes(rbpHaplotype, beagleHaplotype, sample);

                    rbpHaplotype.clear();
                    beagleHaplotype.clear();
                }

                rbpHaplotype.add(rbp);
                beagleHaplotype.add(beagle);
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
}
