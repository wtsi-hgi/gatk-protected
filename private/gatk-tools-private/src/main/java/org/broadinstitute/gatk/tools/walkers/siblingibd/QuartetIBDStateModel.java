package org.broadinstitute.gatk.tools.walkers.siblingibd;


import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.MendelianViolation;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by cwhelan on 5/27/14.
 *
 * Tracks observations of quartet genotypes across the genome and executes the HMM model when
 * appropriate.
 */
public class QuartetIBDStateModel implements IBDStateModel {

    private final Map<SiblingPair, QuartetIBDStateHMM> quartetIBDStateTrackers = new HashMap<>();
    private final PrintStream spFile;
    private final Integer gqThreshold;

    public QuartetIBDStateModel(final List<SiblingPair> siblingPairs, final Integer gqThreshold, final PrintStream spFile) {
        for (final SiblingPair siblingPair : siblingPairs) {
            quartetIBDStateTrackers.put(siblingPair, new QuartetIBDStateHMM());
        }
        this.gqThreshold = gqThreshold;
        this.spFile = spFile;
    }

    @Override
    public List<VariantIBDState> addSite(final VariantContext vc, final SiblingPair siblingPair, final Genotype sib1Gt, final Genotype sib2Gt) {
        final List<VariantIBDState> result = new ArrayList<>();

        final String parent1 = siblingPair.sib1.getPaternalID();
        final String parent2 = siblingPair.sib1.getMaternalID();

        if (vc.hasGenotype(parent1) && vc.hasGenotype(parent2)) {
            final Genotype p1gt = vc.getGenotype(parent1);
            final Genotype p2gt = vc.getGenotype(parent2);
            if (MendelianViolation.isViolation(p1gt, p2gt, sib1Gt) ||
                    MendelianViolation.isViolation(p1gt, p2gt, sib2Gt)) {
                return result;
            }

            if (gqThreshold <= Math.min(Math.min(sib1Gt.getGQ(), sib2Gt.getGQ()), Math.min(p1gt.getGQ(), p2gt.getGQ()))) {
                final QuartetIBDStateHMM ibdStateHmm = quartetIBDStateTrackers.get(siblingPair);
                final IBDObservation observation = IBDObservation.getIBDStateObservation(p1gt, p2gt, sib1Gt, sib2Gt);
                if (observation == IBDObservation.ZERO_OR_ONE_OR_TWO) {
                    return result;
                }
                if (spFile != null) {
                    spFile.print(vc.getChr() + "\t" + vc.getStart() + "\t" + siblingPair.getName() + "\t" + observation.ordinal() + "\n");
                }

                final Iterator<IBDLocus> iterator = ibdStateHmm.addObservation(vc, observation);
                while (iterator.hasNext()) {
                    final IBDLocus locus = iterator.next();
                    result.add(new VariantIBDState(locus.vc, siblingPair, locus.state, locus.statePosteriors));
                }

            }
        }
        return result;
    }

    @Override
    public List<VariantIBDState> finalizeModel() {
        final List<VariantIBDState> finalValues = new ArrayList<>();
        for (final SiblingPair siblingPair : quartetIBDStateTrackers.keySet()) {
            final QuartetIBDStateHMM quartetIBDStateHMM = quartetIBDStateTrackers.get(siblingPair);
            final Iterator<IBDLocus> iterator = quartetIBDStateHMM.runLastChrom();
            while (iterator.hasNext()) {
                final IBDLocus locus = iterator.next();
                final VariantIBDState variantIBDState = new VariantIBDState(locus.vc, siblingPair, locus.state, locus.statePosteriors);
                finalValues.add(variantIBDState);
            }
        }
        return finalValues;
    }
}
