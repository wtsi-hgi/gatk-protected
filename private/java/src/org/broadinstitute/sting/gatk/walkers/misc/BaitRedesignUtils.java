package org.broadinstitute.sting.gatk.walkers.misc;

import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.alignment.bwa.c.BWACAligner;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/3/12
 * Time: 10:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class BaitRedesignUtils {

    public static final double OPT_GC = 0.47;
    private static final double MIXT_COEF = 0.60;
    private static final double ONE_MINUS_MC = 1.0-MIXT_COEF;
    private static final int POP_SIZE = 100;
    private static final int SELECT_SIZE = 15;
    private static final int GENERATIONS = 12;
    private static final double MUTATION_RATE = 0.05;

    private static double optFunction(byte[] seq, BWACAligner aligner, GenomeLoc initialPosition) {
        // count up possible alignments
        int altMapping = 0;
        boolean found = false;
        for (Alignment[] aliSet : aligner.getAllAlignments(seq)) {
            for ( Alignment a : aliSet ) {
                ++altMapping;
                if (! found && a.getContigIndex() == initialPosition.getContigIndex() &&
                        Math.abs(a.getAlignmentStart()-initialPosition.getStart()) < 10 ) {
                    found = true;
                }
            }
        }
        if ( ! found ) {
            // doesn't align anywhere. Not so good.
            return Double.POSITIVE_INFINITY;
        } else {
            --altMapping; // should be one mapping: the current position
        }

        // todo -- don't need to recalculate GC every time, but can cache it and alter it when anything changes
        return MIXT_COEF*altMapping + ONE_MINUS_MC*Math.abs(OPT_GC-calculateGC(seq));
    }

    private static Comparator<byte[]> getSelectionComparator(final BWACAligner aligner, final GenomeLoc position) {
        return new Comparator<byte[]>() {
            @Override
            public int compare(byte[] bytes, byte[] bytes1) {
                return Double.compare(optFunction(bytes,aligner,position),optFunction(bytes1,aligner,position));
            }
        };
    }

    private static double calculateGC(byte[] seq) {
        int tot = 0;
        int gc = 0;
        for ( byte b : seq ) {
            switch (b) {
                case BaseUtils.C:
                    ++gc;
                case BaseUtils.G:
                    ++gc;
                default:
                    ++tot;
            }
        }

        return ((double) gc)/tot;
    }

    public static byte[] getOptimalBases(BWACAligner aligner, byte[] initialSequence, GenomeLoc position) {
        // this is a quasi-genetic algorithm designed to move the initial sequence towards the optimal GC point
        // optimization function to alignability and GC
        List<byte[]> parents = new ArrayList<byte[]>();
        parents.add(initialSequence);
        int generation = 0;
        do {
            List<byte[]> children = reproduce(parents);
            parents = select(children,aligner,position);
            ++generation;
        } while ( generation < GENERATIONS );

        return parents.get(0);
    }

    private static List<byte[]> select(List<byte[]> population, BWACAligner aligner, GenomeLoc position) {
        TreeSet<byte[]> selectionSet = new TreeSet<byte[]>(getSelectionComparator(aligner,position));
        selectionSet.addAll(population);
        List<byte[]> selected = new ArrayList<byte[]>(SELECT_SIZE);
        int num = 0;
        for ( byte[] child : selectionSet ) {
            selected.add(child);
            if ( ++num >= SELECT_SIZE ) {
                break;
            }
        }
        return selected;
    }

    private static List<byte[]> reproduce(List<byte[]> parents) {
        List<byte[]> children = new ArrayList<byte[]>(POP_SIZE);
        // perform recombination to generate a population
        while ( children.size() < POP_SIZE ) {
            byte[] p1 = parents.get(GenomeAnalysisEngine.getRandomGenerator().nextInt(parents.size()));
            byte[] p2 = parents.get(GenomeAnalysisEngine.getRandomGenerator().nextInt(parents.size()));
            if ( p1 == p2 ) {
                children.add(mutate(p1.clone()));
            } else {
                // recombine
                int recomLoc = GenomeAnalysisEngine.getRandomGenerator().nextInt(p1.length);
                byte[] child = p1.clone();
                for ( int i = recomLoc; i < child.length; i++) {
                    child[i] = p2[i];
                }
                // mutate
                children.add(mutate(child));
            }
        }

        return children;
    }

    private static byte[] mutate(byte[] offspring) {
        // every base has a chance to mutate to another base, but biased to increase or decrease GC
        // note that A -> {C,G,T}, so on average we expect GC to go up
        // by symmetry, uniform swapping tends to drive GC -> 50%, and additional bias need not factor in
        for ( int i = 0; i < offspring.length; i++ ) {
            if ( GenomeAnalysisEngine.getRandomGenerator().nextDouble() < MUTATION_RATE ) {
                // mutating this base
                offspring[i] = BaseUtils.baseIndexToSimpleBase(
                        BaseUtils.getRandomBaseIndex(
                                BaseUtils.simpleBaseToBaseIndex(offspring[i])));
            }
        }

        return offspring;
    }
}
