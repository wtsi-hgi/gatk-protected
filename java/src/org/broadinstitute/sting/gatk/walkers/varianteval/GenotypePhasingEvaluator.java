package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.playground.gatk.walkers.phasing.*;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;
import org.broadinstitute.sting.playground.utils.report.utils.TableType;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.apache.log4j.Logger;

import java.util.*;

/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

@Analysis(name = "Genotype Phasing Evaluation", description = "Evaluates the phasing of genotypes in different tracks")
public class GenotypePhasingEvaluator extends VariantEvaluator {
    protected final static Logger logger = Logger.getLogger(GenotypePhasingEvaluator.class);

    // a mapping from sample to stats
    @DataPoint(name = "samples", description = "the phasing statistics for each sample")
    SamplePhasingStatistics samplePhasingStatistics = null;

    SamplePreviousGenotypes samplePrevGenotypes = null;

    public GenotypePhasingEvaluator(VariantEvalWalker parent) {
        super(parent);
        this.samplePhasingStatistics = new SamplePhasingStatistics(getVEWalker().minPhaseQuality);
        this.samplePrevGenotypes = new SamplePreviousGenotypes();
    }

    public String getName() {
        return "GenotypePhasingEvaluator";
    }

    public int getComparisonOrder() {
        return 2;   // we only need to see pairs of (comp, eval)
    }

    public boolean enabled() {
        return true;
    }

    public String toString() {
        return getName() + ": <table>";
    }

    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariantEvalWalker.EvaluationContext group) {
        Reasons interesting = new Reasons();
        if (ref == null)
            return interesting.toString();
        GenomeLoc curLocus = ref.getLocus();

        logger.debug("update2() locus: " + curLocus);
        logger.debug("comp = " + comp + " eval = " + eval);

        Set<String> allSamples = new HashSet<String>();

        Map<String, Genotype> compSampGenotypes = null;
        if (isRelevantToPhasing(comp)) {
            allSamples.addAll(comp.getSampleNames());
            compSampGenotypes = comp.getGenotypes();
        }

        Map<String, Genotype> evalSampGenotypes = null;
        if (isRelevantToPhasing(eval)) {
            allSamples.addAll(eval.getSampleNames());
            evalSampGenotypes = eval.getGenotypes();
        }

        for (String samp : allSamples) {
            logger.debug("sample = " + samp);

            Genotype compSampGt = null;
            if (compSampGenotypes != null)
                compSampGt = compSampGenotypes.get(samp);

            Genotype evalSampGt = null;
            if (evalSampGenotypes != null)
                evalSampGt = evalSampGenotypes.get(samp);

            if (compSampGt == null || evalSampGt == null) {
                // Having a hom site (or an unphased het site) breaks the phasing for the sample - hence, must reset phasing knowledge for both comp and eval [put a null CompEvalGenotypes]:
                if ((compSampGt != null && !permitsTransitivePhasing(compSampGt)) || (evalSampGt != null && !permitsTransitivePhasing(evalSampGt)))
                    samplePrevGenotypes.put(samp, null);
            }
            else { // Both comp and eval have a non-null Genotype at this site:
                Biallele compBiallele = new Biallele(compSampGt);
                Biallele evalBiallele = new Biallele(evalSampGt);

                boolean breakPhasing = false;
                if (!compSampGt.isHet() || !evalSampGt.isHet())
                    breakPhasing = true;
                else { // both are het
                    boolean topMatchesTopAndBottomMatchesBottom = (topMatchesTop(compBiallele, evalBiallele) && bottomMatchesBottom(compBiallele, evalBiallele));
                    boolean topMatchesBottomAndBottomMatchesTop = (topMatchesBottom(compBiallele, evalBiallele) && bottomMatchesTop(compBiallele, evalBiallele));
                    if (!topMatchesTopAndBottomMatchesBottom && !topMatchesBottomAndBottomMatchesTop)
                        breakPhasing = true; // since the 2 VCFs have different diploid genotypes for this sample
                }

                if (breakPhasing) {
                    samplePrevGenotypes.put(samp, null); // nothing to do for this site, AND must remove any history for the future
                }
                else { // comp and eval have the same het Genotype at this site:
                    CompEvalGenotypes prevCompAndEval = samplePrevGenotypes.get(samp);
                    if (prevCompAndEval != null && !prevCompAndEval.getLocus().onSameContig(curLocus)) // exclude curLocus if it is "phased" relative to a different chromosome
                        prevCompAndEval = null;

                    // Replace the previous with current:
                    samplePrevGenotypes.put(samp, curLocus, compSampGt, evalSampGt);

                    if (prevCompAndEval != null) {
                        logger.debug("Potentially phaseable locus: " + curLocus);
                        PhaseStats ps = samplePhasingStatistics.ensureSampleStats(samp);

                        boolean compSampIsPhased = genotypesArePhasedAboveThreshold(compSampGt);
                        boolean evalSampIsPhased = genotypesArePhasedAboveThreshold(evalSampGt);
                        if (compSampIsPhased || evalSampIsPhased) {
                            if (!evalSampIsPhased) {
                                ps.onlyCompPhased++;
                                interesting.addReason("ONLY_COMP", samp, group, "");
                            }
                            else if (!compSampIsPhased) {
                                ps.onlyEvalPhased++;
                                interesting.addReason("ONLY_EVAL", samp, group, "");
                            }
                            else {                        
                                Biallele prevCompBiallele = new Biallele(prevCompAndEval.getCompGenotpye());
                                Biallele prevEvalBiallele = new Biallele(prevCompAndEval.getEvalGenotype());

                                // Sufficient to check only the top of comp, since we ensured that comp and eval have the same diploid genotypes for this sample:
                                boolean topsMatch = (topMatchesTop(prevCompBiallele, prevEvalBiallele) && topMatchesTop(compBiallele, evalBiallele));
                                boolean topMatchesBottom = (topMatchesBottom(prevCompBiallele, prevEvalBiallele) && topMatchesBottom(compBiallele, evalBiallele));

                                if (topsMatch || topMatchesBottom) {
                                    ps.phasesAgree++;
                                }
                                else {
                                    ps.phasesDisagree++;
                                    logger.debug("SWITCHED locus: " + curLocus);
                                    interesting.addReason("SWITCH", samp, group, toString(prevCompBiallele, compBiallele) + " -> " + toString(prevEvalBiallele, evalBiallele));
                                }
                            }
                        }
                        else {
                            ps.neitherPhased++;
                        }
                    }
                }
            }
        }
        logger.debug("\n" + samplePhasingStatistics + "\n");

        return interesting.toString();
    }

    public static boolean isRelevantToPhasing(VariantContext vc) {
        return (vc != null && !vc.isFiltered());
    }

    public boolean genotypesArePhasedAboveThreshold(Genotype gt) {
        if (!gt.genotypesArePhased())
            return false;

        Object pq = gt.getAttributes().get("PQ");
        return (pq == null || (new Double(pq.toString()) >= getVEWalker().minPhaseQuality));
    }

    public boolean permitsTransitivePhasing(Genotype gt) {
        return (gt != null && gt.isHet() && genotypesArePhasedAboveThreshold(gt)); // only a phased het site lets the phase pass through
    }

    public boolean topMatchesTop(Biallele b1, Biallele b2) {
        return b1.getTopAllele().equals(b2.getTopAllele());
    }

    public boolean topMatchesBottom(Biallele b1, Biallele b2) {
        return b1.getTopAllele().equals(b2.getBottomAllele());
    }

    public boolean bottomMatchesTop(Biallele b1, Biallele b2) {
        return topMatchesBottom(b2, b1);
    }

    public boolean bottomMatchesBottom(Biallele b1, Biallele b2) {
        return b1.getBottomAllele().equals(b2.getBottomAllele());
    }

    public String toString(Biallele prev, Biallele cur) {
        return prev.getTopAllele().getBaseString() + "," + cur.getTopAllele().getBaseString() + "|" + prev.getBottomAllele().getBaseString() + "," + cur.getBottomAllele().getBaseString();
    }

    public void finalizeEvaluation() {
    }

    private static class Reasons {
        private StringBuilder sb;

        public Reasons() {
            sb = new StringBuilder();
        }

        public void addReason(String category, String sample, VariantEvalWalker.EvaluationContext evalGroup, String reason) {
             sb.append(category + "(" + sample + " [" + evalGroup.compTrackName + ", " + evalGroup.evalTrackName + "]): " + reason + ";");
        }

        public String toString() {
            if (sb.length() == 0)
                return null;

            return "reasons=" + sb.toString();
        }
    }
}



class CompEvalGenotypes {
    private GenomeLoc loc;
    private Genotype compGt;
    private Genotype evalGt;

    public CompEvalGenotypes(GenomeLoc loc, Genotype compGt, Genotype evalGt) {
        this.loc = loc;
        this.compGt = compGt;
        this.evalGt = evalGt;
    }

    public GenomeLoc getLocus() {
        return loc;
    }

    public Genotype getCompGenotpye() {
        return compGt;
    }
    public Genotype getEvalGenotype() {
        return evalGt;
    }

    public void setCompGenotype(Genotype compGt) {
        this.compGt = compGt;
    }

    public void setEvalGenotype(Genotype evalGt) {
        this.evalGt = evalGt;
    }
}

class SamplePreviousGenotypes {
    private HashMap<String, CompEvalGenotypes> sampleGenotypes = null;

    public SamplePreviousGenotypes() {
        this.sampleGenotypes = new HashMap<String, CompEvalGenotypes>();
    }

    public CompEvalGenotypes get(String sample) {
        return sampleGenotypes.get(sample);
    }

    public void put(String sample, CompEvalGenotypes compEvalGts) {
        sampleGenotypes.put(sample, compEvalGts);
    }

    public void put(String sample, GenomeLoc locus, Genotype compGt, Genotype evalGt) {
        sampleGenotypes.put(sample, new CompEvalGenotypes(locus, compGt, evalGt));
    }
}

class PhaseStats {
    public int neitherPhased;
    public int onlyCompPhased;
    public int onlyEvalPhased;
    public int phasesAgree;
    public int phasesDisagree;

    public PhaseStats() {
        this.neitherPhased = 0;
        this.onlyCompPhased = 0;
        this.onlyEvalPhased = 0;
        this.phasesAgree = 0;
        this.phasesDisagree = 0;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Neither phased: " + neitherPhased + "\tOnly Comp: " + onlyCompPhased + "\tOnly Eval: " + onlyEvalPhased + "\tSame phase: " + phasesAgree + "\tOpposite phase: " + phasesDisagree);
        return sb.toString();
    }

    public static String[] getFieldNamesArray() {
        return new String[]{"total", "neither", "only_comp", "only_eval", "both", "match", "switch", "switch_rate"};
    }

    public Object getField(int index) {
        switch (index) {
            case (0):
                return (neitherPhased + onlyCompPhased + onlyEvalPhased + phasesAgree + phasesDisagree);
            case (1):
                return neitherPhased;
            case (2):
                return onlyCompPhased;
            case (3):
                return onlyEvalPhased;
            case (4):
                return (phasesAgree + phasesDisagree);
            case (5):
                return phasesAgree;
            case (6):
                return phasesDisagree;
            case (7):
                return ((phasesDisagree == 0) ? 0 : ((double) phasesDisagree) / (phasesAgree + phasesDisagree));
            default:
                return -1;
        }
    }
}

/**
 * a table of sample names to genotype phasing statistics
 */
class SamplePhasingStatistics implements TableType {
    private HashMap<String, PhaseStats> sampleStats = null;
    private double minPhaseQuality;

    public SamplePhasingStatistics(double minPhaseQuality) {
        this.sampleStats = new HashMap<String, PhaseStats>();
        this.minPhaseQuality = minPhaseQuality;
    }

    public PhaseStats ensureSampleStats(String samp) {
        PhaseStats ps = sampleStats.get(samp);
        if (ps == null) {
            ps = new PhaseStats();
            sampleStats.put(samp, ps);
        }
        return ps;
    }

    /**
     * @return one row per sample
     */
    public String[] getRowKeys() {
        return sampleStats.keySet().toArray(new String[sampleStats.size()]);
    }

    /**
     * get the column keys
     *
     * @return a list of objects, in this case strings, that are the column names
     */
    public String[] getColumnKeys() {
        return PhaseStats.getFieldNamesArray();
    }

    public Object getCell(int x, int y) {
        String[] rowKeys = getRowKeys();
        PhaseStats ps = sampleStats.get(rowKeys[x]);
        return ps.getField(y);
    }

    public String getName() {
        return "Sample Phasing Statistics (for PQ >= " + minPhaseQuality + ")";
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<String, PhaseStats> sampPhaseStatsEnt : sampleStats.entrySet()) {
            String sample = sampPhaseStatsEnt.getKey();
            PhaseStats ps = sampPhaseStatsEnt.getValue();

            sb.append(sample + "\t" + ps);
        }
        return sb.toString();
    }
}