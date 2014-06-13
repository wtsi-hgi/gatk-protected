package org.broadinstitute.gatk.tools.walkers.siblingibd;

/**
 * Created by cwhelan on 5/6/14.
 *
 * When fed observations loci with their predicted IBD states, emits larger regions
 * across which all IBD predictions are consistent.
 */
public class IBDRegionAccumulator {

    String chr;
    int start;
    int currentMax;
    IBDState state;

    /**
     * Given a locus and IBD state, have we changed state from the previous set of observations?
     *
     * @param chr
     * @param coord
     * @param state
     * @return null if we are still in the same state; an IBDRegion if we have changed state
     */
    public IBDRegion regionChange(final String chr, final int coord, final IBDState state) {
        if (this.chr == null) {
            resetState(chr, coord, state);
            return null;
        }
        final boolean regionChanged = (!chr.equals(this.chr) || state != this.state);
        IBDRegion result = null;
        if (regionChanged) {
            result = new IBDRegion(this.chr, this.start, this.currentMax, this.state);
            resetState(chr, coord, state);
        }
        this.chr = chr;
        this.currentMax = coord;

        return result;
    }

    /**
     * @return the final IBD region in the genome, when there are no more loci to accumulate
     */
    public IBDRegion getFinalRegion() {
        return new IBDRegion(this.chr, this.start, this.currentMax, this.state);
    }

    private void resetState(final String chr, final int coord, final IBDState state) {
        this.chr = chr;
        this.start = coord;
        this.currentMax = coord;
        this.state = state;
    }

    /**
     * Represents a genome interval and its IBD state
     */
    public class IBDRegion {
        String chr;
        int start;
        int end;
        IBDState state;

        public IBDRegion(final String chr, final int start, final int end, final IBDState state) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.state = state;
        }

        @Override
        public String toString() {
            return "IBDRegion{" +
                    "chr='" + chr + '\'' +
                    ", start=" + start +
                    ", end=" + end +
                    ", state=" + state +
                    '}';
        }
    }
}
