package org.broadinstitute.gatk.tools.walkers.siblingibd;

import htsjdk.variant.variantcontext.VariantContext;
import org.mockito.Mockito;
import org.testng.annotations.Test;

import java.util.Iterator;

import static org.junit.Assert.*;

public class QuartetIBDStateHMMTest {

    @Test
    public void testRunModel() {
        final QuartetIBDStateHMM t = new QuartetIBDStateHMM(1.0 / 5.0, 1.0 / 100000000000000000.0);
        t.addObservation(getVC("1", 1), IBDObservation.ONE);
        t.addObservation(getVC("1", 2), IBDObservation.ONE);
        t.addObservation(getVC("1", 3), IBDObservation.ONE_OR_TWO);
        t.addObservation(getVC("1", 4), IBDObservation.ONE);
        t.addObservation(getVC("1", 5), IBDObservation.TWO);
        t.addObservation(getVC("1", 6), IBDObservation.TWO);
        t.addObservation(getVC("1", 7), IBDObservation.TWO);
        t.addObservation(getVC("1", 8), IBDObservation.ONE_OR_TWO);
        t.addObservation(getVC("1", 9), IBDObservation.TWO);
        t.addObservation(getVC("1", 10), IBDObservation.ZERO_OR_ONE_OR_TWO);
        final IBDState[] states = t.runModel().viterbiStates;
        for (int i = 0; i < 4; i++) {
            assertEquals(IBDState.ONE, states[i]);
        }
        for (int i = 4; i < 9; i++) {
            assertEquals(IBDState.TWO, states[i]);
        }

    }

    private VariantContext getVC(final String chr, final int coord) {
        final VariantContext vc = Mockito.mock(VariantContext.class);
        Mockito.when(vc.getChr()).thenReturn(chr);
        Mockito.when(vc.getStart()).thenReturn(coord);
        return vc;
    }

    @Test
    public void testTracking() {
        final QuartetIBDStateHMM t = new QuartetIBDStateHMM(1.0 / 100.0, 1.0 / 100000000000000000.0);
        assertTrue(!t.addObservation(getVC("1", 1), IBDObservation.ZERO_OR_ONE).hasNext());
        assertTrue(!t.addObservation(getVC("1", 2), IBDObservation.ZERO).hasNext());
        assertTrue(!t.addObservation(getVC("1", 3), IBDObservation.ONE_OR_TWO).hasNext());
        assertTrue(!t.addObservation(getVC("1", 4), IBDObservation.ZERO).hasNext());
        assertTrue(!t.addObservation(getVC("1", 5), IBDObservation.ZERO).hasNext());
        assertTrue(!t.addObservation(getVC("1", 6), IBDObservation.ONE_OR_TWO).hasNext());
        assertTrue(!t.addObservation(getVC("1", 7), IBDObservation.ONE).hasNext());
        assertTrue(!t.addObservation(getVC("1", 8), IBDObservation.ONE_OR_TWO).hasNext());
        assertTrue(!t.addObservation(getVC("1", 9), IBDObservation.ONE).hasNext());
        assertTrue(!t.addObservation(getVC("1", 10), IBDObservation.ZERO_OR_ONE_OR_TWO).hasNext());
        Iterator<IBDLocus> iterator = t.addObservation(getVC("2", 1), IBDObservation.ONE);
        assertTrue(iterator.hasNext());

        for (int i = 1; i <= 5; i++) {
            final IBDLocus l = iterator.next();
            assertEquals("1", l.vc.getChr());
            assertEquals(i, l.vc.getStart());
            assertEquals(IBDState.ZERO, l.state);
        }

        for (int i = 6; i <= 10; i++) {
            final IBDLocus l = iterator.next();
            assertEquals("1", l.vc.getChr());
            assertEquals(i, l.vc.getStart());
            assertEquals(IBDState.ONE, l.state);
        }
        assertTrue(!iterator.hasNext());
        iterator = t.addObservation(getVC("2", 2), IBDObservation.ONE);
        assertTrue(!iterator.hasNext());
        iterator = t.runLastChrom();

        for (int i = 1; i <= 2; i++) {
            final IBDLocus l = iterator.next();
            assertEquals("2", l.vc.getChr());
            assertEquals(i, l.vc.getStart());
            assertEquals(IBDState.ONE, l.state);
        }

    }

}