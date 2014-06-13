package org.broadinstitute.gatk.tools.walkers.siblingibd;

import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class SharedMinorAlleleSlidingWindowTest {

    @Test
    public void testSlidingWindow() throws Exception {
        final List<SharedMinorAllele> sharedMinorAlleles = new ArrayList<>();
        sharedMinorAlleles.add(new SharedMinorAllele(1,1));
        sharedMinorAlleles.add(new SharedMinorAllele(1,2));
        sharedMinorAlleles.add(new SharedMinorAllele(2,2));
        sharedMinorAlleles.add(new SharedMinorAllele(0,2));

        final SlidingWindow<SharedMinorAlleleClass> window = new SlidingWindow<>(3);
        assertNull(window.next(1, sharedMinorAlleles.get(0).sharedAlleleClass()));

        assertNull(window.next(10, sharedMinorAlleles.get(1).sharedAlleleClass()));

        SlidingWindow<SharedMinorAlleleClass>.WindowCenterCount next = window.next(20, sharedMinorAlleles.get(2).sharedAlleleClass());
        assertNotNull(next);
        assertEquals(10, next.start);
        assertEquals(.333, next.getCount(SharedMinorAlleleClass.HETVAR_HETVAR), 0.01);
        assertEquals(.333, next.getCount(SharedMinorAlleleClass.HOMVAR_HETVAR), 0.01);
        assertEquals(.333, next.getCount(SharedMinorAlleleClass.HOMVAR_HOMVAR), 0.01);
        assertEquals(0, next.getCount(SharedMinorAlleleClass.HOMVAR_HOMREF), 0.01);
        assertEquals(0, next.getCount(SharedMinorAlleleClass.HETVAR_HOMREF), 0.01);

        next = window.next(23, sharedMinorAlleles.get(3).sharedAlleleClass());
        assertNotNull(next);
        assertEquals(20, next.start);
        assertEquals(0, next.getCount(SharedMinorAlleleClass.HETVAR_HETVAR), 0.01);
        assertEquals(.333, next.getCount(SharedMinorAlleleClass.HOMVAR_HETVAR), 0.01);
        assertEquals(.333, next.getCount(SharedMinorAlleleClass.HOMVAR_HOMVAR), 0.01);
        assertEquals(.333, next.getCount(SharedMinorAlleleClass.HOMVAR_HOMREF), 0.01);
        assertEquals(0, next.getCount(SharedMinorAlleleClass.HETVAR_HOMREF), 0.01);

    }
}