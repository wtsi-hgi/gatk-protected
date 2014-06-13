package org.broadinstitute.gatk.tools.walkers.siblingibd;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

public class IBDRegionAccumulatorTest
{
    @Test
    public void testRegionAccumulator() throws Exception {
        final IBDRegionAccumulator accumulator = new IBDRegionAccumulator();
        assertNull(accumulator.regionChange("1", 10, IBDState.ONE));
        assertNull(accumulator.regionChange("1", 14, IBDState.ONE));
        assertNull(accumulator.regionChange("1", 20, IBDState.ONE));
        IBDRegionAccumulator.IBDRegion region = accumulator.regionChange("1", 21, IBDState.ZERO);
        assertNotNull(region);
        assertEquals("1", region.chr);
        assertEquals(10, region.start);
        assertEquals(20, region.end);
        assertEquals(IBDState.ONE, region.state);
        assertNull(accumulator.regionChange("1", 22, IBDState.ZERO));
        assertNull(accumulator.regionChange("1", 40, IBDState.ZERO));
        assertNull(accumulator.regionChange("1", 45, IBDState.ZERO));
        region = accumulator.regionChange("2", 10, IBDState.ZERO);
        assertNotNull(region);
        assertEquals("1", region.chr);
        assertEquals(21, region.start);
        assertEquals(45, region.end);
        assertEquals(IBDState.ZERO, region.state);

    }

}