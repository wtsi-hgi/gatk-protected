package org.broadinstitute.sting.gatk.walkers.activeregionqc;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionExtension;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 12/8/11
 */

@ActiveRegionExtension(extension=50)
public class CountReadsInActiveRegions extends ActiveRegionWalker<Integer, Integer> {

    boolean coinFlip = false;

    public double isActive( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {
        if( GenomeAnalysisEngine.getRandomGenerator().nextDouble() > 0.9995 ) {
            coinFlip = !coinFlip;
        }
        return ( coinFlip ? 0.9995 : 0.0 );
    }

    public Integer map( final ActiveRegion activeRegion, final RefMetaDataTracker tracker ) {
        return activeRegion.size();
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce( final Integer value, final Integer sum ) {
        return value + sum;
    }
}
