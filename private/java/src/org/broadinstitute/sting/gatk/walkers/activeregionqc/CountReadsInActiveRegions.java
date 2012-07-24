package org.broadinstitute.sting.gatk.walkers.activeregionqc;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionExtension;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 12/8/11
 */

@ActiveRegionExtension(extension=50)
public class CountReadsInActiveRegions extends ActiveRegionWalker<CountReadsInActiveRegions.Datum, GATKReport> {
    @Output
    PrintStream out;

    public static class Datum {
        private final GenomeLoc activeRegionLoc;
        private final GenomeLoc extendedLoc;
        public final boolean isActive;
        public int nReads;

        public Datum(final GenomeLoc activeRegionLoc, final GenomeLoc extendedLoc, final boolean active, final int nReads) {
            this.activeRegionLoc = activeRegionLoc;
            this.extendedLoc = extendedLoc;
            isActive = active;
            this.nReads = nReads;
        }
    }

    boolean coinFlip = false;

    @Override
    public double isActive( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {
        if( GenomeAnalysisEngine.getRandomGenerator().nextDouble() > 0.995 ) {
            coinFlip = !coinFlip;
        }
        return ( coinFlip ? 0.999 : 0.0 );
    }

    @Override
    public Datum map( final org.broadinstitute.sting.utils.activeregion.ActiveRegion activeRegion, final RefMetaDataTracker tracker ) {
        return new Datum(activeRegion.getLocation(), activeRegion.getExtendedLoc(), activeRegion.isActive, activeRegion.size());
    }

    @Override
    public GATKReport reduceInit() {
        return GATKReport.newSimpleReport("CountReadsInActiveRegions", "loc", "extended.loc", "is.active", "n.reads");
    }

    @Override
    public GATKReport reduce( final Datum value, final GATKReport report ) {
        report.addRow(value.activeRegionLoc.toString(), value.extendedLoc.toString(), value.isActive, value.nReads);
        return report;
    }

    @Override
    public void onTraversalDone(final GATKReport report) {
        report.print(out);
    }
}
