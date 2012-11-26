package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.apache.log4j.Priority;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;

/**
 * Run a server process that continually watches the NA12878 db and updates consensus as needed
 */
public class NA12878KnowledgeBaseServer extends NA12878DBWalker {
    @Argument(fullName = "dontRebuildConsensus", shortName = "dontRebuildConsensus", required=false)
    public boolean dontRebuildConsensus = false;

    @Argument(fullName = "updateFrequency", shortName = "updateFrequency", required=false)
    public int updateFrequency = 10000;

    @Argument(fullName = "maxQueriesBeforeFullRebuild", shortName = "maxQueriesBeforeFullRebuild", required = false)
    public int maxQueriesBeforeFullRebuild = 100;

    @Output(fullName = "reviewsFile", shortName = "reviewsFile", required = true)
    public VariantContextWriter reviewsFile;

    /** For testing only */
    protected boolean breakLoop = false;
    /** For testing only */
    protected Object waitForMe = null;
    /** For testing only */
    protected int maxIterations = Integer.MAX_VALUE;

    @Override public boolean isDone() { return true; }

    public void onTraversalDone(Integer result) {
        final SimpleTimer timeSinceStart = new SimpleTimer().start();
        final long maxTimeNano = getToolkit().getRuntimeLimitInNanoseconds();

        if ( ! dontRebuildConsensus ) {
            logger.info("Rebuilding consensus from scratch...");
            db.clearConsensus();
            final int nUpdated = db.updateConsensus(super.makeSiteSelector());
            logger.info("Updated " + nUpdated + " consensus records");
        }

        final NewlyAddedSites newlyAddedSites = new NewlyAddedSites(db);
        int nIterations = 0;
        while ( timeSinceStart.getElapsedTimeNano() < maxTimeNano && ! breakLoop && nIterations++ < maxIterations ) {
            try {
                pauseCycle();
                final SiteSelector updatedSites = newlyAddedSites.getNewlyAddedLocations(getToolkit().getGenomeLocParser(), maxQueriesBeforeFullRebuild);
                if ( updatedSites != null ) {
                    logger.info("Updating sites " + updatedSites + "...");
                    final int nUpdated = db.updateConsensus(updatedSites, Priority.INFO);
                    logger.info("Updated " + nUpdated + " sites");
                }
            } catch ( InterruptedException e ) {
                logger.info("Interrupted, exiting" + e);
                break;
            }
        }

        if ( reviewsFile != null )
            db.writeReviews(reviewsFile, makeSiteSelector());
        logger.info("Server exiting");

        super.onTraversalDone(result);
    }

    private void pauseCycle() throws InterruptedException {
        if ( waitForMe != null )
            waitForMe.wait();
        else if ( updateFrequency > 0 )
            Thread.sleep(updateFrequency);
    }
}