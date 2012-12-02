package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.apache.log4j.Priority;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.ConsensusSummarizer;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.NA12878DBArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.NewlyAddedSites;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.SiteSelector;
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

    @Argument(fullName = "maxIterations", shortName = "maxIterations", required = false)
    protected int maxIterations = Integer.MAX_VALUE;

    @Argument(fullName = "maxQueriesBeforeFullRebuild", shortName = "maxQueriesBeforeFullRebuild", required = false)
    public int maxQueriesBeforeFullRebuild = 100;

    @Output(fullName = "reviewsFile", shortName = "reviewsFile", required = true)
    public VariantContextWriter reviewsFile;

    @Override public boolean isDone() { return true; }

    @Override
    public NA12878DBArgumentCollection.DBType getDefaultDB() {
        return NA12878DBArgumentCollection.DBType.PRODUCTION;
    }

    public void onTraversalDone(Integer result) {
        final SimpleTimer timeSinceStart = new SimpleTimer().start();
        final long maxTimeNano = getToolkit().getRuntimeLimitInNanoseconds();

        if ( ! dontRebuildConsensus ) {
            logger.info("Rebuilding consensus from scratch...");
            db.clearConsensus();
            final ConsensusSummarizer summary = db.updateConsensus(super.makeSiteSelector());
            logger.info("Updated " + summary.getnSites() + " consensus records");
        }

        final NewlyAddedSites newlyAddedSites = new NewlyAddedSites(db);
        int nIterations = 0;
        while ( (timeSinceStart.getElapsedTimeNano() < maxTimeNano || maxTimeNano == -1) && nIterations++ < maxIterations ) {
            try {
                logger.debug("Running cycle " + nIterations);
                pauseCycle();
                final SiteSelector updatedSites = newlyAddedSites.getNewlyAddedLocations(getToolkit().getGenomeLocParser(), maxQueriesBeforeFullRebuild);
                if ( updatedSites != null ) {
                    logger.info("Updating sites " + updatedSites + "...");
                    final int nUpdated = db.updateConsensus(updatedSites, Priority.INFO).getnSites();
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
        if ( updateFrequency > 0 )
            Thread.sleep(updateFrequency);
    }
}