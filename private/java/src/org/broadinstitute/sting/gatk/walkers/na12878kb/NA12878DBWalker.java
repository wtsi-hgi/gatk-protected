package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.NA12878DBArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.NA12878KnowledgeBase;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.SiteSelector;

public abstract class NA12878DBWalker extends RodWalker<Integer, Integer> {
    @ArgumentCollection
    private NA12878DBArgumentCollection dbArgumentCollection = new NA12878DBArgumentCollection();

    protected NA12878KnowledgeBase db;

    public abstract NA12878DBArgumentCollection.DBType getDefaultDB();

    public void initialize() {
        logger.info("Connecting to DB");
        if ( dbArgumentCollection.dbToUse == NA12878DBArgumentCollection.DBType.DEFAULT )
            dbArgumentCollection.dbToUse = getDefaultDB();
        db = new NA12878KnowledgeBase(getToolkit().getGenomeLocParser(), dbArgumentCollection);
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return 0;
    }

    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(Integer value, Integer sum) { return value + sum; }

    @Override
    public void onTraversalDone(Integer result) {
        db.close();
        //db.updateConsensus(makeSiteSelector());
    }

    public SiteSelector makeSiteSelector() {
        final SiteSelector select = new SiteSelector(getToolkit().getGenomeLocParser());

        if ( getToolkit().getIntervals() != null ) {
            select.addIntervals(getToolkit().getIntervals());
        }

        return select;
    }
}