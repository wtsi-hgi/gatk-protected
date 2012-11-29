package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.ConsensusSummarizer;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.NA12878DBArgumentCollection;

import java.io.PrintStream;

public class SummarizeConsensus extends NA12878DBWalker {
    @Output(doc="Output summary here")
    public PrintStream out;

    @Argument(required=false)
    public boolean detailed = false;

    @Override
    public NA12878DBArgumentCollection.DBType getDefaultDB() {
        return NA12878DBArgumentCollection.DBType.PRODUCTION;
    }

    @Override public boolean isDone() { return true; }

    public void onTraversalDone(Integer result) {
        final ConsensusSummarizer summarizer = new ConsensusSummarizer();
        for ( final MongoVariantContext mvc : db.getConsensusSites(makeSiteSelector())) {
            summarizer.add(mvc);
        }

        summarizer.summaryGATKReport(detailed).print(out);
        super.onTraversalDone(result);
    }
}