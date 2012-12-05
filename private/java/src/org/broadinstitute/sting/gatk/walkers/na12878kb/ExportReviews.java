package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.NA12878DBArgumentCollection;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;

public class ExportReviews extends NA12878DBWalker {
    @Output
    public VariantContextWriter out;

    @Override public boolean isDone() { return true; }

    @Override
    public void onTraversalDone(Integer result) {
        db.writeReviews(out, super.makeSiteSelector());
        super.onTraversalDone(result);
    }

    @Override
    public NA12878DBArgumentCollection.DBType getDefaultDB() {
        return NA12878DBArgumentCollection.DBType.PRODUCTION;
    }
}