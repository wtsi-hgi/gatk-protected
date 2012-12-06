package org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors;

import com.mongodb.DBCollection;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.MongoVariantContext;

/**
 * Handle errors by logging them and then removing them from the DBCollection itself
 *
 * User: depristo
 * Date: 11/27/12
 * Time: 8:33 AM
 * To change this template use File | Settings | File Templates.
 */
public class InvalidRecordsRemove<T extends MongoVariantContext> extends InvalidRecordsLogError<T> {
    final DBCollection source;

    public InvalidRecordsRemove(final DBCollection source) {
        super();
        if ( source == null ) throw new IllegalArgumentException("source cannot be null");
        this.source = source;
    }

    @Override
    public void handleFailedRecord(T record, MongoVariantContextException e) {
        super.handleFailedRecord(record, e);
        logger.warn("Removing invalid record " + record + " from " + source);
        source.remove(record); // TODO -- should we check that this is actually deleted?
    }
}
