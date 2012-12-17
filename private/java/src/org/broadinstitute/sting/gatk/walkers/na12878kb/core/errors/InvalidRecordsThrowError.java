package org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors;

import org.broadinstitute.sting.gatk.walkers.na12878kb.core.MongoVariantContext;

/**
 * Handle errors by simply throwing the exception
 *
 * User: depristo
 * Date: 11/27/12
 * Time: 8:33 AM
 */
public class InvalidRecordsThrowError<T extends MongoVariantContext> implements InvalidRecordHandler<T> {
    @Override
    public void handleFailedRecord(T record, MongoVariantContextException e) {
        throw e;
    }
}
