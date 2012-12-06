package org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors;

import org.broadinstitute.sting.gatk.walkers.na12878kb.core.MongoVariantContext;

/**
 * Interface for classes that want to handle their own errors
 *
 * User: depristo
 * Date: 11/27/12
 * Time: 8:27 AM
 */
public interface InvalidRecordHandler<T extends MongoVariantContext> {
    /**
     * Called when an invalid MongoVariantContext object (record) has been detected in the DB
     *
     * @param record the invalid record
     * @param e the exception generated when reading this record
     */
    public void handleFailedRecord(final T record, final MongoVariantContextException e);
}


