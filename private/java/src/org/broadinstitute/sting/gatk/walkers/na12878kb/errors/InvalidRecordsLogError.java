package org.broadinstitute.sting.gatk.walkers.na12878kb.errors;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.walkers.na12878kb.MongoVariantContext;

/**
 * Handle invalid exceptions simply by logging (and ignoring) them
 *
 * User: depristo
 * Date: 11/27/12
 * Time: 8:33 AM
 */
public class InvalidRecordsLogError<T extends MongoVariantContext> implements InvalidRecordHandler<T> {
    final static Logger logger = Logger.getLogger(InvalidRecordsLogError.class);
    int nBad = 0;

    @Override
    public void handleFailedRecord(T record, MongoVariantContextException e) {
        logger.warn("Invalid record detected " + record + " " + e.getMessage());
        nBad++;
    }

    public int getnBad() {
        return nBad;
    }
}
