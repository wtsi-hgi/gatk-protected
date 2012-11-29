package org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors;

import org.broadinstitute.sting.gatk.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

/**
 * Exception specific to MongoVariantContext formatting issues
 * User: depristo
 * Date: 11/21/12
 * Time: 9:45 AM
 */
public class MongoVariantContextException extends ReviewedStingException {
    final MongoVariantContext mongoVariantContext;

    public MongoVariantContextException(String msg, MongoVariantContext mongoVariantContext) {
        super(msg);
        this.mongoVariantContext = mongoVariantContext;
    }

    public MongoVariantContextException(String message, Throwable throwable, MongoVariantContext mongoVariantContext) {
        super(message, throwable);
        this.mongoVariantContext = mongoVariantContext;
    }
}
