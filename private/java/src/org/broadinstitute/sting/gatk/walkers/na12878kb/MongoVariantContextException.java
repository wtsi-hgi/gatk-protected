package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 11/21/12
 * Time: 9:45 AM
 * To change this template use File | Settings | File Templates.
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
