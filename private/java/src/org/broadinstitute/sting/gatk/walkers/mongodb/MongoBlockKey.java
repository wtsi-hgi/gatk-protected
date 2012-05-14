package org.broadinstitute.sting.gatk.walkers.mongodb;

import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 5/11/12
 * Time: 1:34 PM
 *
 * First simple implementation: use contigs as blocks, so this is just a String
 */
final public class MongoBlockKey {
    public final String key;

    public MongoBlockKey(VariantContext vc) {
        key = vc.getChr();
    }

    public MongoBlockKey(GenomeLoc gl) {
        key = gl.getContig();
    }

    @Override
    public boolean equals(Object other) {
        if (other == null || !(other instanceof MongoBlockKey)) {
            return false;
        }

        MongoBlockKey thisTypeOther = (MongoBlockKey) other;
        return new EqualsBuilder().append(key, thisTypeOther.key).isEquals();
    }

    @Override
    public int hashCode() {
        return new HashCodeBuilder().append(key).toHashCode();
    }
}