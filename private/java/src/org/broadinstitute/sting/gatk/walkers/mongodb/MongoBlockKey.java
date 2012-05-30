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
 * Block names are contig|location_modulus
 * e.g. 1:14354654 -> 1|143, Y:12345 -> Y|0
 */
final public class MongoBlockKey {
    private final String key;

    public String getKey() {
        return key;
    }

    // maximum number of bases in a block
    private final int BLOCK_SIZE = 100000;

    public MongoBlockKey(VariantContext vc) {
        key = generateKey(vc.getChr(), vc.getStart());
    }

    public MongoBlockKey(GenomeLoc gl) {
        key = generateKey(gl.getContig(), gl.getStart());
    }

    private String generateKey(String contig, int start) {
        return contig + "|" + (start / BLOCK_SIZE);
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