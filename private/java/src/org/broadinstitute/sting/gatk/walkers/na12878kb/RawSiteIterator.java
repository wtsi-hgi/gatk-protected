package org.broadinstitute.sting.gatk.walkers.na12878kb;

import com.mongodb.DBCursor;
import net.sf.samtools.util.CloseableIterator;

/**
 * Low-level non-public iterator of MongoVariantContexts read from a MongoDB DBCursor
 * User: depristo
 * Date: 11/27/12
 * Time: 9:45 AM
 */
class RawSiteIterator<T extends MongoVariantContext> implements CloseableIterator<T> {
    final DBCursor cursor;

    public RawSiteIterator(DBCursor cursor) {
        this.cursor = cursor;
    }

    @Override public void close() { cursor.close(); }
    @Override public boolean hasNext() { return cursor.hasNext(); }
    @Override public void remove() { throw new UnsupportedOperationException(); }
    @Override public T next() { return (T)cursor.next(); }
}
