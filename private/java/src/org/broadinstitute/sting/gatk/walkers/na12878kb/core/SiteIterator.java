package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import com.mongodb.DBCursor;
import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.InvalidRecordHandler;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.InvalidRecordsThrowError;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.MongoVariantContextException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Generic, functionally rich Java iterator driving a MongoDB cursor containing MongoVariantContext objects
 *
 * User: depristo
 * Date: 11/14/12
 * Time: 10:56 AM
 */
public class SiteIterator<T extends MongoVariantContext> extends PeekableIterator<T> implements Iterable<T>, CloseableIterator<T> {
    final GenomeLocParser parser;

    InvalidRecordHandler<T> errorHandler;

    /** To ensure that records are coming in order */
    GenomeLoc lastLoc = null;

    /**
     * Create a new SiteIterator using parser to make GenomeLocs and reading data from the DBCursor cursor
     * @param parser a GenomeLoc parser creating genome locs, consistent with the data in the cursor
     * @param cursor an initialized DBCursor pointing to data containing values of type T
     */
    public SiteIterator(GenomeLocParser parser, DBCursor cursor) {
        this(parser, cursor, new InvalidRecordsThrowError<T>());
    }

    public SiteIterator(GenomeLocParser parser, DBCursor cursor, final InvalidRecordHandler errorHandler) {
        super(new RawSiteIterator<T>(cursor));

        if ( parser == null ) throw new IllegalArgumentException("Parser cannot be null");
        this.parser = parser;
        setErrorHandler(errorHandler);
    }

    /**
     * Provides this SiteIterator with an InvalidRecordHandler handler to deal with
     * bad records encountered in the DB while iterating.
     *
     * @param errorHandler a non-null error handler
     */
    public void setErrorHandler(InvalidRecordHandler<T> errorHandler) {
        if ( errorHandler == null ) throw new IllegalArgumentException("errorHandler cannot be null");
        this.errorHandler = errorHandler;
    }

    /**
     * Get all of the upcoming T that occur before the contig / start position of loc
     *
     * End of loc is ignored, and must be == start.  Note that by ignoring the start location,
     * we are getting all records with start position < loc.start, including any indels that
     * may span up to (and over) start loc.
     *
     * @param loc the genome loc containing the requested contig and start (stop is ignored)
     * @return a list of all of the records < loc
     */
    @Ensures("result != null")
    public List<T> getSitesBefore(final GenomeLoc loc) {
        if ( loc.size() != 1 ) throw new IllegalArgumentException("GenomeLoc must be a single one bp site but got " + loc);

        final List<T> l = new LinkedList<T>();

        while ( hasNext() ) {
            final T n = peek();
            final GenomeLoc nLoc = n.getLocation(parser);
            if ( nLoc.startsBefore(loc) )
                l.add(next());
            else
                break;
        }

        return l;
    }

    /**
     * Get all of the upcoming T that have the same contig / start position as loc
     *
     * Skips all sites occurring before contig / start

     * @param loc the genome loc containing the requested contig and start (stop is ignored)
     * @return a list of all of the upcoming records
     */
    @Ensures("result != null")
    public List<T> getSitesAtLocation(final GenomeLoc loc) {
        if ( loc.size() != 1 ) throw new IllegalArgumentException("GenomeLoc must be a single one bp site but got " + loc);

        final List<T> l = new LinkedList<T>();

        while ( hasNext() ) {
            final T n = peek();
            final GenomeLoc nLoc = n.getLocation(parser);
            if ( nLoc.isBefore(loc) ) { // skip sites before loc
                next(); // read the value, and continue
            } else if ( nLoc.startsAt(loc) ) { // this site startsAt loc
                l.add(next());
            } else {
                break;
            }
        }

        return l;
    }

    /**
     * Get the list of all the upcoming T that are considered equivalent
     *
     * Equivalent means has the same contig, start, stop, ref and alt alleles
     *
     * @return a list of all of equivalent calls at a site, may be empty if there's nothing left to do
     */
    @Ensures("result != null")
    public List<T> getNextEquivalents() {
        final LinkedList<T> variants = new LinkedList<T>();

        final T first = peek();
        while ( hasNext() ) {
            final T n = peek();
            if ( equivalentVariants(first, n) )
                variants.add(next());
            else
                break;
        }

        return variants;
    }

    /**
     * Note that we require variants to have the same position, start, ref and alt allele
     *
     * @param vc1 first MongoVariantContext
     * @param vc2 second MongoVariantContext
     * @return true if vc1 and vc2 are equivalent, false otherwise
     */
    @Requires({"vc1 != null", "vc2 != null"})
    private boolean equivalentVariants(final T vc1, final T vc2) {
        return vc1.getChr().equals(vc2.getChr())
                && vc1.getStart() == vc2.getStart()
                && vc1.getRef().equals(vc2.getRef())
                && vc1.getAlt().equals(vc2.getAlt());
    }

    /**
     * Get all of the records in this sites as a list
     * @return a list containing all of the records in this iterator
     */
    public List<T> toList() {
        final List<T> l = new LinkedList<T>();

        while ( hasNext() ) {
            l.add(next());
        }

        return l;
    }

    @Override
    public boolean hasNext() {
        while ( super.hasNext() ) {
            final T n = super.peek();
            try {
                n.validate(parser);
                return true;
            } catch (MongoVariantContextException e) {
                handleError(n, e);
                super.next();
            }
        }

        return false;
    }

    private void handleError(final T record, final MongoVariantContextException e) {
        errorHandler.handleFailedRecord(record, e);
    }

    @Override
    public T next() {
        final T n = super.next();
        final GenomeLoc nLoc = n.getLocation(parser);
        if ( lastLoc != null && nLoc.isBefore(lastLoc) )
            throw new IllegalStateException("Records appearing out of order.  Current location is " + nLoc + " but last location was " + lastLoc);
        lastLoc = nLoc;
        return n;
    }

    @Override
    public Iterator<T> iterator() {
        return this;
    }
}

