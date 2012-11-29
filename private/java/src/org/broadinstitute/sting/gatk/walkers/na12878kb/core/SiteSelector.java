package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import com.mongodb.BasicDBObject;
import com.mongodb.DBObject;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * Create queries on a MongoDB from DBCollections containing MongoVariantContexts
 *
 * Add intervals, types, sets, and isReview requirements, and then request
 * a DBObject suitable for a MongoDB.find() call to return just sites satisfying
 * those requirements.
 *
 * Follows the builder model, so add() functions return the SiteSelector itself
 *
 * User: depristo
 * Date: 11/4/12
 * Time: 8:16 AM
 */
public class SiteSelector {
    private static Logger logger = Logger.getLogger(SiteSelector.class);

    final GenomeLocParser parser;
    GenomeLocSortedSet intervals = null;
    final Set<CallSet> setsToInclude = new HashSet<CallSet>();
    final Set<TruthStatus> typesToInclude = new HashSet<TruthStatus>();
    boolean onlyReviewed = false;

    /**
     * Create a new SiteSelector initialized to select by default all records
     *
     * @param parser A GenomeLocParser to make GenomeLocs as necessary
     */
    public SiteSelector(final GenomeLocParser parser) {
        if ( parser == null )
            throw new IllegalArgumentException("GenomeLocParser cannot be null");
        this.parser = parser;
    }

    /**
     * Select only records within the intervals specified by locs
     *
     * Note that this function is cumulative, so subsequent calls *add* intervals
     * for selection.  So adding [1, 2] and [10-11] will result in a query
     * that returns records within either interval.
     *
     * @param locs the locations we wish to see records within
     * @return this SiteSelector
     */
    public SiteSelector addIntervals(final GenomeLocSortedSet locs) {
        if ( locs == null )
            throw new IllegalArgumentException("Locs cannot be null");

        for ( final GenomeLoc loc : locs )
            addInterval(loc);
        return this;
    }

    /**
     * Select only records within the interval specified by loc
     *
     * Note that this function is cumulative, so subsequent calls *add* intervals
     * for selection.  So adding [1, 2] and [10-11] will result in a query
     * that returns records within either interval.
     *
     * @param loc the location we wish to see records within
     * @return this SiteSelector
     */
    public SiteSelector addInterval(final GenomeLoc loc) {
        if ( loc == null )
            throw new IllegalArgumentException("Loc cannot be null");

        if ( intervals == null )
            intervals = new GenomeLocSortedSet(parser);

        intervals.addRegion(loc);
        return this;
    }

    /**
     * Convenience function equivalent to #addInterval
     * @param contig
     * @param start
     * @param stop
     * @return
     */
    public SiteSelector addInterval(final String contig, final int start, final int stop) {
        return addInterval(parser.createGenomeLoc(contig, start, stop));
    }

    /**
        * Only include calls within the callset set
        *
        * Note that this function is cumulative, so subsequent calls *add* sets
        * for selection.  So adding X and Y will result in a query
        * that returns records within either set.
        *
        * @param set the call set we wish to include
        * @return this SiteSelector
        */
    public SiteSelector addSetToInclude(final CallSet set) {
        setsToInclude.add(set);
        return this;
    }

    /**
     * Only include calls having TruthStatus type
     *
     * Note that this function is cumulative, so subsequent calls *add* sets
     * for selection.  So adding X and Y will result in a query
     * that returns records within either set.
     *
     * @param type the type of call we wish to include
     * @return this SiteSelector
     */
    public SiteSelector addTypeToInclude(final TruthStatus type) {
        typesToInclude.add(type);
        return this;
    }

    /**
     * Only include reviewed sites
     *
     * @return this SiteSelector
     */
    public SiteSelector onlyReviewed() {
        onlyReviewed = true;
        return this;
    }

    public GenomeLocSortedSet getIntervals() {
        return intervals;
    }

    public Set<CallSet> getSetsToInclude() {
        return setsToInclude;
    }

    public Set<TruthStatus> getTypesToInclude() {
        return typesToInclude;
    }

    /**
     * Return a newly allocated DBObject selecting only sites consistent with the current selections
     *
     * This DBObject can be handed to collection.find(this.toQuery()) to return only the sites in
     * collection meeting the currently active selections (given by addX functions).  Collection
     * must contain MongoDB encoded MongoVariantContext objects, as this function assumes
     * that we are selecting MongoVariantContext objects.
     *
     * @return a DBObject suitable for mongodb.collection.find()
     */
    public DBObject toQuery() {
        // a list of all of the conditions to add to the query
        final List<DBObject> conditions = new LinkedList<DBObject>();

        // add the interval restrictions
        if ( intervals != null ) {
            final List<DBObject> regionsToOr = new LinkedList<DBObject>();

            for ( final GenomeLoc interval : intervals.toList() ) {
                final DBObject StartRange = new BasicDBObject("$gte", interval.getStart()).append("$lte", interval.getStop());
                regionsToOr.add(new BasicDBObject("Chr", interval.getContig()).append("Start", StartRange));
            }

            conditions.add(new BasicDBObject("$or", regionsToOr));
        }

        // add the only reviewed restriction
        if ( onlyReviewed ) {
            conditions.add(new BasicDBObject("Reviewed", true));
        }

        if ( ! getSetsToInclude().isEmpty() )
            throw new UnsupportedOperationException("SetsToInclude not yet implemented");

        if ( ! getTypesToInclude().isEmpty() )
            throw new UnsupportedOperationException("TypesToInclude not yet implemented");

        final DBObject query = conditions.isEmpty() ? new BasicDBObject() : new BasicDBObject("$and", conditions);
        if ( logger.isDebugEnabled() ) logger.debug("Query " + query);
        return query;
    }
}
