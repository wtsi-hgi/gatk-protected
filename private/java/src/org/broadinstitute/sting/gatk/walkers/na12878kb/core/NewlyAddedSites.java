package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import com.mongodb.BasicDBObject;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.bson.types.ObjectId;

/**
 * Object that allows us to continually request newly added sites to the NA12878 Knowledge base
 *
 * Keeps track of the last object added to the sites collection, and provides get functions
 * that will return all of the objects added to the collection since the last object was added.
 * These get functions all update the pointer to the last object, so subsequent calls will
 * return new values.
 *
 * User: depristo
 * Date: 11/4/12
 * Time: 8:16 AM
 */
public class NewlyAddedSites {
    private static final DBObject SORT_BY_ID = new BasicDBObject("_id", -1);

    private static Logger logger = Logger.getLogger(NewlyAddedSites.class);

    /**
     * The knowledge base we are using to get sites / location
     */
    final NA12878KnowledgeBase kb;

    /**
     * The ObjectID of the last added record in the db.  May be null, indicating
     * that no objects were in the collection when this object was created.  Updated
     * each time a getNew* function is called, so that we keep returning
     * new objects from the collection upon subsequent queries to this object
     */
    ObjectId lastObjectID = null;

    /**
     * Create a NewlyAddedSites object that can be used to fetch newly added
     * records from the NA12878KnowledgeBase upon subsequent calls to get*
     * function here.
     *
     * @param kb the NA12878KnowledgeBase to get records from
     * @param lastObjectID the max ObjectID of the sites collection in kb, null if all
     *                     records should be treated as new
     */
    public NewlyAddedSites(final NA12878KnowledgeBase kb, final ObjectId lastObjectID) {
        if ( kb == null )
            throw new IllegalArgumentException("NA12878KnowledgeBase cannot be null");

        this.kb = kb;
        this.lastObjectID = lastObjectID;
    }

    /**
     * Create a NewlyAddedSites object that can be used to fetch newly added
     * records from the NA12878KnowledgeBase upon subsequent calls to get*
     * function here.
     *
     * The first newly added sites will be computed from the max ObjectID
     * of the sites in kb at the instant of this objects creation.
     *
     * @param kb the NA12878KnowledgeBase to get records from
     */
    public NewlyAddedSites(final NA12878KnowledgeBase kb) {
        this(kb, getLastObjectID(kb.getSitesCollection()));
    }

    /**
     * Get the newly added records in the knowledge base
     *
     * Note this only returns the actual newly added records, not all of the records
     * with the same chr/start as the newly added records
     *
     * @return a DBCursor iterating over all newly added records in the database
     */
    public DBCursor getNewlyAddedRecords() {
        final DBObject query = makeQuery();
        final DBCursor cursor = kb.getSitesCollection().find(query).sort(kb.sitesOrder());

        logger.debug("Query for newly added sites with lastObjectID " + lastObjectID + " returned " + cursor.size() + " records");

        // TODO -- might skip some sites if something is added between the above iteration and this query
        lastObjectID = getLastObjectID(kb.getSitesCollection());
        return cursor;
    }

    /**
     * Get a SiteSelector returning the sites of the newly added records
     *
     * Suitable for issuing a site selector query to return all of the records at sites
     * with newly added records.
     *
     * If the total number of queries in the SiteSelector would exceed maxItemizedQueries, then
     * getNewlyAddedLocations may return a SiteSelector returning *all* records in the db, as an
     * unrestricted query + processing may be faster than processing, for example, 1M range queries
     *
     * If there are no newly added records, return value is null
     *
     * @param parser the genome loc parser needed to create new GenomeLocs for the SiteSelector
     * @param maxItemizedQueries if the number of locations that we'd need to query to return all
     *                           newly added sites is greater than maxItemizedQueries, then all of the
     *                           records in the db will be selected instead.  This is essential as a
     *                           query of >100K locations will be slower than just updating the entire
     *                           DB in one go.  If -1 any number of location will be query.
     * @return null if no new locations have been added, or a SiteSelector that will query for all of the changed locations
     */
    public SiteSelector getNewlyAddedLocations(final GenomeLocParser parser, final int maxItemizedQueries) {
        if ( parser == null )
            throw new IllegalArgumentException("GenomeLocParser cannot be null");

        final DBCursor cursor = getNewlyAddedRecords();

        if ( cursor.size() == 0 ) {
            return null;
        } else if ( cursor.size() > maxItemizedQueries && maxItemizedQueries != -1 ) {
            return new SiteSelector(parser);
        } else {
            final SiteSelector selector = new SiteSelector(parser);
            for ( final MongoVariantContext mvc : new SiteIterator<MongoVariantContext>(parser, cursor) ) {
                selector.addInterval(mvc.getLocation(parser));
            }

            return selector;
        }
    }

    /**
     * Make a MongoDB query that will return all objects with ObjectID > lastObjectID
     *
     * Note that if lastObjectID == null (i.e., none was provided an there were no objects
     * in the sites collection) then this returns a query that returns all objects in the collection.
     *
     * @return a Query suitable to get all MongoDB objects with ObjectID > lastObjectID
     */
    private DBObject makeQuery() {
        if ( lastObjectID == null )
            return new BasicDBObject(); // find all objects
        else
            return new BasicDBObject("_id", new BasicDBObject("$gt", lastObjectID));
    }

    /**
     * Get the max object id for the DBCollection sites, or null if there are no records
     *
     * @param sites the collection to get the max object id of
     * @return the max ObjectID in collection, or null if there are no records in sites
     */
    private static ObjectId getLastObjectID(final DBCollection sites) {
        final DBCursor cursor = sites.find().sort(SORT_BY_ID).limit(1);
        if ( cursor.size() == 0 )
            return null;
        else
            return (ObjectId)cursor.next().get("_id");
    }

}
