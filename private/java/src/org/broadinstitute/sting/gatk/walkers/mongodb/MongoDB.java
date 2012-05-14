package org.broadinstitute.sting.gatk.walkers.mongodb;

import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.Mongo;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.net.UnknownHostException;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 4/26/12
 * Time: 3:01 PM
 * Handles Mongo DB connections
 */
final public class MongoDB {
    private final static String DB_HOST = "couchdb.broadinstitute.org";
    private final static Integer DB_PORT = 43054;
    private final static String DB_NAME = "nexus";
    private final static String SITES_COLLECTION = "sites";
    private final static String SAMPLES_COLLECTION = "samples";

    protected Mongo mongo;
    protected DBCollection sites;
    protected DBCollection samples;

    final private static MongoDB INSTANCE = new MongoDB();

    protected static DBCollection getSitesCollection() {
        return INSTANCE.sites;
    }

    protected static DBCollection getSamplesCollection() {
        return INSTANCE.samples;
    }

    protected static void close() {
        INSTANCE.mongo.close();
    }

    private MongoDB() {
        try {
            mongo = new Mongo(DB_HOST, DB_PORT);
            DB mongoDb = mongo.getDB(DB_NAME);
            sites = mongoDb.getCollection(SITES_COLLECTION);
            samples = mongoDb.getCollection(SAMPLES_COLLECTION);
        } catch (UnknownHostException e) {
            throw new StingException(e.getMessage(), e);
        }
    }
}
