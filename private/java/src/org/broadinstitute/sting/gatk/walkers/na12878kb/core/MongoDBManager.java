package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.Mongo;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.net.UnknownHostException;
import java.util.HashMap;
import java.util.Map;

/**
 * Manages connections to Mongo databases.
 * Connections are created and indexed with Locator instances.
 * @author: jacob
 * @date: Nov 27, 2012
 */
final public class MongoDBManager {
    private final static Logger logger = Logger.getLogger(MongoDBManager.class);


    private static Map<Locator, DBWrapper> connections = new HashMap<Locator, DBWrapper>(5);

    public static DBWrapper getDB(Locator locator){
        DBWrapper wrapper = connections.get(locator);
        if(wrapper == null){
            wrapper = new DBWrapper(locator);
            connections.put(locator, wrapper);
        }
        return wrapper;
    }

    protected static void closeAll(){
        for(DBWrapper wrapper: connections.values()){
            wrapper.close();
        }
    }

    public static class DBWrapper{
        protected Mongo mongo;
        protected DBCollection sites;
        protected DBCollection callsets;
        protected DBCollection consensus;

        protected Mongo getMongo(){
            return mongo;
        }

        protected DBCollection getSites() {
            return sites;
        }

        protected DBCollection getCallsets() {
            return callsets;
        }
        
        protected DBCollection getConsensus() {
            return consensus;
        }

        protected void close() {
            mongo.close();
        }

        private DBWrapper(Locator locator) {
            try {
                logger.info("Connecting to MongoDB host=" + locator.host + " port=" + locator.port + " name=" + locator.name);
                mongo = new Mongo(locator.host, locator.port);
                DB mongoDB = mongo.getDB(locator.name);
                sites = mongoDB.getCollection(locator.sitesCollection);
                callsets = mongoDB.getCollection(locator.callsetsCollection);
                consensus = mongoDB.getCollection(locator.consensusCollection);
            } catch (UnknownHostException e) {
                throw new StingException(e.getMessage(), e);
            }
        }
    }

    public static class Locator{
        public final String host;
        public final Integer port;
        public final String name;
        public final String sitesCollection;
        public final String callsetsCollection;
        public final String consensusCollection;

        public Locator(String host, int port, String name, String sitesCollection,
                       String callsetsCollection, String consensusCollection){
            this.host = host;
            this.port = port;
            this.name = name;
            this.sitesCollection = sitesCollection;
            this.callsetsCollection = callsetsCollection;
            this.consensusCollection = consensusCollection;
        }

        @Override
        public String toString() {
            return "Locator{" +
                    "host='" + host + '\'' +
                    ", port=" + port +
                    ", name='" + name + '\'' +
                    '}';
        }
    }
}
