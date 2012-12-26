package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import com.google.gson.Gson;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.io.*;
import java.util.Random;

/**
 * Standard arguments for interacting with the NA12878 DB
 *
 * User: depristo
 * Date: 11/15/12
 */
public class NA12878DBArgumentCollection {

    public static final String DEFAULT_SPEC_PATH = "resources/NA12878kb.json";
    public static final String LOCALHOST_SPEC_PATH = "resources/NA12878kb_local.json";

    @Argument(fullName = "useLocal", shortName = "useLocal", doc = "If true, the localhost MongoDB will be used; for testing only", required=false)
    protected boolean useLocal = false;

    /** Lazy loaded to work with GATK style argument value injection */
    String dbSpecPath = null;

    public MongoDBManager.Locator getLocator(){
        if ( dbSpecPath == null )
            dbSpecPath = getDBSpecPath(useLocal);
        InputStream is = getClass().getResourceAsStream(dbSpecPath);
        if(is == null){
            try {
                is = new FileInputStream(dbSpecPath);
            } catch (FileNotFoundException e) {
                throw new StingException("db spec path not found", e);
            }
        }

        Reader reader = new InputStreamReader(is);
        MongoDBManager.Locator tmpLocator = (new Gson()).fromJson(reader, MongoDBManager.Locator.class);
        String dbName = tmpLocator.name + dbToUse.getExtension();
        try {
            reader.close();
            return new MongoDBManager.Locator(tmpLocator.host, tmpLocator.port, dbName, tmpLocator.sitesCollection,
                    tmpLocator.callsetsCollection, tmpLocator.consensusCollection);
        } catch ( IOException e ) {
            throw new RuntimeException("Failed to close json reader for " + dbSpecPath, e);
        }
    }

    public enum DBType {
        /** The production database, contains the stable data for analysis */
        PRODUCTION("_production"),
        /** A persistent development database, for playing with the KB itself */
        DEV("_development"),
        /** For unit and integration tests, not persistent */
        TEST("_test"),
        /** default one */
        DEFAULT("_NA");

        private String extension;

        // Don't use the GATK fixed random seed -- we want to avoid namespace collisions between multiple GATK instances
        private static final Random rand = new Random();

        private DBType(String extension) {
            this.extension = extension;
        }

        public String getExtension() {
            // If using a test database, make an over-the-top effort to avoid naming collisions with
            // concurrently-running instances of the test suite
            return extension.equals("_test") ? String.format("%s_%d_%d", extension, System.currentTimeMillis(), rand.nextInt()) :
                                               extension;
        }
    }

    @Argument(fullName = "dbToUse", shortName = "dbToUse", doc = "Which database should we connect to?", required=false)
    public DBType dbToUse = DBType.DEFAULT;

    @Argument(shortName = "reset", required=false)
    public boolean resetDB = false;

    public NA12878DBArgumentCollection(){
        this(false);
    }

    private static String getDBSpecPath(boolean useLocal){
        return useLocal ? LOCALHOST_SPEC_PATH : DEFAULT_SPEC_PATH;
    }

    /**
     * Convenience constructor for choosing between local or remote
     * NA12878KB database
     * @param useLocal
     */
    public NA12878DBArgumentCollection(boolean useLocal){
        this.useLocal = useLocal;
    }

    /**
     * Uses the provided file to load database host/name/port/etc.
     * Must be appropriate JSON format
     *
     * This *forces* the dbToUse to be the production database
     *
     * @param dbSpecPath
     */
    public NA12878DBArgumentCollection(String dbSpecPath){
        this.dbSpecPath = dbSpecPath;
        this.dbToUse = DBType.PRODUCTION;
    }
}
