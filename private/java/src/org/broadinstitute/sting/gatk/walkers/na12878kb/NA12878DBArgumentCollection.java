package org.broadinstitute.sting.gatk.walkers.na12878kb;

import com.google.gson.Gson;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.io.*;

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

    @Argument(fullName = "dbSpecPath", shortName = "dbSpecPath", doc = "Path to JSON file specifying database connection parameters", required=false)
    protected String dbSpecPath = getDBSpecPath(useLocal);

    private MongoDBManager.Locator locator;

    public MongoDBManager.Locator getLocator(){
        return locator;
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

        private DBType(String extension) {
            this.extension = extension;
        }

        public String getExtension() {
            return extension;
        }
    }

    @Argument(fullName = "dbToUse", shortName = "dbToUse", doc = "Which database should we connect to?", required=false)
    protected DBType dbToUse = DBType.DEFAULT;

    @Argument(shortName = "reset", required=false)
    protected boolean resetDB = false;

    /**
     * Calls {@link org.broadinstitute.sting.gatk.walkers.na12878kb.NA12878DBArgumentCollection#NA12878DBArgumentCollection(false)}
     */
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
        this.dbSpecPath = getDBSpecPath(useLocal);
        initLocator(dbSpecPath);
    }

    /**
     * Uses the provided file to load database host/name/port/etc.
     * Must be appropriate JSON format
     * @param dbSpecPath
     */
    public NA12878DBArgumentCollection(String dbSpecPath){
        this.dbSpecPath = dbSpecPath;
        initLocator(dbSpecPath);

    }

    private void initLocator(String dbSpecPath){
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
        this.locator = new MongoDBManager.Locator(tmpLocator.host, tmpLocator.port, dbName, tmpLocator.sitesCollection,
                tmpLocator.callsetsCollection, tmpLocator.consensusCollection);
    }

}
