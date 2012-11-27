package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.commandline.Argument;

/**
 * Standard arguments for interacting with the NA12878 DB
 *
 * User: depristo
 * Date: 11/15/12
 */
public class NA12878DBArgumentCollection {
    @Argument(fullName = "useLocal", shortName = "useLocal", doc = "If true, the localhost MongoDB will be used; for testing only", required=false)
    protected boolean useLocal = false;

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
}
