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

    @Argument(fullName = "useTest", shortName = "useTest", doc = "If true, we will use the test knowledge base, suitable for development and unit testing", required=false)
    protected boolean useTest = false;

    @Argument(shortName = "reset", required=false)
    protected boolean resetDB = false;
}
