package org.broadinstitute.sting.gatk.walkers.na12878kb;

import junit.framework.Assert;
import org.broadinstitute.sting.BaseTest;
import org.testng.annotations.Test;

/**
 * User: jacob
 * Date: 2012-Nov-27
 */
public class NA12878DBArgumentCollectionUnitTest extends BaseTest {

    //We have two spec files saved. Test that we can read them, and
    //the fields are as expected
    @Test
    public void testCompareLocalRemoteLocators() throws Exception {
        NA12878DBArgumentCollection args = new NA12878DBArgumentCollection(true);
        MongoDBManager.Locator localLocator = args.getLocator();

        Assert.assertNotNull(localLocator);
        Assert.assertNotNull(localLocator.host);
        Assert.assertNotNull(localLocator.port);
        Assert.assertNotNull(localLocator.name);
        Assert.assertNotNull(localLocator.callsetsCollection);
        Assert.assertNotNull(localLocator.consensusCollection);
        Assert.assertNotNull(localLocator.sitesCollection);

        NA12878DBArgumentCollection args1 = new NA12878DBArgumentCollection(false);
        MongoDBManager.Locator remoteLocator = args1.getLocator();


        Assert.assertNotSame(localLocator, remoteLocator);
        Assert.assertNotSame(localLocator.host, remoteLocator.host);
        Assert.assertEquals(localLocator.name, remoteLocator.name);
        Assert.assertEquals(localLocator.port, remoteLocator.port);
        Assert.assertEquals(localLocator.callsetsCollection, remoteLocator.callsetsCollection);
        Assert.assertEquals(localLocator.consensusCollection, remoteLocator.consensusCollection);
        Assert.assertEquals(localLocator.callsetsCollection, remoteLocator.sitesCollection);
    }
}
