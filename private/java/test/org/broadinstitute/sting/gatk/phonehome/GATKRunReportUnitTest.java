/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.phonehome;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.activeregionqc.CountReadsInActiveRegions;
import org.broadinstitute.sting.gatk.walkers.qc.CountLoci;
import org.broadinstitute.sting.gatk.walkers.qc.CountRODs;
import org.broadinstitute.sting.gatk.walkers.qc.CountReads;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.jets3t.service.S3Service;
import org.jets3t.service.S3ServiceException;
import org.jets3t.service.ServiceException;
import org.jets3t.service.model.S3Object;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GATKRunReportUnitTest extends BaseTest {
    private final static boolean DEBUG = false;
    // WARNING WARNING WARNING WARNING WARNING WARNING WARNING -- do not distribute this code
    // WARNING WARNING WARNING WARNING WARNING WARNING WARNING -- do not distribute this code
    // WARNING WARNING WARNING WARNING WARNING WARNING WARNING -- do not distribute this code
    private final static String downloaderAccessKey = "AKIAJ6L5JRU4FZWN3YKQ";
    private final static String downloaderSecretKey = "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";
    // WARNING WARNING WARNING WARNING WARNING WARNING WARNING -- do not distribute this code
    // WARNING WARNING WARNING WARNING WARNING WARNING WARNING -- do not distribute this code
    // WARNING WARNING WARNING WARNING WARNING WARNING WARNING -- do not distribute this code

    private static final long S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING = 30 * 1000;

    private Walker walker;
    private Exception exception;
    private GenomeAnalysisEngine engine;


    @BeforeClass
    public void setup() {
        walker = new CountReads();
        exception = new IllegalArgumentException("javaException");
        engine = new GenomeAnalysisEngine();
        engine.setArguments(new GATKArgumentCollection());
    }

    @Test(enabled = ! DEBUG)
    public void testAWSKeysAreValid() {
        // throws an exception if they aren't
        GATKRunReport.checkAWSAreValid();
    }

    @Test(enabled = ! DEBUG)
    public void testAccessKey() throws Exception {
        testAWSKey(GATKRunReport.getAWSUploadAccessKey(), GATKRunReport.AWS_ACCESS_KEY_MD5);
    }

    @Test(enabled = ! DEBUG)
    public void testSecretKey() throws Exception {
        testAWSKey(GATKRunReport.getAWSUploadSecretKey(), GATKRunReport.AWS_SECRET_KEY_MD5);
    }

    private void testAWSKey(final String accessKey, final String expectedMD5) throws Exception {
        Assert.assertNotNull(accessKey, "AccessKey should not be null");
        final String actualmd5 = Utils.calcMD5(accessKey);
        Assert.assertEquals(actualmd5, expectedMD5);
    }

    @DataProvider(name = "GATKReportCreationTest")
    public Object[][] makeGATKReportCreationTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Walker readWalker = new CountReads();
        final Walker lociWalker = new CountLoci();
        final Walker rodWalker = new CountRODs();
        final Walker artWalker = new CountReadsInActiveRegions();

        final Exception noException = null;
        final Exception javaException = new IllegalArgumentException("javaException");
        final Exception stingException = new ReviewedStingException("StingException");
        final Exception userException = new UserException("userException");

        final GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        engine.setArguments(new GATKArgumentCollection());

        for ( final Walker walker : Arrays.asList(readWalker, lociWalker, rodWalker, artWalker) ) {
            for ( final Exception exception : Arrays.asList(noException,  javaException, stingException, userException) ) {
                tests.add(new Object[]{walker, exception, engine});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "GATKReportCreationTest")
    public void testGATKReportCreationReadingAndWriting(final Walker walker, final Exception exception, final GenomeAnalysisEngine engine) throws Exception {
        final GATKRunReport report = new GATKRunReport(walker, exception, engine, GATKRunReport.PhoneHomeOption.STANDARD);
        final ByteArrayOutputStream captureStream = new ByteArrayOutputStream();
        final boolean succeeded = report.postReportToStream(captureStream);
        Assert.assertTrue(succeeded, "Failed to write report to stream");
        Assert.assertFalse(report.exceptionOccurredDuringPost(), "Post succeeded but report says it failed");
        Assert.assertNull(report.getErrorMessage(), "Post succeeded but there was an error message");
        Assert.assertNull(report.getErrorThrown(), "Post succeeded but there was an error message");
        final InputStream readStream = new ByteArrayInputStream(captureStream.toByteArray());

        GATKRunReport deserialized = null;
        try {
            deserialized = GATKRunReport.deserializeReport(readStream);
        } catch ( Exception e ) {
            final String reportString = new String(captureStream.toByteArray());
            Assert.fail("Failed to deserialize GATK report " + reportString + " with exception " + e);
        }

        if ( deserialized != null )
            Assert.assertEquals(report, deserialized);
    }

    @DataProvider(name = "GATKAWSReportMode")
    public Object[][] makeGATKAWSReportMode() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final GATKRunReport.AWSMode mode : GATKRunReport.AWSMode.values() ) {
            tests.add(new Object[]{mode});
        }

        return tests.toArray(new Object[][]{});
    }

    // Will fail with timeout if AWS time out isn't working
    // Will fail with exception if AWS doesn't protect itself from errors
    @Test(enabled = ! DEBUG, dataProvider = "GATKAWSReportMode", timeOut = S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING * 2)
    public void testAWS(final GATKRunReport.AWSMode awsMode) {
        logger.warn("Starting testAWS mode=" + awsMode);

        // Use a shorter timeout than usual when we're testing GATKRunReport.AWSMode.TIMEOUT
        final long thisTestS3Timeout = awsMode == GATKRunReport.AWSMode.TIMEOUT ? 30 * 1000 : S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING;
        final GATKRunReport report = new GATKRunReport(walker, exception, engine, GATKRunReport.PhoneHomeOption.STANDARD, thisTestS3Timeout);
        report.sendAWSToTestBucket();
        report.setAwsMode(awsMode);
        final S3Object s3Object = report.postReportToAWSS3();

        if ( awsMode == GATKRunReport.AWSMode.NORMAL ) {
            Assert.assertNotNull(s3Object, "Upload to AWS failed, s3Object was null. error was " + report.formatError());
            Assert.assertFalse(report.exceptionOccurredDuringPost(), "The upload should have succeeded but the report says it didn't.  Error was " + report.formatError());
            Assert.assertNull(report.getErrorMessage(), "Report succeeded but an error message was found");
            Assert.assertNull(report.getErrorThrown(), "Report succeeded but an thrown error was found");
            try {
                final GATKRunReport deserialized = GATKRunReport.deserializeReport(downloaderAccessKey, downloaderSecretKey, report.getS3ReportBucket(), s3Object);
                Assert.assertEquals(report, deserialized);
                deleteFromS3(report);
            } catch ( Exception e ) {
                Assert.fail("Failed to read, deserialize, or delete GATK report " + s3Object.getName() + " with exception " + e);
            }
        } else {
            Assert.assertNull(s3Object, "AWS upload should have failed for mode " + awsMode + " but got non-null s3 object back " + s3Object + " error was " + report.formatError());
            Assert.assertTrue(report.exceptionOccurredDuringPost(), "S3 object was null but the report says that the upload succeeded");
            Assert.assertNotNull(report.getErrorMessage(), "Report succeeded but an error message wasn't found");
            if ( awsMode == GATKRunReport.AWSMode.FAIL_WITH_EXCEPTION )
                Assert.assertNotNull(report.getErrorThrown());
        }
    }

    private void deleteFromS3(final GATKRunReport report) throws Exception {
        final S3Service s3Service = GATKRunReport.initializeAWSService(downloaderAccessKey, downloaderSecretKey);
        // Retrieve the whole data object we created previously
        s3Service.deleteObject(report.getS3ReportBucket(), report.getReportFileName());
    }

    @DataProvider(name = "PostReportByType")
    public Object[][] makePostReportByType() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final GATKRunReport.PhoneHomeOption et : GATKRunReport.PhoneHomeOption.values() ) {
            tests.add(new Object[]{et});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = ! DEBUG, dataProvider = "PostReportByType", timeOut = S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING * 2)
    public void testPostReportByType(final GATKRunReport.PhoneHomeOption type) {
        final GATKRunReport report = new GATKRunReport(walker, exception, engine, GATKRunReport.PhoneHomeOption.STANDARD, S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING);
        Assert.assertFalse(report.exceptionOccurredDuringPost(), "An exception occurred during posting the report");
        final boolean succeeded = report.postReport(type);

        if ( type == GATKRunReport.PhoneHomeOption.NO_ET )
            Assert.assertFalse(succeeded, "NO_ET option shouldn't write a report");
        else {
            Assert.assertTrue(succeeded, "Any non NO_ET option should succeed in writing a report");

            if ( type == GATKRunReport.PhoneHomeOption.STDOUT ) {
                // nothing to do
            } else if ( type == GATKRunReport.PhoneHomeOption.STANDARD && ! report.wentToAWS()) {
                final boolean wasDeleted = report.getLocalReportFullPath().delete();
                Assert.assertTrue(wasDeleted, "Couldn't delete a supposedly written local GATK report");
            } else {
                // must have gone to AWS
                try {
                    Assert.assertTrue(report.wentToAWS(), "The report should have gone to AWS but the report says it wasn't");
                    deleteFromS3(report);
                } catch ( Exception e ) {
                    Assert.fail("Failed delete GATK report " + report.getReportFileName() + " with exception " + e);
                }
            }
        }
    }

    public interface S3Op {
        public void apply() throws ServiceException;
    }

    // Will fail with timeout if AWS time out isn't working
    // Will fail with exception if AWS doesn't protect itself from errors
    @Test(timeOut = S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING * 2)
    public void testAWSPublicKeyHasAccessControls() throws Exception {
        final GATKRunReport report = new GATKRunReport(walker, exception, engine, GATKRunReport.PhoneHomeOption.STANDARD, S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING);
        report.sendAWSToTestBucket();
        final S3Object s3Object = report.postReportToAWSS3();
        Assert.assertNotNull(s3Object, "Upload to AWS failed, s3Object was null. error was " + report.formatError());

        // create a service with the public key, and make sure it cannot list or delete
        final S3Service s3Service = GATKRunReport.initializeAWSService(GATKRunReport.getAWSUploadAccessKey(), GATKRunReport.getAWSUploadSecretKey());
        assertOperationNotAllowed("listAllBuckets", new S3Op() {
            @Override
            public void apply() throws S3ServiceException {
                s3Service.listAllBuckets();
            }
        });
        assertOperationNotAllowed("listBucket", new S3Op() {
            @Override
            public void apply() throws S3ServiceException { s3Service.listObjects(report.getS3ReportBucket()); }
        });
        assertOperationNotAllowed("createBucket", new S3Op() {
            @Override
            public void apply() throws S3ServiceException { s3Service.createBucket("ShouldNotCreate"); }
        });
        assertOperationNotAllowed("deleteObject", new S3Op() {
            @Override
            public void apply() throws ServiceException { s3Service.deleteObject(report.getS3ReportBucket(), report.getReportFileName()); }
        });
    }

    private void assertOperationNotAllowed(final String name, final S3Op op) {
        try {
            op.apply();
            // only gets here if the operation was successful
            Assert.fail("Operation " + name + " ran successfully but we expected to it fail");
        } catch ( ServiceException e ) {
            Assert.assertEquals(e.getErrorCode(), "AccessDenied");
        }
    }
}
