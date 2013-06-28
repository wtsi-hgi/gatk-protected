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

package org.broadinstitute.sting.utils.pairhmm;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class BandedLoglessPairHMMUnitTest extends BaseTest {
    private final static boolean DEBUG = false;
    private final File FILENAME = new File(privateTestDir + "likelihoods_2x101_100kb.txt.gz");
//    private final File FILENAME = new File("likelihoods_2x250_1mb.txt.gz");
    private final List<PairHMMTestData> results = new LinkedList<>();

    @BeforeClass
    public void setUp() throws Exception {
        for ( final List<PairHMMTestData> data : PairHMMTestData.readLikelihoods(FILENAME).values() )
            results.addAll(data);
    }

    @DataProvider(name = "InfoData")
    public Object[][] makeInfoData() {
        List<Object[]> tests = new ArrayList<>();

        final List<Integer> infosToTake = Arrays.asList();
        final int maxInfo = 10000000;
        int counter = 0;
        for ( final PairHMMTestData info : results ) {
            if ( counter < maxInfo && ( infosToTake.isEmpty() || infosToTake.contains(counter)))
                tests.add(new Object[]{counter, info});
            counter++;
        }

        return tests.toArray(new Object[][]{});
    }

    // this value (-20) depends on the band size and mle distance of our parameters
    private final static double SMALLEST_LIKELIHOOD_THAT_WILL_BE_EQUAL_TO_FULL_LOGLESS = -20;
    @Test(dataProvider = "InfoData", enabled = !DEBUG)
    public void testBandedLoglessHMM(final int i, final PairHMMTestData info) {
        //logger.warn("Testing " + info.ref + " against " + info.read);
        final BandedLoglessPairHMM hmm = new BandedLoglessPairHMM(5, 1e-20);
        final double myL = info.runHMM(hmm);
        logger.warn("hmm " + hmm);
        if ( info.log10l >= SMALLEST_LIKELIHOOD_THAT_WILL_BE_EQUAL_TO_FULL_LOGLESS )
            Assert.assertEquals(myL, info.log10l, 1e-3);
    }



    @DataProvider(name = "ConstructedData")
    public Object[][] makeConstructedData() {
        List<Object[]> tests = new ArrayList<>();

        final int bandSize = 5;
        final String allCommonBases = "ACGTAACCGGTTAAACCCGGGTTTAAAACCCCGGGGGTTTT";

        for ( int commonSize = 1; commonSize < allCommonBases.length(); commonSize++ ) {
            final String commonBases = allCommonBases.substring(0, commonSize);
            tests.add(new Object[]{"Read starts in middle of haplotype, goes to end",
                    bandSize,
                    new PairHMMTestData(Utils.dupString("A", 5*bandSize) + commonBases, commonBases, (byte)30)});
            tests.add(new Object[]{"Read starts in start of haplotype, doesn't have tail",
                    bandSize,
                    new PairHMMTestData(commonBases + Utils.dupString("A", 5*bandSize), commonBases, (byte)30)});
            tests.add(new Object[]{"Read starts in middle of haplotype, has both header and tail",
                    bandSize,
                    new PairHMMTestData(Utils.dupString("C", 5*bandSize) + commonBases + Utils.dupString("A", 5*bandSize), commonBases, (byte)30)});
        }

//        final int commonSize = allCommonBases.length();
//        final String commonBases = allCommonBases.substring(0, commonSize);
//        tests.add(new Object[]{"Read starts in middle of haplotype, goes to end",
//                    bandSize,
//                    new PairHMMTestData(Utils.dupString("A", 5*bandSize) + commonBases, commonBases, (byte)30)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ConstructedData", enabled = !DEBUG)
    public void testConstructedExample(final String name, final int bandSize, final PairHMMTestData info) {
        //logger.warn("Testing " + info.ref + " against " + info.read);
        final BandedLoglessPairHMM hmm = new BandedLoglessPairHMM(bandSize);
        final double myL = info.runHMM(hmm);
        final double standardL = info.runHMM(new LoglessPairHMM());
        logger.warn("hmm " + hmm);
        Assert.assertEquals(myL, standardL, 1e-3);
    }

    final static double DONT_KEEP_VALUE = 1e60;
    final static double MLE_VALUE = 1e100;

    @DataProvider(name = "UpdateBandData")
    public Object[][] makeUpdateBandData() {
        List<Object[]> tests = new ArrayList<>();

        final int readLen = 80;
        final int hapLen = 100;
        final byte[] hapBases = Utils.dupBytes((byte)'A', hapLen);
        final byte[] readBases = Utils.dupBytes((byte)'A', readLen);
        final byte[] quals = Utils.dupBytes((byte)30, readLen);

        final int bandSize = 5;
        boolean alreadyGeneratedMultiBandTest = false;
//        for ( final int readI : Arrays.asList(50) ) {
        for ( int readI = 1; readI < readLen; readI++ ) {
            final BandedLoglessPairHMM hmm = new BandedLoglessPairHMM(bandSize);
            hmm.initialize(readLen, hapLen);
            // necessary to ensure all of the internal state is set up correctly
            hmm.computeReadLikelihoodGivenHaplotypeLog10(hapBases, readBases, quals, quals, quals, quals, true);

            final Bands fullBand = new Bands();
            fullBand.addBand(1, hapLen+1);

            // test single band
            for ( int hapPos = 1; hapPos < hapLen; hapPos++ ) {
//            for ( final int hapPos : Arrays.asList(11) ) {
                final double[] bandData = new double[hapLen+1];
                Arrays.fill(bandData, DONT_KEEP_VALUE);
                bandData[hapPos] = MLE_VALUE;

                if ( readI < hmm.firstRowsBandSize() ) {
                    tests.add(new Object[]{hmm, bandData, readI, fullBand, fullBand});
                } else {
                    final Bands expected = new Bands();
                    expected.addPaddedBand(hapPos, hapPos, bandSize, hapLen+1);
                    tests.add(new Object[]{hmm, bandData, readI, fullBand, expected});
                }
            }

            // test multiple bands
            if ( readI >= hmm.firstRowsBandSize() && ! alreadyGeneratedMultiBandTest ) {
//                for ( int i : Arrays.asList(50) ) {
//                    for ( int j : Arrays.asList(60)) {
                for ( int i = 1; i < hapLen; i++ ) {
                    for ( int j = i + 1; j < hapLen; j++ ) {
                        final double[] bandData = new double[hapLen+1];
                        Arrays.fill(bandData, DONT_KEEP_VALUE);
                        bandData[i] = MLE_VALUE;
                        bandData[j] = MLE_VALUE;

                        final Bands expected = new Bands();
                        if ( i + bandSize >= j - bandSize ) {
                            expected.addPaddedBand(i, j, bandSize, hapLen+1);
                        } else {
                            expected.addPaddedBand(i, i, bandSize, hapLen+1);
                            expected.addPaddedBand(j, j, bandSize, hapLen+1);
                        }
                        tests.add(new Object[]{hmm, bandData, readI, fullBand, expected});
                    }
                }

                // only do this expense test once, since it doesn't depend on the hapPos
                alreadyGeneratedMultiBandTest = true;
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "UpdateBandData", enabled = true)
    public void testUpdateBand(final BandedLoglessPairHMM hmm,
                               final double[] bandData,
                               final int readPos,
                               final Bands currentBands,
                               final Bands expectedBands) {
        // set the band data for use int he updateBands calculation
        hmm.curRow.clear();
        for ( int j = 0; j < bandData.length; j++) {
            hmm.curRow.match[j] = bandData[j];
        }

        hmm.updateBands(currentBands, readPos, MLE_VALUE);
        Assert.assertEquals(hmm.currBands, expectedBands, "Expected bands " + expectedBands + " but found " + hmm.currBands);
    }
}
