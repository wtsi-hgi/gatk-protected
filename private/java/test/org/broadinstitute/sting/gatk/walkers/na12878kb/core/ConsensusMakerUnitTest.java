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

package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class ConsensusMakerUnitTest extends BaseTest {
    final static MongoVariantContext template = new MongoVariantContext(Arrays.asList("x"), "20", 1, 1, "A", "C", TruthStatus.TRUE_POSITIVE, new MongoGenotype(0, 0), new Date(), false);

    ConsensusMaker maker;

    MongoVariantContext reviewedX, reviewedX2, reviewedX3, reviewedA, notReviewedY, notReviewedZ, notReviewedZ2;
    MongoVariantContext reviewedTPX, reviewedFPX, reviewedFPY, reviewedTPZ, reviewedSUSPECT, notReviewedTP, notReviewedFP;
    static long date = 10;

    @BeforeMethod
    public void setUp() throws Exception {
        maker = new ConsensusMaker();

        reviewedX = copyTemplate("x", true);
        reviewedX2 = copyTemplate("x", true);
        reviewedX3 = copyTemplate("x", true);
        reviewedA = copyTemplate("A", true);
        notReviewedY = copyTemplate("y", false);
        notReviewedZ = copyTemplate("z", false);
        notReviewedZ2 = copyTemplate("z", false);

        // ***** IMPORTANT *****
        // the ordering of these lines is important because I am trying to replicate an issue we encountered [EB]
        reviewedTPZ = copyTemplate("z", true, TruthStatus.TRUE_POSITIVE);
        reviewedFPY = copyTemplate("Y", true, TruthStatus.FALSE_POSITIVE);
        reviewedTPX = copyTemplate("x", true, TruthStatus.TRUE_POSITIVE);
        reviewedFPX = copyTemplate("x", true, TruthStatus.FALSE_POSITIVE);

        reviewedSUSPECT = copyTemplate("z", true, TruthStatus.SUSPECT);
        notReviewedTP = copyTemplate("z", false, TruthStatus.TRUE_POSITIVE);
        notReviewedFP = copyTemplate("z", false, TruthStatus.FALSE_POSITIVE);
    }

    private static MongoVariantContext copyTemplate(final String name, final boolean reviewed, final TruthStatus status) {
        try {
            final MongoVariantContext mvc = template.clone();
            mvc.setSupportingCallSets(Arrays.asList(name));
            mvc.setReviewed(reviewed);
            mvc.setDate(new Date(date++));
            mvc.setTruth(status);
            mvc.setGt(new MongoGenotype(0, 1));
            return mvc;
        } catch ( CloneNotSupportedException e ) {
            throw new RuntimeException(e);
        }
    }

    private static MongoVariantContext copyTemplate(final String name, final boolean reviewed) {
        return copyTemplate(name, reviewed, TruthStatus.TRUE_POSITIVE);
    }

    @DataProvider(name = "MakerData")
    public Object[][] makeMakerData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final MongoGenotype het = new MongoGenotype(0, 1);
        final MongoGenotype noCall = new MongoGenotype(-1, -1);

        tests.add(new Object[]{TruthStatus.TRUE_POSITIVE, het, true, Arrays.asList(reviewedTPX)});
        tests.add(new Object[]{TruthStatus.TRUE_POSITIVE, het, true, Arrays.asList(reviewedTPX, notReviewedTP)});
        tests.add(new Object[]{TruthStatus.TRUE_POSITIVE, het, true, Arrays.asList(reviewedTPX, notReviewedFP)});
        tests.add(new Object[]{TruthStatus.SUSPECT, noCall, true, Arrays.asList(reviewedTPX, reviewedSUSPECT)});
        tests.add(new Object[]{TruthStatus.FALSE_POSITIVE, noCall, true, Arrays.asList(reviewedFPX, notReviewedFP)});
        tests.add(new Object[]{TruthStatus.FALSE_POSITIVE, noCall, true, Arrays.asList(reviewedFPX, notReviewedTP)});

        // special case -- takes FP because FP was created after TP and they are the same call set
        tests.add(new Object[]{TruthStatus.FALSE_POSITIVE, noCall, true, Arrays.asList(reviewedTPX, reviewedFPX)});
        tests.add(new Object[]{TruthStatus.DISCORDANT, noCall, true, Arrays.asList(reviewedTPX, reviewedFPY)});
        tests.add(new Object[]{TruthStatus.DISCORDANT, noCall, true, Arrays.asList(reviewedTPZ, reviewedFPY, reviewedTPX, reviewedFPX)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MakerData")
    public void testMake(final TruthStatus expectedTruth, final MongoGenotype expectedGT, final boolean expectedReviewed, final List<MongoVariantContext> supporting) throws Exception {
        final MongoVariantContext consensus = maker.makeConsensus(supporting);
        Assert.assertEquals(consensus.getType(), expectedTruth);
        Assert.assertEquals(consensus.isReviewed(), expectedReviewed, "isReviewed failed");
        Assert.assertEquals(consensus.getGt(), expectedGT);
    }

    @Test
    public void testConsensusGT() throws Exception {

    }

    @Test
    public void testCallsForConsensus() throws Exception {
        Assert.assertEquals(maker.selectCallsForConsensus(Arrays.asList(reviewedX)), Arrays.asList(reviewedX));
        Assert.assertEquals(maker.selectCallsForConsensus(Arrays.asList(notReviewedY)), Arrays.asList(notReviewedY));
        Assert.assertEquals(maker.selectCallsForConsensus(Arrays.asList(reviewedX, notReviewedY)), Arrays.asList(reviewedX));
        Assert.assertEquals(maker.selectCallsForConsensus(Arrays.asList(reviewedX, notReviewedY, reviewedA)), Arrays.asList(reviewedX, reviewedA));
        Assert.assertEquals(maker.selectCallsForConsensus(Arrays.asList(reviewedX, reviewedX2)), Arrays.asList(reviewedX2), "Should take the most recent calls for any call set");
        Assert.assertEquals(maker.selectCallsForConsensus(Arrays.asList(reviewedX, reviewedX2, reviewedX3)), Arrays.asList(reviewedX3), "Should take the most recent calls for any call set");
        Assert.assertEquals(maker.selectCallsForConsensus(Arrays.asList(reviewedA, reviewedX, reviewedX2)), Arrays.asList(reviewedA, reviewedX2), "Should take the most recent calls for any call set");
    }

    @Test
    public void testIsReviewed() throws Exception {
        Assert.assertEquals(maker.isReviewed(Arrays.asList(reviewedX)), true);
        Assert.assertEquals(maker.isReviewed(Arrays.asList(notReviewedY)), false);
        Assert.assertEquals(maker.isReviewed(Arrays.asList(notReviewedY, notReviewedZ)), false);
        Assert.assertEquals(maker.isReviewed(Arrays.asList(notReviewedY, reviewedX)), true);
        Assert.assertEquals(maker.isReviewed(Arrays.asList(notReviewedY, notReviewedZ, reviewedX)), true);
        Assert.assertEquals(maker.isReviewed(Arrays.asList(notReviewedY, reviewedX, notReviewedZ)), true);
        Assert.assertEquals(maker.isReviewed(Arrays.asList(notReviewedY, reviewedX, reviewedX2)), true);
    }

    @DataProvider(name = "TestConsensusGT")
    public Object[][] makeTestConsensusGT() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final List<Allele> alleles = Arrays.asList(Allele.create("A", true), Allele.create("C"));

        final MongoVariantContext tpHet = new MongoVariantContext(Arrays.asList("x"), "20", 1, 1, "A", "C",
                TruthStatus.TRUE_POSITIVE, new MongoGenotype(0, 1), new Date(), false);

        final MongoVariantContext tpNoCall = new MongoVariantContext(Arrays.asList("x"), "20", 1, 1, "A", "C",
                TruthStatus.TRUE_POSITIVE, new MongoGenotype(-1, -1), new Date(), false);

        final MongoVariantContext tpHet2 = new MongoVariantContext(Arrays.asList("y"), "20", 1, 1, "A", "C",
                TruthStatus.TRUE_POSITIVE, new MongoGenotype(0, 1), new Date(), false);

        final MongoVariantContext tpHetReviewed = new MongoVariantContext(Arrays.asList("z"), "20", 1, 1, "A", "C",
                TruthStatus.TRUE_POSITIVE, new MongoGenotype(0, 1), new Date(), true);

        final MongoVariantContext tpHomVar = new MongoVariantContext(Arrays.asList("x"), "20", 1, 1, "A", "C",
                TruthStatus.TRUE_POSITIVE, new MongoGenotype(1, 1), new Date(), false);

        final MongoVariantContext fpHet = new MongoVariantContext(Arrays.asList("x"), "20", 1, 1, "A", "C",
                TruthStatus.FALSE_POSITIVE, new MongoGenotype(0, 1), new Date(), false);

        final Genotype hetGT = tpHet.getGt().toGenotype(alleles);

        for ( final List<MongoVariantContext> l : Utils.makePermutations(Arrays.asList(tpHet, tpHet2, tpHomVar, fpHet), 2, false) ) {
            // false positive -> gt is no call
            tests.add(new Object[]{TruthStatus.FALSE_POSITIVE, l, MongoGenotype.NO_CALL});
        }

        for ( final MongoVariantContext o : Arrays.asList(tpHet, tpHet2, tpHomVar, tpHetReviewed) ) {
            final List<MongoVariantContext> l = Arrays.asList(o, tpNoCall);
            tests.add(new Object[]{TruthStatus.TRUE_POSITIVE, l, o.getGt().toGenotype(alleles)});
            final List<MongoVariantContext> lrev = Arrays.asList(tpNoCall, o);
            tests.add(new Object[]{TruthStatus.TRUE_POSITIVE, lrev, o.getGt().toGenotype(alleles)});
        }

        // hets are combined correctly
        tests.add(new Object[]{TruthStatus.TRUE_POSITIVE, Arrays.asList(tpHet, tpHet2), hetGT});

        // het + hom-var -> discordant
        final Genotype discordant = MongoGenotype.createDiscordant(hetGT);
        tests.add(new Object[]{TruthStatus.TRUE_POSITIVE, Arrays.asList(tpHet, tpHomVar), discordant});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TestConsensusGT")
    final void testConsensusGT(final TruthStatus truthStatus,
                               final Collection<MongoVariantContext> individualCalls,
                               final Genotype expectedGT) {
        final List<Allele> alleles = Arrays.asList(Allele.create("A", true), Allele.create("C"));
        final Genotype actualGT = maker.consensusGT(truthStatus, PolymorphicStatus.POLYMORPHIC, alleles, individualCalls);

        if ( expectedGT.getGQ() == MongoGenotype.DISCORDANT_GQ )
            Assert.assertEquals(actualGT.getGQ(), MongoGenotype.DISCORDANT_GQ, "Expected GT was discordant but didn't get a discordant result");
        else
            Assert.assertEquals(
                    new MongoGenotype(alleles, actualGT),
                    new MongoGenotype(alleles, expectedGT),
                    "Failed to create expected consensus GT");
    }
}
