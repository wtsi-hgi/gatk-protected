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

package org.broadinstitute.sting.gatk.walkers.na12878kb.assess;

import net.sf.samtools.SAMFileReader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.TruthStatus;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AssessorUnitTest extends BaseTest {
    private final static boolean DEBUG = false;

    private VariantContext makeVC(final int start, final String ... alleles) {
        return GATKVariantContextUtils.makeFromAlleles("vcf", "20", start, Arrays.asList(alleles));
    }

    private MongoVariantContext makeMVC(final VariantContext vc, final TruthStatus status) {
        final Genotype het = GenotypeBuilder.create("NA12878", vc.getAlleles());
        final MongoVariantContext mvc = MongoVariantContext.create("kb", vc, status, het);
        return mvc;
    }

    // ------------------------------------------------------------
    // Tests for assessing a site
    // ------------------------------------------------------------

    @DataProvider(name = "AssessSiteData")
    public Object[][] makeAssessSiteData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final VariantContext vcAC10 = makeVC(10, "A", "C");
        final VariantContext vcAT10 = makeVC(10, "A", "T");
        final VariantContext vcAG10 = makeVC(10, "A", "G");
        final VariantContext vcACT10 = makeVC(10, "A", "C", "T");
        final VariantContext vcNoVariation = makeVC(10, "A");

        final MongoVariantContext mvcAC10TP = makeMVC(vcAC10, TruthStatus.TRUE_POSITIVE);
        final MongoVariantContext mvcAT10TP = makeMVC(vcAT10, TruthStatus.TRUE_POSITIVE);
        final MongoVariantContext mvcAG10TP = makeMVC(vcAG10, TruthStatus.TRUE_POSITIVE);
        final MongoVariantContext mvcAC10FP = makeMVC(vcAC10, TruthStatus.FALSE_POSITIVE);
        final MongoVariantContext mvcAC10SUS = makeMVC(vcAC10, TruthStatus.SUSPECT);

        final List<MongoVariantContext> noKB = Collections.emptyList();

        final Assessment emptyAssessment = new Assessment(AssessmentType.DETAILED_ASSESSMENTS);
        final Assessment oneNovel = new Assessment(AssessmentType.DETAILED_ASSESSMENTS, AssessmentType.CALLED_NOT_IN_DB_AT_ALL);
        final Assessment oneTP = new Assessment(AssessmentType.DETAILED_ASSESSMENTS, AssessmentType.TRUE_POSITIVE);
        final Assessment oneFP = new Assessment(AssessmentType.DETAILED_ASSESSMENTS, AssessmentType.FALSE_POSITIVE_SITE_IS_FP);
        final Assessment oneFN = new Assessment(AssessmentType.DETAILED_ASSESSMENTS, AssessmentType.FALSE_NEGATIVE_NOT_CALLED_AT_ALL);
        final Assessment oneTN = new Assessment(AssessmentType.DETAILED_ASSESSMENTS, AssessmentType.CORRECTLY_UNCALLED);
        final Assessment oneNotRel = new Assessment(AssessmentType.DETAILED_ASSESSMENTS, AssessmentType.NOT_RELEVANT);

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{Arrays.asList(vcAC10), noKB, oneNovel, emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcAC10), Arrays.asList(mvcAC10TP), oneTP, emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcAC10), Arrays.asList(mvcAC10FP), oneFP, emptyAssessment});

        // false negatives
        tests.add(new Object[]{Arrays.asList(), Arrays.asList(mvcAC10TP), oneFN, emptyAssessment});
        tests.add(new Object[]{Arrays.asList(), Arrays.asList(mvcAC10FP), oneTN, emptyAssessment});

        // suspect sites
        tests.add(new Object[]{Arrays.asList(vcAC10), Arrays.asList(mvcAC10SUS), oneNotRel, emptyAssessment});
        tests.add(new Object[]{Arrays.asList(), Arrays.asList(mvcAC10SUS), oneNotRel, emptyAssessment});

        // make sure no variation is handled corrected
        tests.add(new Object[]{Arrays.asList(vcNoVariation), noKB, emptyAssessment, emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcNoVariation), Arrays.asList(mvcAC10TP), oneFN, emptyAssessment});

        // make sure multiple events are handled properly
        tests.add(new Object[]{Arrays.asList(vcAC10, vcAT10), Arrays.asList(mvcAC10TP), oneNovel.add(oneTP), emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcAC10, vcAT10), Arrays.asList(mvcAC10FP), oneNovel.add(oneFP), emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcAC10, vcAT10), Arrays.asList(mvcAC10FP, mvcAT10TP), oneTP.add(oneFP), emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcAC10, vcAT10), Arrays.asList(mvcAC10TP, mvcAT10TP), oneTP.add(oneTP), emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcAC10, vcAT10, vcAG10), Arrays.asList(mvcAC10TP, mvcAT10TP), oneTP.add(oneTP).add(oneNovel), emptyAssessment});

        tests.add(new Object[]{Arrays.asList(vcAT10), Arrays.asList(mvcAC10TP), oneNovel.add(oneFN), emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcAT10), Arrays.asList(mvcAC10TP, mvcAT10TP), oneFN.add(oneTP), emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcAT10), Arrays.asList(mvcAC10FP, mvcAT10TP), oneTN.add(oneTP), emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcAT10), Arrays.asList(mvcAC10FP, mvcAG10TP), oneTN.add(oneNovel).add(oneFN), emptyAssessment});

        // multi-allelic tests for SNPs
        tests.add(new Object[]{Arrays.asList(vcACT10), Arrays.asList(mvcAC10TP), oneNovel.add(oneTP), emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcACT10), Arrays.asList(mvcAC10TP, mvcAT10TP), oneTP.add(oneTP), emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcACT10), Arrays.asList(mvcAC10FP, mvcAT10TP), oneTP.add(oneFP), emptyAssessment});
        tests.add(new Object[]{Arrays.asList(vcACT10), Arrays.asList(mvcAC10SUS, mvcAT10TP), oneTP.add(oneNotRel), emptyAssessment});

        final VariantContext vcAC_C10 = makeVC(10, "AC", "A");
        final VariantContext vcA_AC10 = makeVC(10, "A", "AC");
        final VariantContext vcA_ACC10 = makeVC(10, "A", "ACC");
        final VariantContext vcAC_A_ACC = makeVC(10, "AC", "A", "ACC");
        final VariantContext vcAC_A_ACC_ACCC = makeVC(10, "AC", "A", "ACC", "ACCC");
        final MongoVariantContext mvcAC_C10_TP = makeMVC(vcAC_C10, TruthStatus.TRUE_POSITIVE);
        final MongoVariantContext mvcA_AC10_FP = makeMVC(makeVC(10, "A", "AC"), TruthStatus.FALSE_POSITIVE);
        final MongoVariantContext mvcA_ACC10_SUS = makeMVC(makeVC(10, "A", "ACC"), TruthStatus.SUSPECT);

        tests.add(new Object[]{Arrays.asList(vcAC_C10), Arrays.asList(), emptyAssessment, oneNovel});
        tests.add(new Object[]{Arrays.asList(vcAC_C10), Arrays.asList(mvcAC_C10_TP), emptyAssessment, oneTP});
        tests.add(new Object[]{Arrays.asList(vcA_AC10), Arrays.asList(mvcA_AC10_FP), emptyAssessment, oneFP});
        tests.add(new Object[]{Arrays.asList(vcA_ACC10), Arrays.asList(mvcA_ACC10_SUS), emptyAssessment, oneNotRel});

        tests.add(new Object[]{Arrays.asList(vcAC_C10), Arrays.asList(mvcA_AC10_FP), emptyAssessment, oneNovel.add(oneTN)});
        tests.add(new Object[]{Arrays.asList(vcAC_C10), Arrays.asList(mvcAC_C10_TP, mvcA_AC10_FP), emptyAssessment, oneTP.add(oneTN)});

        tests.add(new Object[]{Arrays.asList(vcAC_A_ACC), Arrays.asList(mvcA_AC10_FP), emptyAssessment, oneFP.add(oneNovel)});
        tests.add(new Object[]{Arrays.asList(vcAC_A_ACC), Arrays.asList(mvcA_AC10_FP, mvcAC_C10_TP), emptyAssessment, oneFP.add(oneTP)});

        tests.add(new Object[]{Arrays.asList(vcAC_A_ACC_ACCC), Arrays.asList(mvcAC_C10_TP), emptyAssessment, oneTP.add(oneNovel).add(oneNovel)});
        tests.add(new Object[]{Arrays.asList(vcAC_A_ACC_ACCC), Arrays.asList(mvcAC_C10_TP, mvcA_AC10_FP), emptyAssessment, oneTP.add(oneFP).add(oneNovel)});
        tests.add(new Object[]{Arrays.asList(vcAC_A_ACC_ACCC), Arrays.asList(mvcAC_C10_TP, mvcA_AC10_FP, mvcA_ACC10_SUS), emptyAssessment, oneTP.add(oneFP).add(oneNotRel)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AssessSiteData")
    public void testAssessSite(final List<VariantContext> fromVCF, final List<MongoVariantContext> fromKB, Assessment expectedSNPs, Assessment expectedIndels) {
        final Assessor assessor = new Assessor("test");
        assessor.accessSite(fromVCF, fromKB, false);
        Assert.assertEquals(assessor.getSNPAssessment(), expectedSNPs);
        Assert.assertEquals(assessor.getIndelAssessment(), expectedIndels);
    }

    @Test
    public void testSubsetMultiAllelic() {
        // This ensures that subsetting down multi-allelic to NA12878 works properly

        final Allele threeCopies = Allele.create("GTTTTATTTTATTTTA", true);
        final Allele twoCopies = Allele.create("GTTTTATTTTA", true);
        final Allele zeroCopies = Allele.create("G", false);
        final Allele oneCopies = Allele.create("GTTTTA", false);

        final VariantContextBuilder b = new VariantContextBuilder().source("foo").chr("20").start(10).stop(25);
        b.alleles(Arrays.asList(threeCopies, zeroCopies, oneCopies));
        b.genotypes(
                new GenotypeBuilder("NA12878", Arrays.asList(threeCopies, oneCopies)).make(),
                new GenotypeBuilder("NA12891", Arrays.asList(threeCopies, zeroCopies)).make());

        final VariantContext truthVC = new VariantContextBuilder("truth", "20", 10, 20, Arrays.asList(twoCopies, zeroCopies)).make();
        final MongoVariantContext truthMVC = makeMVC(truthVC, TruthStatus.TRUE_POSITIVE);

        final Assessment oneTP = new Assessment(AssessmentType.DETAILED_ASSESSMENTS, AssessmentType.TRUE_POSITIVE);

        final Assessor assessor = new Assessor("test");
        assessor.accessSite(Arrays.asList(b.make()), Arrays.asList(truthMVC), false);
        Assert.assertEquals(assessor.getIndelAssessment(), oneTP);
    }

    // ------------------------------------------------------------
    // Tests for assessing a site
    // ------------------------------------------------------------

    @DataProvider(name = "FilteringSiteData")
    public Object[][] makeFilteringSiteData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final VariantContext vcAC10 = makeVC(10, "A", "C");
        final VariantContext vcACC10 = makeVC(10, "AC", "C");

        for ( final VariantContext vcRoot : Arrays.asList(vcAC10, vcACC10) ) {
            final VariantContextBuilder b = new VariantContextBuilder(vcRoot);
            final MongoVariantContext mvcTP = makeMVC(vcRoot, TruthStatus.TRUE_POSITIVE);
            final MongoVariantContext mvcFP = makeMVC(vcRoot, TruthStatus.FALSE_POSITIVE);

            tests.add(new Object[]{b.make(), null, AssessmentType.CALLED_NOT_IN_DB_AT_ALL});
            tests.add(new Object[]{b.make(), mvcTP, AssessmentType.TRUE_POSITIVE});
            tests.add(new Object[]{b.make(), mvcFP, AssessmentType.FALSE_POSITIVE_SITE_IS_FP});

            tests.add(new Object[]{new VariantContextBuilder(vcRoot).filter("Filtered").make(), mvcTP, AssessmentType.FALSE_NEGATIVE_CALLED_BUT_FILTERED});
            tests.add(new Object[]{new VariantContextBuilder(vcRoot).filter("Filtered").make(), mvcFP, AssessmentType.CORRECTLY_FILTERED});
            tests.add(new Object[]{new VariantContextBuilder(vcRoot).filter("Filtered").make(), null, AssessmentType.NOT_RELEVANT});
            tests.add(new Object[]{new VariantContextBuilder(vcRoot).attribute("FS", 1000.0).make(), mvcTP, AssessmentType.TRUE_POSITIVE});
            tests.add(new Object[]{new VariantContextBuilder(vcRoot).attribute("FS", 1000.0).make(), mvcFP, AssessmentType.REASONABLE_FILTERS_WOULD_FILTER_FP_SITE});
            tests.add(new Object[]{new VariantContextBuilder(vcRoot).attribute("QD", 0.0).make(), mvcTP, AssessmentType.TRUE_POSITIVE});
            tests.add(new Object[]{new VariantContextBuilder(vcRoot).attribute("QD", 0.0).make(), mvcFP, AssessmentType.REASONABLE_FILTERS_WOULD_FILTER_FP_SITE});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "FilteringSiteData")
    public void testFilteringSites(final VariantContext vc, final MongoVariantContext mvc, final AssessmentType expectedType) {
        final Assessor assessor = new Assessor("test");
        assessor.accessSite(Collections.singletonList(vc), mvc == null ? Collections.<MongoVariantContext>emptyList() : Collections.singletonList(mvc), false);
        final Assessment actual = vc.isSNP() ? assessor.getSNPAssessment() : assessor.getIndelAssessment();
        final Assessment expected = new Assessment(AssessmentType.DETAILED_ASSESSMENTS, expectedType);
        Assert.assertEquals(actual, expected);
    }

    // ------------------------------------------------------------
    // Tests for assessing a site
    // ------------------------------------------------------------

    // java -jar dist/GenomeAnalysisTK.jar -T DepthOfCoverage -I private/testdata/reduced.readNotFullySpanningDeletion.bam -R ~/Desktop/broadLocal/localData/human_g1k_v37.fasta -L 1:167,022,605-167,025,904 | grep "1:" | grep -v "-" | awk '{print $1, $2}' | awk -F ":" '{print $1, $2}' > private/testdata/reduced.readNotFullySpanningDeletion.doc
    private final static File DoC = new File(privateTestDir + "reduced.readNotFullySpanningDeletion.doc");
    private final static File BAM = new File(privateTestDir + "reduced.readNotFullySpanningDeletion.bam");

    @DataProvider(name = "DoCSites")
    public Object[][] makeDocSites() throws Exception {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final String line : new XReadLines(DoC).readLines() ) {
            final String[] parts = line.split(" ");
            final String chr = parts[0];
            final int pos = Integer.valueOf(parts[1]);
            final int doc = Integer.valueOf(parts[2]);
            tests.add(new Object[]{BAM, chr, pos, doc});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "DoCSites")
    public void testFilteringSites(final File bam, final String chr, final int pos, final int expectedDoC) {
        final SAMFileReader bamReader = Assessor.makeSAMFileReaderForDoCInBAM(bam);
        final Assessor assessor = new Assessor("test", AssessNA12878.TypesToInclude.BOTH, Collections.<String>emptySet(), BadSitesWriter.NOOP_WRITER, bamReader, 5);
        final int actualDoC = assessor.getDepthAtLocus(chr, pos);
        Assert.assertEquals(actualDoC, expectedDoC, "Depth of coverage at " + chr + ":" + pos + " had unexpected depth");
    }
}
