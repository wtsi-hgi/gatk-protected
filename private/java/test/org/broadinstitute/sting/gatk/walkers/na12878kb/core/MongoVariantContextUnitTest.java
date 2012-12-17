package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.MongoVariantContextException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextTestProvider;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class MongoVariantContextUnitTest extends NA12878KBUnitTestBase {
    @DataProvider(name = "MVCBasicTest")
    public Object[][] makeMVCBasicTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final VariantContext vc : testVCs ) {
            final Genotype NO_CALL = MongoGenotype.NO_CALL;
            final Genotype HOMREF = MongoGenotype.create(vc, 0, 0);
            final Genotype HET = MongoGenotype.create(vc, 0, 1);
            final Genotype HOMVAR = MongoGenotype.create(vc, 1, 1);

            // date is passively tested before we get different dates each time
            for ( final String name : Arrays.asList("x", "y", "z") )
                tests.add(new Object[]{vc, new MongoVariantContext(name, vc, TruthStatus.TRUE_POSITIVE, new Date(), HET, true)});

            for ( final TruthStatus status : TruthStatus.values() )
                tests.add(new Object[]{vc, new MongoVariantContext("x", vc, status, new Date(), HET, true)});

            for ( final Genotype genotype : Arrays.asList(NO_CALL, HOMREF, HET, HOMVAR) )
                tests.add(new Object[]{vc, new MongoVariantContext("x", vc, TruthStatus.TRUE_POSITIVE, new Date(), genotype, true)});

            for ( final boolean reviewed : Arrays.asList(true, false) )
                tests.add(new Object[]{vc, new MongoVariantContext("x", vc, TruthStatus.TRUE_POSITIVE, new Date(), HET, reviewed)});
        }

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "TestVCProvider")
    public Object[][] testVCProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( final VariantContext vc : testVCs )
            tests.add(new Object[]{vc});
        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "TestVCProvider")
    public void testVCToMVC(final VariantContext vc) {
        final MongoVariantContext mvc = MongoVariantContext.create("x", vc, TruthStatus.UNKNOWN, MongoGenotype.NO_CALL);
        Assert.assertTrue(mvc.isSingleCallset());
        Assert.assertEquals(mvc.getCallSetName(), "x");
        Assert.assertEquals(mvc.getSupportingCallSets().size(), 1);
        Assert.assertEquals(mvc.getChr(), vc.getChr());
        Assert.assertEquals(mvc.getStart(), vc.getStart());
        Assert.assertEquals(mvc.getStop(), vc.getEnd());
        Assert.assertEquals(mvc.getRefAllele(), vc.getReference());
        Assert.assertEquals(mvc.getAltAllele(), vc.getAlternateAllele(0));
    }

    @Test(enabled = true, dataProvider = "TestVCProvider")
    public void testMVCMatching(final VariantContext vc) {
        final MongoVariantContext mvc = MongoVariantContext.create("x", vc, TruthStatus.UNKNOWN, MongoGenotype.NO_CALL);
        for ( final VariantContext vc2 : testVCs ) {
            final boolean expectedEquals = vc == vc2;
            Assert.assertEquals(mvc.matches(vc2), expectedEquals, "MVC match returned unexpected value for mvc=" + mvc + " vs vc2=" + vc2);
        }
    }

    @Test(enabled = true, dataProvider = "MVCBasicTest")
    public void testMVCBasic(final VariantContext originalVC, final MongoVariantContext mvc) {
        try {
            setupBeforeMethod();
            db.addCall(mvc);
            final MongoVariantContext fromDB = readOneMVCFromDB();

            Assert.assertEquals(fromDB, mvc, "Input MongoVariantContext not the same as the one read from DB");
            VariantContextTestProvider.assertEquals(fromDB.getVariantContext(), mvc.getVariantContext());
        }
        finally {
            teardownMethod();
        }
    }

    final MongoVariantContext readOneMVCFromDB() {
        final SiteIterator<MongoVariantContext> it = db.getCalls();
        Assert.assertTrue(it.hasNext(), "Expected at least 1 call in the db");
        final MongoVariantContext mvc = it.next();
        Assert.assertFalse(it.hasNext(), "Only expected 1 call in the db but saw > 1");
        it.close();
        return mvc;
    }

    final static MongoVariantContext good = new MongoVariantContext(Arrays.asList("x"), "20", 1, 1, "A", "C", TruthStatus.TRUE_POSITIVE, new MongoGenotype(0, 0), new Date(), false);
    private static MongoVariantContext makeBad(final List<MongoVariantContext> bads) throws CloneNotSupportedException {
        final MongoVariantContext mvc = good.clone();
        bads.add(mvc);
        return mvc;
    }

    public static List<MongoVariantContext> makeBadMVCs() {
        try {
            final List<MongoVariantContext> bads = new LinkedList<MongoVariantContext>();
            makeBad(bads).setSupportingCallSets(new ArrayList<String>());
            final ArrayList<String> l = new ArrayList<String>();
            makeBad(bads).setSupportingCallSets(l);
            final ArrayList<String> l2 = new ArrayList<String>();
            l2.add(null);
            makeBad(bads).setSupportingCallSets(l2);
            makeBad(bads).setChr("chr20");
            makeBad(bads).setChr("-1");
            makeBad(bads).setChr(null);
            makeBad(bads).setStart(-1);
            makeBad(bads).setStop(-1);
            makeBad(bads).setRef(null);
            makeBad(bads).setRef("");
            makeBad(bads).setRef("X");
            makeBad(bads).setRef("a");
            makeBad(bads).setAlt(null);
            makeBad(bads).setAlt("");
            makeBad(bads).setAlt("X");
            makeBad(bads).setAlt("a");
            makeBad(bads).setGt(new MongoGenotype(-1, -2));
            makeBad(bads).setGt(new MongoGenotype(-2, -1));
            makeBad(bads).setGt(new MongoGenotype(0, -1));
            makeBad(bads).setGt(new MongoGenotype(-1, 0));
            makeBad(bads).setGt(new MongoGenotype(2, 0));
            makeBad(bads).setGt(new MongoGenotype(0, 2));
            makeBad(bads).setGt(new MongoGenotype(0, 0, -1, -2));
            makeBad(bads).setGt(new MongoGenotype(0, 0, -2, -1));
            return bads;
        } catch ( CloneNotSupportedException e ) {
            throw new ReviewedStingException("Failed to make BADMVCs", e);
        }
    }

    @DataProvider(name = "BadMVCs")
    public Object[][] makeBadMVCsProvider() throws CloneNotSupportedException {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final MongoVariantContext bad : makeBadMVCs() )
            tests.add(new Object[]{bad});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BadMVCs", expectedExceptions = MongoVariantContextException.class)
    public void testMVCValidate(final MongoVariantContext mvc) {
        mvc.validate(parser);
    }

    final static MongoVariantContext DUP = new MongoVariantContext(Arrays.asList("x"), "20", 1, 1, "A", "C", TruthStatus.TRUE_POSITIVE, new MongoGenotype(0, 1), new Date(), false);
    private static MongoVariantContext makeNoDup(final List<MongoVariantContext> nonDups) throws CloneNotSupportedException {
        final MongoVariantContext mvc = DUP.clone();
        nonDups.add(mvc);
        return mvc;
    }

    @DataProvider(name = "Duplicates")
    public Object[][] makeDuplicates() throws CloneNotSupportedException {
        List<Object[]> tests = new ArrayList<Object[]>();

        final List<MongoVariantContext> nonDups = new LinkedList<MongoVariantContext>();
        makeNoDup(nonDups).setSupportingCallSets(Arrays.asList("x", "y"));
        makeNoDup(nonDups).setSupportingCallSets(Arrays.asList("y"));
        makeNoDup(nonDups).setChr("21");
        makeNoDup(nonDups).setStart(2);
        makeNoDup(nonDups).setStop(2);
        makeNoDup(nonDups).setAlt("G");
        makeNoDup(nonDups).setTruth(TruthStatus.FALSE_POSITIVE);
        makeNoDup(nonDups).setReviewed(true);
        makeNoDup(nonDups).setGt(new MongoGenotype(0, 0));
        makeNoDup(nonDups).setGt(new MongoGenotype(1, 0));
        makeNoDup(nonDups).setGt(new MongoGenotype(1, 1));
        makeNoDup(nonDups).setGt(new MongoGenotype(-1, -1));
        makeNoDup(nonDups).setGt(new MongoGenotype(0, 0, 1, -1));
        makeNoDup(nonDups).setGt(new MongoGenotype(0, 0, -1, 1));

        tests.add(new Object[]{DUP, DUP, true});

        final MongoVariantContext mvc1 = DUP.clone();
        tests.add(new Object[]{DUP, mvc1, true});

        final MongoVariantContext mvc2 = DUP.clone();
        Calendar cal = Calendar.getInstance();
        cal.setTime(mvc2.getDate());
        mvc2.setDate(cal.getTime());
        tests.add(new Object[]{DUP, mvc2, true});

        final MongoVariantContext mvc3 = DUP.clone();
        Calendar cal2 = Calendar.getInstance();
        cal2.setTime(mvc3.getDate());
        cal2.add(Calendar.HOUR, 1);
        mvc3.setDate(cal2.getTime());
        tests.add(new Object[]{DUP, mvc3, true});

        for ( final MongoVariantContext nondup : nonDups )
            tests.add(new Object[]{DUP, nondup, false});


        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "Duplicates")
    public void testDuplicates(final MongoVariantContext mvc1, final MongoVariantContext mvc2, final boolean isDup) {
        Assert.assertEquals(mvc1.isDuplicate(mvc2), isDup, "MVCs " + mvc1 + " is dup of " + mvc2 + " returned expected value");
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
        try {
            setupBeforeMethod();

            final List<Allele> alleles = Arrays.asList(Allele.create("A", true), Allele.create("C"));
            final Genotype actualGT = db.consensusGT(truthStatus, PolymorphicStatus.POLYMORPHIC, alleles, individualCalls);

            if ( expectedGT.getGQ() == MongoGenotype.DISCORDANT_GQ )
                Assert.assertEquals(actualGT.getGQ(), MongoGenotype.DISCORDANT_GQ, "Expected GT was discordant but didn't get a discordant result");
            else
                Assert.assertEquals(
                        new MongoGenotype(alleles, actualGT),
                        new MongoGenotype(alleles, expectedGT),
                        "Failed to create expected consensus GT");
        }
        finally {
            teardownMethod();
        }
    }
}