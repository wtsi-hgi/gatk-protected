package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.gatk.walkers.na12878kb.errors.MongoVariantContextException;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextTestProvider;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
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

    @BeforeMethod
    public void setup() {
        setupBeforeMethod();
    }

    @AfterMethod
    public void teardown() {
        teardownMethod();
    }

    @Test(dataProvider = "TestVCProvider")
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

    @Test(dataProvider = "TestVCProvider")
    public void testMVCMatching(final VariantContext vc) {
        final MongoVariantContext mvc = MongoVariantContext.create("x", vc, TruthStatus.UNKNOWN, MongoGenotype.NO_CALL);
        for ( final VariantContext vc2 : testVCs ) {
            final boolean expectedEquals = vc == vc2;
            Assert.assertEquals(mvc.matches(vc2), expectedEquals, "MVC match returned unexpected value for mvc=" + mvc + " vs vc2=" + vc2);
        }
    }

    @Test(dataProvider = "MVCBasicTest")
    public void testMVCBasic(final VariantContext originalVC, final MongoVariantContext mvc) {
        db.addCall(mvc);
        final MongoVariantContext fromDB = readOneMVCFromDB();

        Assert.assertEquals(fromDB, mvc, "Input MongoVariantContext not the same as the one read from DB");
        VariantContextTestProvider.assertEquals(fromDB.getVariantContext(), mvc.getVariantContext());
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

}