package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.AfterSuite;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class NA12878KBUnitTestBase extends BaseTest {
    private final static Logger logger = Logger.getLogger(NA12878KBUnitTestBase.class);
    protected final static String NA12878_KB_TESTFILES = privateTestDir + "/na12878kb/";
    protected final static File testVCF = new File(NA12878_KB_TESTFILES + "test.vcf");
    protected final static List<VariantContext> testVCs;
    protected final static GenomeLocParser parser;

    protected NA12878KnowledgeBase db;

    static {
        try {
            testVCs = Collections.unmodifiableList(VCFUtils.readVCF(testVCF).getSecond());
            final IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(new File(b37KGReference));
            parser = new GenomeLocParser(fasta.getSequenceDictionary());
        } catch ( IOException e ) {
            throw new SkipException("Couldn't read test VCF so skipping NA12878KB unit tests", e);
        }
    }

    @BeforeSuite
    protected void setupBeforeMethod() {
        try {
            final NA12878DBArgumentCollection args = new NA12878DBArgumentCollection();
            //args.useLocal = true;
            args.dbToUse = NA12878DBArgumentCollection.DBType.TEST;
            args.resetDB = true;
            db = new NA12878KnowledgeBase(parser, args);
        } catch ( Exception e ) {
            throw new SkipException("Failed to setup DB connection", e);
        }
    }

    @AfterSuite
    protected void teardownMethod() {
        if ( db != null )
            db.close();
    }

    @DataProvider(name = "TestVCProvider")
    public Object[][] testVCProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( final VariantContext vc : testVCs )
            tests.add(new Object[]{vc});
        return tests.toArray(new Object[][]{});
    }

    @Test
    final MongoVariantContext readOneMVCFromDB() {
        final SiteIterator<MongoVariantContext> it = db.getCalls();
        Assert.assertTrue(it.hasNext(), "Expected at least 1 call in the db");
        final MongoVariantContext mvc = it.next();
        Assert.assertFalse(it.hasNext(), "Only expected 1 call in the db but saw > 1");
        it.close();
        return mvc;
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
        final Genotype actualGT = db.consensusGT(truthStatus, PolymorphicStatus.POLYMORPHIC, alleles, individualCalls);

        if ( expectedGT.getGQ() == MongoGenotype.DISCORDANT_GQ )
            Assert.assertEquals(actualGT.getGQ(), MongoGenotype.DISCORDANT_GQ, "Expected GT was discordant but didn't get a discordant result");
        else
            Assert.assertEquals(
                    new MongoGenotype(alleles, actualGT),
                    new MongoGenotype(alleles, expectedGT),
                    "Failed to create expected consensus GT");
    }
}