package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.testng.Assert;
import org.testng.SkipException;
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

    protected void teardownMethod() {
        if ( db != null ) {
            db.delete();
            db.close();
            db = null;
        }
    }
}