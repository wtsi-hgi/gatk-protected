package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import net.sf.picard.reference.FastaSequenceFile;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.io.File;

public class NA12878KBMain {
    public static void main(final String[] args) throws Exception {
        final FastaSequenceFile fasta = new FastaSequenceFile(new File(args[0]), false);
        final GenomeLocParser parser = new GenomeLocParser(fasta.getSequenceDictionary());
        final NA12878DBArgumentCollection dbArgumentCollection = new NA12878DBArgumentCollection();
        final NA12878KnowledgeBase db = new NA12878KnowledgeBase(parser, dbArgumentCollection);

        System.out.printf("Printing consensuss%n");
        int n = 0;
        for ( final MongoVariantContext mvc : db.getConsensusSites(new SiteSelector(parser))) {
            if ( n++ % 10 == 0 )
                System.out.printf("mvc + " + mvc.getChr() + ":" + mvc.getStart() + "%n");
            if ( n > 100 ) break;
        }

        db.close();
    }
}