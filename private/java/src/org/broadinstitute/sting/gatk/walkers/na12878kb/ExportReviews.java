package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;

public class ExportReviews extends NA12878DBWalker {
    @Output
    public VariantContextWriter out;

    public void initialize() {
        super.initialize();

        final VCFHeader header = new VCFHeader();

        for ( final CallSet callSet : db.getCallSets() ) {
            if ( callSet.isReviewer() )
                header.addMetaDataLine(callSet.asVCFHeaderLine());
        }

        for ( final VCFHeaderLine line : MongoVariantContext.reviewHeaderLines() )
            header.addMetaDataLine(line);

        out.writeHeader(header);
    }

    @Override public boolean isDone() { return true; }

    @Override
    public void onTraversalDone(Integer result) {
        db.writeReviews(out, super.makeSiteSelector());
        super.onTraversalDone(result);
    }
}