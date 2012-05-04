/*
 * Copyright (c) 2012, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.bcf2;

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Collection;

public class BCF2ToVCFWalker extends RodWalker<Integer, Integer> {

    @Input(fullName="bcf", shortName="bcf", doc="BCF2 file to convert", required=true)
    private RodBinding<VariantContext> bcf;

    @Output(doc="VCF File to create",required=true)
    private VCFWriter vcfWriter = null;

    @Override
    public void initialize() {
        VCFHeader bcfHeader = null;

        for ( ReferenceOrderedDataSource source : getToolkit().getRodDataSources() ) {
            if ( source.getName().equals(bcf.getName()) && source.getRecordType().equals(VariantContext.class) ) {
                bcfHeader = (VCFHeader)source.getHeader();
            }
        }

        if ( bcfHeader == null ) {
            throw new UserException("Failed to read BCF2 header");
        }

        // TODO: prevent VCF writer from adding duplicate contig declarations in the header
        vcfWriter.writeHeader(bcfHeader);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) {
            return 0;
        }

        Collection<VariantContext> vcs = tracker.getValues(bcf, context.getLocation());

        for ( VariantContext vc : vcs ) {
            vcfWriter.add(vc);
        }

        return vcs.size();
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce( Integer value, Integer sum ) {
        return value + sum;
    }

    @Override
    public void onTraversalDone( Integer result ) {
        super.onTraversalDone(result);
    }
}
