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

package org.broadinstitute.sting.utils.bcf2;

import org.broad.tribble.Feature;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.*;
import java.util.*;

/**
 * Testing BCF2
 *
 * @author Mark DePristo
 * @since 2012
 */
public class BCF2TestWalker extends RodWalker<Integer, Integer> {
    /**
     * Variants from this VCF file are used by this tool as input.
     * The file must at least contain the standard VCF header lines, but
     * can be empty (i.e., no variants are contained in the file).
     */
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public List<RodBinding<VariantContext>> variants;

    @Output(doc="File to which results should be written",required=true)
    protected File bcfFile;

    private final List<VariantContext> vcs = new ArrayList<VariantContext>();
    protected OutputStream mapOut;
    protected SimpleBCFEncoder encoder = new SimpleBCFEncoder();

    @Override
    public void initialize() {
        try {
            mapOut = new FileOutputStream(bcfFile);
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(bcfFile, "bad user!");
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        try {
            for ( VariantContext vc : tracker.getValues(variants, context.getLocation())) {
                encoder.encode(vc, mapOut);
                vcs.add(vc);
            }
        } catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile(bcfFile, "bad user!");
        }

        return 1;
    }

    //
    // default reduce -- doesn't do anything at all
    //
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer counter, Integer sum) { return counter + sum; }

    public void onTraversalDone(Integer sum) {
        try {
            mapOut.close();

            // read in the BCF records
            SimpleBCFDecoder codec = new SimpleBCFDecoder();
            PositionalBufferedStream pbs = new PositionalBufferedStream(new FileInputStream(bcfFile));
            codec.readHeader(pbs);

            Iterator<VariantContext> it = vcs.iterator();
            while ( ! pbs.isDone() ) {
                Feature loc = codec.decodeLoc(pbs);
                VariantContext expected = it.next();

                System.out.printf("bcf = %s %d%n", loc.getChr(), loc.getStart());
                System.out.printf("vcf = %s %d%n", expected.getChr(), expected.getStart());
                System.out.printf("%n");
            }

        } catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile(bcfFile, "bad user!");
        }
    }
}
