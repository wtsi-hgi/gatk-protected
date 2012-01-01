/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.*;

/**
 * Code to explore http://arxiv.org/abs/1112.1528
 *
 * @author Mark DePristo
 * @since 2011
 */
@Reference(window=@Window(start=-Yamagishi.WINDOW_SIZE,stop=Yamagishi.WINDOW_SIZE))
public class Yamagishi extends RodWalker<Pair<String, String>, Map<String, Yamagishi.KMerCounter>> {
    public final static int WINDOW_SIZE = 50;

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(shortName = "k", fullName = "kmerSize", doc="KMer size", required=true)
    protected int kmerSize;

    @Argument(shortName = "debug", fullName = "debug", doc="If true print out a lot of info", required=true)
    protected boolean DEBUG;

    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    @Argument(fullName="selectTypeToInclude", shortName="selectType", doc="Select only a certain type of variants from the input file. Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION. Can be specified multiple times", required=false)
    private List<VariantContext.Type> TYPES_TO_INCLUDE = Collections.emptyList();

    EnumSet<VariantContext.Type> includedTypes;

    @Override
    public void initialize() {
        if ( kmerSize != 3 )
            throw new UserException.BadArgumentValue("kmerSize", "Only k=3 is currently implemented");

        // if user specified types to include, add these, otherwise, add all possible variant context types to list of vc types to include
        includedTypes = TYPES_TO_INCLUDE.isEmpty() ? EnumSet.allOf(VariantContext.Type.class) : EnumSet.copyOf(TYPES_TO_INCLUDE);
    }

    private final List<String> ALL_NAMES = Arrays.asList("ref", "alt"); // , "snp", "indel");

    public class KMerCounter {
        final int k;
        final Map<String, Integer> counts = new HashMap<String, Integer>();
        final String name;

        public KMerCounter(final int k, final String name) {
            this.k = k;
            this.name = name;
        }

        public void inc(final String bases) {
            for ( int i = 0; i < bases.length() - k; i++ ) {
                final String kmer = bases.substring(i, i+k);
                final int prev = counts.containsKey(kmer) ? counts.get(kmer) : 0;
                counts.put(kmer, prev + 1);
            }
        }

        public void add(final KMerCounter o) {
            for ( Map.Entry<String, Integer> entry : o.counts.entrySet() ) {
                final String key = entry.getKey();
                final int by = entry.getValue();
                final int prev = counts.containsKey(key) ? counts.get(key) : 0;
                counts.put(key, prev + by);
            }
        }

        public Iterable<Map.Entry<String, Integer>> countsInOrder() {
            return new TreeMap<String, Integer>(counts).entrySet();
        }
    }

    @Override
    public Pair<String, String> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return null;

        // snp or indel if possible
        final VariantContext vc = tracker.getFirstValue(variantCollection.variants, ref.getLocus());

        if ( vc == null || ! vc.isBiallelic() || ! includedTypes.contains(vc.getType())) {
            return null;
        }

        final String refBases = new String(ref.getBases());
        final String leftRef = refBases.substring(0, WINDOW_SIZE);
        String vcBases = "";
        String rightRef = refBases.substring(WINDOW_SIZE+1, 2*WINDOW_SIZE+1);

        if ( vc.isSNP() ) {
            vcBases = vc.getAlternateAllele(0).getBaseString();
        } else if (vc.isIndel()) {
            vcBases = new String(new byte[]{vc.getReferenceBaseForIndel()});
            if ( vc.isSimpleInsertion() ) {
                vcBases = vcBases + vc.getAlternateAllele(0).getBaseString();
            } else if ( vc.isSimpleDeletion() ) {
                final int delSize = Math.abs(vc.getIndelLengths().get(0));
                if ( delSize > WINDOW_SIZE )
                    return null;
                rightRef = refBases.substring(WINDOW_SIZE+delSize+1, 2*WINDOW_SIZE+1);
            } else {
                return null;
                //throw new ReviewedStingException("What the hell is this VCF: " + vc);
            }
        } else {
            return null;
        }

        final String altBases = leftRef + vcBases + rightRef;

        if ( DEBUG ) {
            final String offset = Utils.dupString(' ', leftRef.length());
            logger.info("vc.ref  : " + vc.getReference());
            logger.info("vc.alt  : " + vc.getAlternateAllele(0));
            logger.info("left    : " + leftRef);
            logger.info("vc      : " + offset + vcBases);
            logger.info("right   : " + offset + rightRef);
            logger.info("refBases: " + refBases);
            logger.info("alt     : " + altBases);
        }

        return new Pair<String, String>(refBases, altBases);
    }

    @Override
    public Map<String, KMerCounter> reduceInit() {
        Map<String, KMerCounter> tables = new HashMap<String, KMerCounter>();
        for ( String name : ALL_NAMES )
            tables.put(name, new KMerCounter(kmerSize, name));
        return tables;
    }

    @Override
    public Map<String, KMerCounter> reduce(Pair<String, String> refAndAlt, Map<String, KMerCounter> sum) {
        if ( refAndAlt != null ) {
            sum.get("ref").inc(refAndAlt.getFirst());
            sum.get("alt").inc(refAndAlt.getSecond());
        }
        return sum;
    }

    @Override
    public void onTraversalDone(Map<String, KMerCounter> sum) {
        out.printf("name\tkmer\tcount%n");
        for ( KMerCounter table : new TreeMap<String, KMerCounter>(sum).values() ) {
            for ( Map.Entry<String, Integer> entry : table.countsInOrder() ) {
                out.printf("%s\t%s\t%d%n", table.name, entry.getKey(), entry.getValue());
            }
        }
    }
}