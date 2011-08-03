/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.*;

import static org.broadinstitute.sting.utils.codecs.vcf.VCFUtils.getVCFHeadersFromRods;

/**
 * Walks along all variant ROD loci and annotates inherited alleles .
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})

@ReadFilters({MappingQualityZeroReadFilter.class})
// Filter out all reads with zero mapping quality

public class AnnotateTrioPhasingInheritanceNoRecombinationWalker extends RodWalker<Integer, Integer> {
    public final static String TRIO_ROD_NAME = "trio";

    private final static int DIPLOID = 2;

    public static final String INHERITANCE_KEY = "HP";

    @Output(doc = "File to which trio-phased variants should be written", required = true)
    protected VCFWriter writer = null;

    @Argument(shortName = "f", fullName = "familyPattern", required = true, doc = "Pattern for the family structure (usage: mom+dad=child)")
    public String familyStr = null;

    private String SAMPLE_NAME_MOM;
    private String SAMPLE_NAME_DAD;
    private String SAMPLE_NAME_CHILD;

    public void initialize() {
        initializeVcfWriter();
    }

    private void initializeVcfWriter() {
        String[] pieces = familyStr.split("[\\+\\=]");

        SAMPLE_NAME_MOM = pieces[0];
        SAMPLE_NAME_DAD = pieces[1];
        SAMPLE_NAME_CHILD = pieces[2];

        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(TRIO_ROD_NAME);

        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(rodNameToHeader, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        if (vcfSamples.size() != 3) {
            throw new UserException("File to annotate with trio phasing contains more than three samples.  This walker only" +
                    "accepts VCFs with three samples.");
        }

        if (!vcfSamples.contains(SAMPLE_NAME_MOM) || !vcfSamples.contains(SAMPLE_NAME_DAD) || !vcfSamples.contains(SAMPLE_NAME_CHILD)) {
            throw new UserException("One or more of the samples specified in the familyPattern argument is not present" +
                    "in this file.  Please supply a VCF file that contains only three samples: the" +
                    "mother, the father, and the child");
        }

        Set<String> samples = new HashSet<String>();
        samples.add(SAMPLE_NAME_MOM);
        samples.add(SAMPLE_NAME_DAD);
        samples.add(SAMPLE_NAME_CHILD);

        // setup the header fields:
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));
        hInfo.addAll(VCFUtils.getHeaderFields(this.getToolkit()));

        hInfo.add(new VCFFormatHeaderLine(INHERITANCE_KEY, 1, VCFHeaderLineType.String, "Source of inherited allele"));

        writer.writeHeader(new VCFHeader(hInfo, samples));
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public Integer reduceInit() {
        return 0;
    }

    /**
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        GenomeLoc loc = ref.getLocus();
        VariantContext trioVc = tracker.getFirstValue(VariantContext.class, TRIO_ROD_NAME, loc);
        if (trioVc == null)
            return null;

        int annotated = 0;
        if (!trioVc.isFiltered()) {
            Genotype father = trioVc.getGenotype(SAMPLE_NAME_DAD);
            Genotype mother = trioVc.getGenotype(SAMPLE_NAME_MOM);
            Genotype child = trioVc.getGenotype(SAMPLE_NAME_CHILD);

            if (child.isHet() && father.isCalled() && father.isNotFiltered() && mother.isCalled() && mother.isNotFiltered() && (father.isHom() || mother.isHom())) {
                Genotype hom;
                Genotype other;
                if (father.isHom()) {
                    hom = father;
                    other = mother;
                }
                else if (mother.isHom()) {
                    hom = mother;
                    other = father;
                }
                else
                    throw new ReviewedStingException("LOGICAL ERROR: at least one parent is hom!");

                Allele homAllele = hom.getAllele(0);
                Set<Allele> otherAlleles = new TreeSet<Allele>(other.getAlleles());

                if (child.getPloidy() > DIPLOID)
                    throw new UserException("Can only trio phase DIPLOID genotypes!");

                String[] alleleSources = new String[child.getPloidy()];
                int ind = 0;
                for (Allele allele : child.getAlleles()) {
                    if (allele.equals(homAllele)) {
                        alleleSources[ind] = getParentTitle(hom.getSampleName());
                    }
                    else if (otherAlleles.contains(allele)) {
                        alleleSources[ind] = getParentTitle(other.getSampleName());
                    }
                    else {
                        logger.warn("CANNOT trio phase, due to de novo appearance of alleles at: " + VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), trioVc));
                        alleleSources = null;
                        break;
                    }
                    ind++;
                }

                if (alleleSources != null) {
                    StringBuffer sb = new StringBuffer();
                    for (int i = 0; i < alleleSources.length; i++) {
                        sb.append(alleleSources[i]);
                        if (i < alleleSources.length - 1)
                            sb.append("|");
                    }

                    Map<String, Object> childInfo = new HashMap<String, Object>(child.getAttributes());
                    childInfo.put(INHERITANCE_KEY, sb.toString());
                    child = new Genotype(child.getSampleName(), child.getAlleles(), child.getNegLog10PError(), child.getFilters(), childInfo, child.isPhased());

                    Map<String, Genotype> genotypes = trioVc.getGenotypes();
                    genotypes.put(SAMPLE_NAME_CHILD, child);
                    trioVc = VariantContext.modifyGenotypes(trioVc, genotypes);

                    annotated++;
                }
            }
        }

        WriteVCF.writeVCF(trioVc, writer, logger);

        return annotated;
    }

    public String getParentTitle(String sample) {
        if (sample.equals(SAMPLE_NAME_DAD))
            return "P";
        else if (sample.equals(SAMPLE_NAME_MOM))
            return "M";

        throw new ReviewedStingException("Logical error: should only pass father or mother");
    }

    public Integer reduce(Integer addIn, Integer runningCount) {
        if (addIn == null)
            addIn = 0;

        return runningCount + addIn;
    }

    /**
     * @param result the number of annotated sites.
     */
    public void onTraversalDone(Integer result) {
        System.out.println("Annotated " + result + " sites.");
    }
}