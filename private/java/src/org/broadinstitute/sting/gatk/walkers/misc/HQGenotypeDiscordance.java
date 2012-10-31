package org.broadinstitute.sting.gatk.walkers.misc;

import com.google.common.base.Function;
import com.google.common.collect.Collections2;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 9/11/12
 * Time: 9:49 AM
 * To change this template use File | Settings | File Templates.
 */
public class HQGenotypeDiscordance extends RodWalker<Integer,Long> {

    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public List<RodBinding<VariantContext>> variants;


    @Argument(doc="The genotype quality threshold",required=false,fullName="genotypeQualityThreshold")
    public int genotypeQualityThreshold = 50;

    @Output
    public PrintStream out;

    private Set<String> sampleOverlap;

    public void initialize() {
        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit());
        sampleOverlap = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.PRIORITIZE);
    }

    public Long reduceInit() {
        return 0L;
    }

    public Long reduce(Integer numDiscord, Long total) {
        return ( total + numDiscord );
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        logger.info("I am in map!!");
        int numDiscordant = 0;
        int numSamples = 0;
        if ( tracker != null ) {
            List<VariantContext> boundVCs = removeFiltered(tracker.getValues(variants));
            if ( boundVCs.size() > 1 ) {
                for ( final String sample : sampleOverlap ) {
                    Collection<Genotype> sampleGenotypes = extractGenotypes(boundVCs,sample);
                    Genotype previousGenotype = null;
                    for ( Genotype genotype : sampleGenotypes ) {
                        if ( previousGenotype != null ) {
                            if ( canUseGenotypes(previousGenotype,genotype) ) {
                                if ( genotypesDiscordant(previousGenotype,genotype)) {
                                    numDiscordant += 1;
                                }
                                numSamples += 1;
                            }
                        }

                        previousGenotype = genotype;
                    }
                }
            }
        }
        double discordantProp = numSamples == 0 ? 0.0 : ((double) numDiscordant)/numSamples;
        VariantContext vc = tracker.getValues(variants).get(0);
        logger.info(String.format("%s\t%d\t%d\t%d\t%.3f%n",vc.getChr(),vc.getStart(),numDiscordant,numSamples,discordantProp));
        return numDiscordant;
    }

    private Collection<Genotype> extractGenotypes(final List<VariantContext> variantContexts, final String sample) {
        return Collections2.transform(variantContexts, new Function<VariantContext, Genotype>() {
                    @Override
                    public Genotype apply(VariantContext variantContext) {
                        return variantContext.getGenotype(sample);  //To change body of implemented methods use File | Settings | File Templates.
                    }
                });
    }

    private List<VariantContext> removeFiltered(List<VariantContext> contexts) {
        ArrayList<VariantContext> cleanVCs = new ArrayList<VariantContext>(contexts.size());
        for ( VariantContext vc : contexts ) {
            if ( ! vc.isFiltered() ) {
                cleanVCs.add(vc);
            }
        }

        return cleanVCs;
    }

    private boolean canUseGenotypes(Genotype a, Genotype b) {
        return ! a.isNoCall() && ! b.isNoCall() && a.hasGQ() && b.hasGQ() &&
                a.getGQ() > genotypeQualityThreshold && b.getGQ() > genotypeQualityThreshold;
    }

    private boolean genotypesDiscordant(Genotype a, Genotype b) {
        return a.isHomRef() == b.isHomRef() && a.isHet() == b.isHet();
    }
}
