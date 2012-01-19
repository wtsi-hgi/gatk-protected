package org.broadinstitute.sting.gatk.walkers.misc;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/18/12
 * Time: 4:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class MultiplyLikelihoods extends RodWalker<Integer,Integer> {
    @Input(shortName="V",fullName="Variants",required = true,doc="A set of variant contexts to merge genotypes. Samples will be intersected.")
    List<RodBinding<VariantContext>> variants;

    @Argument(shortName="Q",fullName="ChipQual",required=false,doc="The chip genotype quality.")
    int chipQual = 30;

    @Output
    VCFWriter out;

    Set<String> sampleIntersection;

    double[] HOM_REF;
    double[] HET;
    double[] HOM_VAR;
    List<Allele> NO_CALL = new ArrayList<Allele>(Arrays.asList(new Allele[]{Allele.NO_CALL}));

    public void initialize() {
        List<String> rodNames = new ArrayList<String>(16);
        for ( RodBinding<VariantContext> vcrb : variants ) {
            rodNames.add(vcrb.getName());
        }
        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> samples = new HashSet<String>(3200);
        boolean init = false;
        for ( Map.Entry<String,VCFHeader> header : vcfRods.entrySet() ) {
            if ( ! init ) {
                samples.addAll(header.getValue().getGenotypeSamples());
                init = true;
            } else {
                Set<String> notPres = new HashSet<String>((header.getValue().getGenotypeSamples()));
                notPres.removeAll(samples);
                samples.removeAll(notPres);
            }
        }
        sampleIntersection = new HashSet<String>(samples);

        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);
        out.writeHeader(new VCFHeader(headerLines,samples));

        double p1 =-((double) chipQual)/10;
        double p2 = p1+p1;

        HOM_REF = new double[]{0,p1,p2};
        HET = new double[]{p1,0,p1};
        HOM_VAR = new double[]{p2,p1,0};
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer map, Integer red) {
        return red + map;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) {
            return 0;
        }
        VariantContextBuilder builder = null;
        GenotypesContext genotypes = null;
        Set<Allele> alleles = null;
        boolean first = true;
        if ( tracker.getValues(variants).size() > 2 ) {
            logger.debug("foo");
        }
        for ( VariantContext vc : tracker.getValues(variants) ) {
            if ( first ) {
                builder = new VariantContextBuilder(vc);
                genotypes = getGenotypes(vc,sampleIntersection);
                alleles = new HashSet<Allele>(vc.getAlleles());
                first = false;
            } else {
                // for all other contexts check if alleles match
                if ( ! alleles.containsAll(vc.getAlleles()) ) {
                    logger.warn("Alleles do not match between variant contexts at location "+ref.getLocus().toString());
                    continue;
                }
                // if so proceed in adding genotype likelihoods
                ArrayList<Genotype> newGeno = new ArrayList<Genotype>(genotypes.size());
                for ( Genotype geno : genotypes ) {
                    newGeno.add(addGenotypeLikelihoods(geno, vc.getGenotype(geno.getSampleName())));
                }
                genotypes = GenotypesContext.create(newGeno);
            }
        }

        if ( builder == null ) {
            return 0;
        }
        builder.genotypes(genotypes); // reset the genotypes
        builder.attributes(new HashMap<String,Object>()); // reset the attributes
        out.add(builder.make());
        return 1;
    }

    private GenotypesContext getGenotypes(VariantContext context, Set<String> samples) {
        GenotypesContext genotypes = context.getGenotypes(samples);
        ArrayList<Genotype> newG = new ArrayList<Genotype>(genotypes.size());
        for ( Genotype g : genotypes ) {
            newG.add(Genotype.modifyAttributes(Genotype.modifyAlleles(g,NO_CALL),new HashMap<String,Object>()));
        }
        return GenotypesContext.create(newG);
    }

    private Genotype addGenotypeLikelihoods(Genotype g1, Genotype g2) {
        if ( g1.isNoCall() && ! g1.hasLikelihoods() ) {
            return g2;
        } else if ( g2.isNoCall() && ! g2.hasLikelihoods() ) {
            return g1;
        }

        double[] l1 = getLikelihoods(g1);
        double[] l2 = getLikelihoods(g2);

        double[] l3 = l1.clone();

        for ( int o = 0; o < l1.length; o++ ) {
            l3[o] += l2[o];
        }

        // todo -- recalculate GQ if we want it
        List<Allele> al = new ArrayList<Allele>(2);
        // todo -- generalize for non bi-allelic
        al.add(Allele.NO_CALL);
        return new Genotype(g1.getSampleName(),al, Genotype.NO_LOG10_PERROR ,new HashSet<String>(),new HashMap<String,Object>(),g1.isPhased(),l3);
    }

    protected double[] getLikelihoods(Genotype g) {
        if ( g.hasLikelihoods() ) {
            return g.getLikelihoods().getAsVector();
        }

        if ( g.isHomRef() ) {
            return HOM_REF;
        } else if ( g.isHet() ) {
            return HET;
        } else {
            return HOM_VAR;
        }
    }
}
