package org.broadinstitute.sting.gatk.walkers.misc;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.GenotypeConcordance;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

/**
 * A walker for reassigning ref and alt alleles for an improperly-created VCF by using concordance metrics to some
 * orthogonal data (sequencing). If you feel the need to use this walker, contact the creator of your VCF and berate
 * him or her until he or she makes it correctly. Do not use this walker. In fact, the walker will try to stop you
 * if you are not me.
 */
@By(DataSource.REFERENCE_BASES)
public class FixAllelesByConcordance extends RodWalker<Integer,Integer> {

    @Input(fullName="chip",doc="You shouldn't be using this walker.",required=true)
    public RodBinding<VariantContext> chip = null;

    @Input(fullName="seq",doc="You shoulnd't be using this walker.",required=true)
    public RodBinding<VariantContext> seq = null;

    @Output
    public VCFWriter out;

    Set<String> samples;

    public void initialize() {
        /*if ( ! System.getProperty("user.name").equals("chartl") ) {
            class UserShouldNotIgnoreRepeatedWarningsException extends UserException {
                public UserShouldNotIgnoreRepeatedWarningsException(String msg) {
                    super(msg);
                }
            }

            throw new UserShouldNotIgnoreRepeatedWarningsException("You shouldn't be using this walker.");
        }*/

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit());
        Map<String,VCFHeader> onlyChip = new HashMap<String,VCFHeader>();
        onlyChip.put("chip",vcfRods.get("chip"));
        samples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.UNSORTED);
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);
        out.writeHeader(new VCFHeader(headerLines,samples));
    }

    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer a, Integer b) { return b + a;}

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) {
            return 0;
        }

        if ( ! tracker.hasValues(chip) ) {
            return 0;
        }

        Allele.create(ref.getBase(),true);

        VariantContext chipContext = tracker.getFirstValue(chip);
        if ( BaseUtils.basesAreEqual(chipContext.getReference().getBases()[0], ref.getBase()) ) {
            if ( VariantContextUtils.isTransversion(chipContext) ) {
                out.add(chipContext);
            } else {
                VariantContext revC = possiblyReverseComplement(chipContext, ref, tracker);
                if ( revC == null ) {
                    VariantContextBuilder b = new VariantContextBuilder(chipContext);
                    Set<String> filts = new HashSet<String>(chipContext.getFilters());
                    filts.add("NoComparisonForStrandCheck");
                    b.filters(filts);
                    out.add(b.make());
                    return 0;
                }
                out.add(revC);
            }
        } else {
            VariantContext chipFlip = flipRefAndAlt(chipContext,ref);
            if ( VariantContextUtils.isTransversion(chipFlip) ) {
                out.add(chipFlip);
            } else {
                VariantContext revC = possiblyReverseComplement(chipFlip, ref, tracker);
                if ( revC == null ) {
                    VariantContextBuilder b = new VariantContextBuilder(chipFlip);
                    Set<String> filts = new HashSet<String>(chipFlip.getFilters());
                    filts.add("NoComparisonForStrandCheck");
                    b.filters(filts);
                    out.add(b.make());
                    return 0;
                }
                out.add(revC);
            }
        }

        return 1;
    }

    private VariantContext flipRefAndAlt(VariantContext vc, ReferenceContext ref) {
        VariantContextBuilder vcb = new VariantContextBuilder(vc);
        List<Allele> newAlleles = new ArrayList<Allele>(2);
        Allele newRef = Allele.create(ref.getBase(),true);
        Allele newAlt = Allele.create(vc.getReference().getBases()[0],false);
        newAlleles.add(newRef);
        newAlleles.add(newAlt);
        List<Allele> HOM_REF = Arrays.asList(newRef,newRef);
        List<Allele> HOM_VAR = Arrays.asList(newAlt,newAlt);
        List<Allele> HET = Arrays.asList(newRef,newAlt);
        vcb.alleles(newAlleles);
        ArrayList<Genotype> newGenotypes = new ArrayList<Genotype>(vc.getNSamples());
        for ( Genotype g : vc.getGenotypes() ) {
            Genotype ng;

            if ( g.isHomRef() ) {
                ng = Genotype.modifyAlleles(g,HOM_VAR);
            } else if ( g.isHomVar() ) {
                ng = Genotype.modifyAlleles(g,HOM_REF);
            } else if ( g.isHet() ) {
                ng = Genotype.modifyAlleles(g,HET);
            } else {
                ng = g;
            }

            newGenotypes.add(ng);
        }

        vcb.genotypes(newGenotypes);

        return vcb.make();
    }

    private VariantContext possiblyReverseComplement(VariantContext vc, ReferenceContext ref, RefMetaDataTracker tracker) {
        if ( ! tracker.hasValues(seq) ) {
            return null;
        }
        VariantContext seqContext = tracker.getFirstValue(seq);
        double initConc = getHomConcordance(vc,seqContext);
        if ( initConc > 0.95 ) {
            return vc;
        }
        VariantContext revComp = reverseComplement(vc);
        double revConc = getHomConcordance(revComp,seqContext);
        if ( revConc > 0.95 ) {
            return revComp;
        }

        VariantContextBuilder b = new VariantContextBuilder(vc);
        Set<String> filt = new HashSet<String>(vc.getFilters());
        filt.add("PotentialStrandIssue(Bias)");
        b.filters(filt);
        return b.make();
    }

    private VariantContext reverseComplement(VariantContext vc) {
        // swap hom ref and hom var
        List<Allele> HOM_REF = Arrays.asList(vc.getReference(),vc.getReference());
        List<Allele> HOM_VAR = Arrays.asList(vc.getAlternateAllele(0),vc.getAlternateAllele(0));
        ArrayList<Genotype> newGenotypes = new ArrayList<Genotype>();
        for ( Genotype g : vc.getGenotypes() ) {
            Genotype ng;
            if ( g.isHomRef() ) {
                ng = Genotype.modifyAlleles(g,HOM_VAR);
            } else if ( g.isHomVar() ) {
                ng = Genotype.modifyAlleles(g,HOM_REF);
            } else {
                ng = g;
            }
            newGenotypes.add(ng);
        }
        VariantContextBuilder vcb = new VariantContextBuilder(vc);
        vcb.genotypes(newGenotypes);
        return vcb.make();
    }

    private double getHomConcordance(VariantContext v1, VariantContext v2) {
        int same = 0;
        int tot = 0;
        for ( String s : samples ) {
            Genotype g1 = v1.getGenotype(s);
            Genotype g2 = v2.getGenotype(s);
            if ( g1 != null && g2 != null && g1.isHom() && g2.isHom() ) {
                ++tot;
                same = same + (g1.getAllele(0).equals(g2.getAllele(0)) ? 1 : 0);
            }
        }
        return ( (double) same )/tot;
    }
}
