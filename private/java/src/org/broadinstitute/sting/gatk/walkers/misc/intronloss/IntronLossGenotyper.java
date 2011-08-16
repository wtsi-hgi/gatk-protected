package org.broadinstitute.sting.gatk.walkers.misc.intronloss;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.refseq.RefSeqFeature;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * [Short one sentence description of this walker]
 *
 * <p>
 * [Functionality of this walker]
 * </p>
 *
 * <h2>Input</h2>
 * <p>
 * [Input description]
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * [Output description]
 * </p>
 *
 * <h2>Examples</h2>
 * PRE-TAG
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T $WalkerName
 * PRE-TAG
 *
 * @author Your Name
 * @since Date created
 */
public class IntronLossGenotyper extends LocusWalker<VariantCallContext,Integer> {
    // as a development note, this is implemented as a LocusWalker for now, and does a test at a specific locus,
    // in order to focus on genotype likelihood implementation and not on aggregation of information across
    // a window or a whole gene

    //@Argument(shortName="r",fullName="refSeq",required=true,doc="The RefSeq Gene definition track")
    public RodBinding<RefSeqFeature> refSeqRodBinding;

    @Argument(shortName="H",fullName="insertHistogram",required=true,doc="The insert size histogram per read group")
    public File readGroupInsertHistogram;

    @ArgumentCollection
    UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    @Output(doc="Output file")
    public VCFWriter writer;

    private UnifiedGenotyperEngine UG_engine = null;
    private IntronLossGenotypeLikelihoodCalculationModel ilglcm;

    private Map<String,byte[]> insertQualsByRG = new HashMap<String,byte[]>();

    public void initialize() {
        Set<Sample> samples = getToolkit().getSAMFileSamples();
        Set<String> sampleStr = new HashSet<String>(samples.size());
        for ( Sample s : samples ) {
            sampleStr.add(s.getId());
        }

        try {

            for ( String entry : new XReadLines(readGroupInsertHistogram) ) {
                String[] split1 = entry.split("\\t");
                String id = split1[0];
                String[] histogram = split1[1].split(";");
                byte[] quals = new byte[histogram.length];
                int idx = 0;
                for ( String histEntry : histogram ) {
                    quals[idx++] = Byte.parseByte(histEntry);
                }

                insertQualsByRG.put(id,quals);
            }

        } catch (FileNotFoundException e) {
            throw new UserException("Histogram file not found",e);
        }

        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();
        headerInfo.add(new VCFInfoHeaderLine(VCFConstants.ALLELE_FREQUENCY_KEY, -1, VCFHeaderLineType.Float, "Allele Frequency, for each ALT allele, in the same order as listed"));
        headerInfo.add(new VCFInfoHeaderLine(VCFConstants.ALLELE_COUNT_KEY, -1, VCFHeaderLineType.Integer, "Allele count in genotypes, for each ALT allele, in the same order as listed"));
        headerInfo.add(new VCFInfoHeaderLine(VCFConstants.ALLELE_NUMBER_KEY, 1, VCFHeaderLineType.Integer, "Total number of alleles in called genotypes"));
        headerInfo.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
        headerInfo.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1, VCFHeaderLineType.Float, "Genotype Quality"));
        headerInfo.add(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Read Depth (only filtered reads used for calling)"));
        headerInfo.add(new VCFFormatHeaderLine(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, 3, VCFHeaderLineType.Float, "Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic"));


        // initialize the header
        writer.writeHeader(new VCFHeader(headerInfo, sampleStr));

        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, sampleStr);
        ilglcm = new IntronLossGenotypeLikelihoodCalculationModel(refSeqRodBinding,getToolkit().getGenomeLocParser(),insertQualsByRG);
    }

    public Integer reduceInit() { return null; }

    public boolean filter(SAMRecord record) {
        // todo -- for now filtering unmapped reads, they should actually be an event of their own
        if ( record.getMateUnmappedFlag() ) {
            return false;
        }

        transform(record);
        return true;
    }

    public VariantCallContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
         if ( tracker == null )
            return null;
        if ( ! tracker.hasValues(refSeqRodBinding) ) {
            return null;
        }
        RefSeqFeature refSeq = tracker.getFirstValue(refSeqRodBinding);
        ilglcm.resetRefSeqData(refSeq,ref.getLocus());

        // initialize the genotype likelihoods
        ArrayList<Allele> refAltAlleles = new ArrayList<Allele>(2);
        refAltAlleles.add(Allele.create(ref.getBase()));
        refAltAlleles.add(Allele.NO_CALL);
        Map<String,MultiallelicGenotypeLikelihoods> genotypeLikelihoods = new HashMap<String,MultiallelicGenotypeLikelihoods>();
        Map<String,AlignmentContext> splitContext = AlignmentContextUtils.splitContextBySampleName(context);
        for ( Map.Entry<String,AlignmentContext> entry : splitContext.entrySet() ) {
            genotypeLikelihoods.put(entry.getKey(),new MultiallelicGenotypeLikelihoods(entry.getKey(),refAltAlleles,new double[IntronLossGenotypeLikelihoodCalculationModel.PLOIDY],0));
        }

        ilglcm.getLikelihoods(tracker,ref, splitContext,
                AlignmentContextUtils.ReadOrientation.COMPLETE,null,genotypeLikelihoods,Allele.NO_CALL,false);

        // convert likelihoods into genotypes
        Map<String,Genotype> sampleGenotypes = new HashMap<String,Genotype>(genotypeLikelihoods.size());
        for ( Map.Entry<String,MultiallelicGenotypeLikelihoods> entry : genotypeLikelihoods.entrySet() ) {
            HashMap<String,Object> attributes = new HashMap<String,Object>(1);
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY,entry.getValue());
            Genotype g = new Genotype(entry.getKey(), refAltAlleles, Genotype.NO_NEG_LOG_10PERROR, null, attributes, false);
            sampleGenotypes.put(entry.getKey(),g);
        }

        // todo -- when generalizing ploidy here, switch variant context to calling ExactAF directly
        VariantContext vc = new VariantContext("ILG",ref.getLocus().getContig(),ref.getLocus().getStart(),ref.getLocus().getStop(),refAltAlleles,sampleGenotypes,VariantContext.NO_NEG_LOG_10PERROR,null,null);
        return UG_engine.calculateGenotypes(tracker, ref, context, vc);
    }

    public Integer reduce(VariantCallContext value, Integer sum) {
        if ( value == null )
            return sum;

        try {
            Map<String, Object> attrs = new HashMap<String, Object>(value.getAttributes());
            VariantContextUtils.calculateChromosomeCounts(value, attrs, true);
            writer.add(VariantContext.modifyAttributes(value, attrs));
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(e.getMessage() + "; this is often caused by using the --assume_single_sample_reads argument with the wrong sample name");
        }

        return sum + 1;
    }


    private void transform(SAMRecord read) {
        // add in read information (p[insert size] etc)
    }
}
