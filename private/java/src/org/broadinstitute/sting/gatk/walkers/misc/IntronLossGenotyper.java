package org.broadinstitute.sting.gatk.walkers.misc;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.refseq.RefSeqFeature;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.MultiallelicGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
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
    }

    public Integer reduceInit() { return null; }

    public boolean filter(SAMRecord record) {
        transform(record);
        return true;
    }

    public VariantCallContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
         if ( tracker == null )
            return null;
        RefSeqFeature refSeq = tracker.getValues(refSeqRodBinding).get(0);
        int numExons = refSeq.getNumExons();
        int thisExonNum = refSeq.getExonNumbersInInterval(ref.getLocus()).get(0);
        GenomeLoc prevExonLoc = thisExonNum > 1 ? refSeq.getExonLocation(thisExonNum-1) : null;
        GenomeLoc nextExonLoc = thisExonNum < numExons ? refSeq.getExonLocation(thisExonNum+1) : null;
        GenomeLoc exonLoc = refSeq.getExonLocation(thisExonNum);

        VariantContext vc = calculateGenotypeLikelihoods(context,ref, exonLoc, prevExonLoc,nextExonLoc);

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

    private VariantContext calculateGenotypeLikelihoods(AlignmentContext context, ReferenceContext ref, GenomeLoc thisEx, GenomeLoc prev, GenomeLoc next ) {
        Map<String,MultiallelicGenotypeLikelihoods> likelihoods = new HashMap<String,MultiallelicGenotypeLikelihoods>(512);
        for ( PileupElement e : context.getBasePileup() ) {
            if ( ! likelihoods.containsKey(e.getRead().getReadGroup().getSample()) ) {
                likelihoods.put(e.getRead().getReadGroup().getSample(),new MultiallelicGenotypeLikelihoods(e.getRead().getReadGroup().getSample(), getFakeAlt(ref), new double[3],0));
            }

            updateLikelihoods(likelihoods.get(e.getRead().getReadGroup().getSample()).getLikelihoods(),e.getRead(),thisEx,prev,next);
            // todo :: probably need multinomial coefficients, which requires some restructuring
        }

        return createVariantContextFromLikelihoods(ref,getFakeAlt(ref).get(0),likelihoods);
    }

    private void updateLikelihoods(double[] lik, SAMRecord read, GenomeLoc thisEx, GenomeLoc prev, GenomeLoc next) {
        // three genotype states: FF FE EE. What does the read contribute to each of these?
        // depends on the state of the read, we have:
        // --- PRETEND THESE DON'T HAPPEN TODO -- MAKE THEM HAPPEN ---
        // I - I
        // X - Other Contig
        // X with clipping
        // X - Unmapped
        // --- PRETEND THESE DO HAPPEN ---
        // X - I
        // X ----- X
        // X - X

        GenomeLoc alignment = getReadMateAlignment(read);

        boolean jumpsP = prev != null && alignment.overlapsP(prev);
        boolean jumpsN = next != null && alignment.overlapsP(next);

        if ( jumpsP && jumpsN ) { throw new ReviewedStingException("Read pair cannot occupy three locations: spanning previous and next exon"); }

        if ( jumpsP || jumpsN ) {
            // note for now we assume anchor read is in the exon, so this can only be
            // X ------------ X
            int delSize = jumpsP ? thisEx.distance(prev) : thisEx.distance(next);

            // FF - P(X1 not good) + P(X2 not good) + P(good and good)*P(insert size)
            double x1ng = Math.log10(QualityUtils.qualToErrorProb((byte) read.getMappingQuality()));
            double x2ng = Math.log10(QualityUtils.qualToErrorProb(Byte.parseByte((String) read.getAttribute("MA"))));
            double x1gx2g = Math.log10(QualityUtils.qualToErrorProb((byte)read.getMappingQuality())) + Math.log10(QualityUtils.qualToErrorProb(Byte.parseByte((String) read.getAttribute("MQ"))));
            double mod1 = MathUtils.log10sumLog10(new double[]{x1ng,x2ng,x1gx2g+Math.log10(QualityUtils.qualToProb(getInsertQual(read)))});
            lik[0] += mod1;

            // EE - P(good and good)*P(insert size-intron size)
            double mod2 = x1gx2g + Math.log10(QualityUtils.qualToProb(getInsertQual(read, delSize)));
            lik[2] += mod2;

            // FE - mixed
            lik[1] += mod1 + mod2 + -2*Math.log10(2.0);

        } else if ( alignment.overlapsP(thisEx) ) {
            // note for now we assume anchor read is in the exon, so this can only be
            // X -- X :: effectively uninformative
            // FF - P(good good)
            // EF - P(good good)
            // EE - P(good good)
            // todo -- see notes below on the mapping quality
            double val = Math.log10(QualityUtils.qualToErrorProb((byte)read.getMappingQuality())) + Math.log10(QualityUtils.qualToErrorProb(Byte.parseByte((String) read.getAttribute("MQ"))));
            lik[0] += val;
            lik[1] += val;
            lik[2] += val;
        } else {
            // note for now we assume anchor read is in the exon, so this can only be
            // X -- I


            // FF - P(X good map, I good map)
            // todo -- P(X good map) != 1-prob(mapping qual) but something close. Also P(X good, I good) != P(X good)P(I good) due to paired alignments.
            double mod1 = Math.log10(QualityUtils.qualToErrorProb((byte)read.getMappingQuality())) + Math.log10(QualityUtils.qualToErrorProb(Byte.parseByte((String) read.getAttribute("MQ"))));

            // EE - P(X good map, I bad map) + P(X bad map, I bad map) = P(I bad map)
            // todo -- P(I bad map) != prob(Mate Mapping Qual) -- but something close. Perhaps do this the right way.
            double mod2 = Math.log10(QualityUtils.qualToProb(Byte.parseByte( (String) read.getAttribute("MQ"))));

            lik[0] += mod1;

            // EF // 1/2 of EE and FF models
            lik[1] += -2*Math.log10(2.0) + MathUtils.log10sumLog10(new double[]{mod1,mod2} );

            lik[2] += mod2; //Math.log10(QualityUtils.qualToProb(Byte.parseByte( (String) read.getAttribute("MQ"))));
        }
    }

    // stolen from UGE
    private VariantContext createVariantContextFromLikelihoods(ReferenceContext refContext, Allele refAllele, Map<String, MultiallelicGenotypeLikelihoods> GLs) {
        // no-call everyone for now
        List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        Set<Allele> alleles = new LinkedHashSet<Allele>();
        alleles.add(refAllele);
        boolean addedAltAlleles = false;

        HashMap<String, Genotype> genotypes = new HashMap<String, Genotype>();
        for ( MultiallelicGenotypeLikelihoods GL : GLs.values() ) {
            if ( !addedAltAlleles ) {
                addedAltAlleles = true;
                // ordering important to maintain consistency
                for (Allele a: GL.getAlleles()) {
                    alleles.add(a);
                }
            }

            HashMap<String, Object> attributes = new HashMap<String, Object>();
            //GenotypeLikelihoods likelihoods = new GenotypeLikelihoods(GL.getLikelihoods());
            GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(GL.getLikelihoods());
            attributes.put(VCFConstants.DEPTH_KEY, GL.getDepth());
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, likelihoods);

            genotypes.put(GL.getSample(), new Genotype(GL.getSample(), noCall, Genotype.NO_NEG_LOG_10PERROR, null, attributes, false));
        }

        GenomeLoc loc = refContext.getLocus();
        int endLoc = 1+loc.getStop();

        return new VariantContext("UG_call",
                loc.getContig(),
                loc.getStart(),
                endLoc,
                alleles,
                genotypes,
                VariantContext.NO_NEG_LOG_10PERROR,
                null,
                null,
                refContext.getBase());
    }

    private ArrayList<Allele> getFakeAlt(ReferenceContext ref) {
        ArrayList<Allele> aal = new ArrayList<Allele>(1);
        aal.add(Allele.create(ref.getBase() == BaseUtils.A ? BaseUtils.C : BaseUtils.A));
        return aal;
    }

    private GenomeLoc getReadMateAlignment(SAMRecord read) {
        return getToolkit().getGenomeLocParser().createGenomeLoc(read.getMateReferenceName(),read.getMateAlignmentStart(), read.getMateAlignmentStart()+read.getReadLength());
    }

    private byte getInsertQual(SAMRecord read, int shift) {
        int insert = Math.abs(read.getInferredInsertSize())-shift;
        byte[] quals = insertQualsByRG.get(read.getReadGroup().getId());
        if ( insert > quals.length ) {
            return (byte) 70;
        }

        return quals[insert];
    }

    private byte getInsertQual(SAMRecord read) { return getInsertQual(read,0); }
}
