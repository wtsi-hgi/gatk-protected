/*
 * Copyright (c) 2011 The Broad Institute
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.fragments.FragmentCollection;
import org.broadinstitute.sting.utils.fragments.FragmentUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Call SNPs and indels simultaneously via local de-novo assembly. Haplotype likelihoods are evaluated with an affine gap penalty HMM.
 *
 * <h2>Input</h2>
 * <p>
 * None! Although on-the-fly base quality score recalibration tables are highly recommended.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * VCF file with raw, unrecalibrated SNP and indel calls.
 * </p>
 *
 * <h2>Examples</h2>
 * PRE-TAG
 *   java
 *     -jar GenomeAnalysisTK.jar
 *     -T HaplotypeCaller
 *     -R reference/human_g1k_v37.fasta
 *     -I input.bam
 *     -BQSR input.kmer.8.recal_data.csv [[optional, but recommended]]
 *     -o output.raw.snps.indels.vcf
 * PRE-TAG
 *
 * @author rpoplin
 * @since 8/22/11
 */

@PartitionBy(PartitionType.INTERVAL)
@ActiveRegionExtension(extension=30)
public class HaplotypeCaller extends ActiveRegionWalker<Integer, Integer> {

    /**
     * A raw, unfiltered, highly specific callset in VCF format.
     */
    @Output(doc="File to which variants should be written", required = true)
    protected VCFWriter vcfWriter = null;

    @Output(fullName="graphOutput", shortName="graph", doc="File to which debug assembly graph information should be written", required = false)
    protected PrintStream graphWriter = null;

    @Argument(fullName = "assembler", shortName = "assembler", doc = "Assembler to use; currently only SIMPLE_DE_BRUIJN is available.", required = false)
    protected LocalAssemblyEngine.ASSEMBLER ASSEMBLER_TO_USE = LocalAssemblyEngine.ASSEMBLER.SIMPLE_DE_BRUIJN;

    @Argument(fullName="gopHMM", shortName="gopHMM", doc="gopHMM", required = false)
    protected double gopHMM = 45.0;

    @Argument(fullName="gcpHMM", shortName="gcpHMM", doc="gcpHMM", required = false)
    protected double gcpHMM = 10.0;

    @Argument(fullName="gopSW", shortName="gopSW", doc="gopSW", required = false)
    protected double gopSW = 30.0;

    @Argument(fullName="gcpSW", shortName="gcpSW", doc="gcpSW", required = false)
    protected double gcpSW = 1.4;

    @Argument(fullName="downsampleRegion", shortName="dr", doc="coverage per sample to downsample each region to", required = false)
    protected int DOWNSAMPLE_PER_SAMPLE_PER_REGION = 1000;

    @ArgumentCollection
    private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // the calculation arguments
    private UnifiedGenotyperEngine UG_engine = null;
    private UnifiedGenotyperEngine UG_engine_simple_genotyper = null;
    
    @Argument(fullName="debug", shortName="debug", doc="If specified print out very verbose debug information about each triggering interval", required = false)
    protected boolean DEBUG;

    // the assembly engine
    LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    LikelihoodCalculationEngine likelihoodCalculationEngine = null;

    // the genotyping engine
    GenotypingEngine genotypingEngine = null;

    // fasta reference reader to supplement the edges of the reference sequence
    private IndexedFastaSequenceFile referenceReader;

    // reference base padding size
    private static final int REFERENCE_PADDING = 120;

    // bases with quality less than or equal to this value are trimmed off the tails of the reads
    private static final byte MIN_TAIL_QUALITY = 8;

    private ArrayList<String> samplesList = new ArrayList<String>();
    private final static double LOG_ONE_HALF = -Math.log10(2.0);
    private final static double LOG_ONE_THIRD = -Math.log10(3.0);
    private boolean GENOTYPE_GIVEN_ALLELES_MODE = false;
    private final ArrayList<VariantContext> allelesToGenotype = new ArrayList<VariantContext>();

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        super.initialize();

        // get all of the unique sample names
        Set<String> samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        samplesList.addAll( samples );
        // initialize the UnifiedGenotyper Engine which is used to call into the exact model
        UAC.MAX_ALTERNATE_ALLELES = 4;
        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, samples);
        UAC.STANDARD_CONFIDENCE_FOR_CALLING = 0.01; // low values used for isActive determination only, default/user-specified values used for actual calling
        UAC.STANDARD_CONFIDENCE_FOR_EMITTING = 0.01; // low values used for isActive determination only, default/user-specified values used for actual calling
        UG_engine_simple_genotyper = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, samples);
        // initialize the output VCF header
        vcfWriter.writeHeader(new VCFHeader(new HashSet<VCFHeaderLine>(), samples));

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile, e);
        }

        assemblyEngine = makeAssembler(ASSEMBLER_TO_USE, referenceReader);
        likelihoodCalculationEngine = new LikelihoodCalculationEngine(gopHMM, gcpHMM, DEBUG, true, false);
        genotypingEngine = new GenotypingEngine( DEBUG, gopSW, gcpSW );
        if( UAC.alleles != null ) { GENOTYPE_GIVEN_ALLELES_MODE = true; }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // isActive
    //
    //---------------------------------------------------------------------------------------------------------------

    // enable deletions in the pileup
    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable non primary reads in the active region
    @Override
    public boolean wantsNonPrimaryReads() { return true; }

    @Override
    public double isActive( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {
        if( GENOTYPE_GIVEN_ALLELES_MODE ) {
            allelesToGenotype.addAll(tracker.getValues(UAC.alleles));
            return ( tracker.getValues(UAC.alleles).size() == 0 ? 0.0 : 1.0 ); 
        }
        if( context == null ) { return 0.0; }

        final List<Allele> noCall = new ArrayList<Allele>(); // used to noCall all genotypes until the exact model is applied
        noCall.add(Allele.NO_CALL);

        final Map<String, AlignmentContext> splitContexts = AlignmentContextUtils.splitContextBySampleName(context);
        final GenotypesContext genotypes = GenotypesContext.create(splitContexts.keySet().size());
        for( String sample : splitContexts.keySet() ) {
            final double[] genotypeLikelihoods = new double[3]; // ref versus non-ref (any event)
            genotypeLikelihoods[0] = 0.0;
            genotypeLikelihoods[1] = 0.0;
            genotypeLikelihoods[2] = 0.0;

            for( PileupElement p : context.getBasePileup() ) {
                byte qual = p.getQual();
                if (qual > QualityUtils.MIN_USABLE_Q_SCORE ) {
                    int AA = 0; int AB = 1; int BB = 2;
                    if( p.getBase() != ref.getBase() || p.isDeletion() || p.isBeforeInsertion() || p.isNextToSoftClip() || (p.getRead().getReadPairedFlag() && p.getRead().getMateUnmappedFlag()) || BadMateFilter.hasBadMate(p.getRead()) ) {
                        AA = 2;
                        BB = 0;
                        qual = (byte) ((int)qual + 6); // be overly permissive so as to not miss any slight signal of variation
                    } 
                    genotypeLikelihoods[AA] += Math.log10(QualityUtils.qualToProb(qual));
                    genotypeLikelihoods[AB] += MathUtils.softMax(Math.log10(QualityUtils.qualToProb(qual)) + LOG_ONE_HALF, Math.log10(QualityUtils.qualToErrorProb(qual)) + LOG_ONE_THIRD + LOG_ONE_HALF);
                    genotypeLikelihoods[BB] += Math.log10(QualityUtils.qualToErrorProb(qual)) + LOG_ONE_THIRD;
                }
            }

            final HashMap<String, Object> attributes = new HashMap<String, Object>();
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods((genotypeLikelihoods)));
            genotypes.add(new Genotype(sample, noCall, Genotype.NO_LOG10_PERROR, null, attributes, false));
        }

        final ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add( Allele.create("A", true) ); // fake ref Allele
        alleles.add( Allele.create("G", false) ); // fake alt Allele
        final VariantCallContext vcOut = UG_engine_simple_genotyper.calculateGenotypes(new VariantContextBuilder("HCisActive!", context.getContig(), context.getLocation().getStart(), context.getLocation().getStop(), alleles).genotypes(genotypes).make(), GenotypeLikelihoodsCalculationModel.Model.SNP);
        return ( vcOut == null ? 0.0 : QualityUtils.qualToProb( vcOut.getPhredScaledQual() ) );
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Integer map( final ActiveRegion activeRegion, final RefMetaDataTracker metaDataTracker ) {
        if( GENOTYPE_GIVEN_ALLELES_MODE ) {
            final ArrayList<VariantContext> vcsToRemove = new ArrayList<VariantContext>();
            for( final VariantContext vc : allelesToGenotype ) {
                if( activeRegion.getLocation().overlapsP( getToolkit().getGenomeLocParser().createGenomeLoc(vc) ) ) {
                    vcsToRemove.add(vc); // BUGBUG: do something with these VCs during GGA
                }
            }
            allelesToGenotype.removeAll( vcsToRemove );
        }
        
        if( !activeRegion.isActive ) { return 0; } // Not active so nothing to do!
        if ( activeRegion.size() == 0 ) { return 0; } // No reads here so nothing to do!

        if( DEBUG ) { System.out.println("Assembling " + activeRegion.getLocation() + " with " + activeRegion.size() + " reads:"); }
        finalizeActiveRegion( activeRegion ); // merge overlapping fragments, clip adapter and low qual tails
        final ArrayList<Haplotype> haplotypes = assemblyEngine.runLocalAssembly( activeRegion.getReads(), new Haplotype(activeRegion.getReference(referenceReader) ) );
        if( DEBUG ) { System.out.println("Found " + haplotypes.size() + " candidate haplotypes to evaluate"); }
        genotypingEngine.createEventDictionaryAndFilterBadHaplotypes( haplotypes, activeRegion.getReference(referenceReader, REFERENCE_PADDING), getPaddedLoc( activeRegion ), activeRegion.getLocation() );
        if( DEBUG ) { System.out.println(haplotypes.size() + " candidate haplotypes remain after filtering"); }

        if( haplotypes.size() == 0 ) {
            if( DEBUG ) { System.out.println("WARNING! No haplotypes created during assembly!"); }
            return 0;
        }

        final HashMap<String, Double[][]> haplotypeLikehoodMatrixMap = new HashMap<String, Double[][]>();
        final ArrayList<Haplotype> bestTwoHaplotypesPerSample = new ArrayList<Haplotype>();
        filterNonPassingReads( activeRegion );

        final HashMap<String, ArrayList<GATKSAMRecord>> readListMap = splitReadsBySample( activeRegion.getReads() );

        for( final String sample : readListMap.keySet() ) {
            if( DEBUG ) { System.out.println("Evaluating sample " + sample + " with " + readListMap.get( sample ).size() + " passing reads"); }
            likelihoodCalculationEngine.computeLikelihoods( haplotypes, readListMap.get( sample ) );
            for( final Haplotype bestHaplotype : likelihoodCalculationEngine.chooseBestHaplotypes(haplotypes) ) {
                if( !bestTwoHaplotypesPerSample.contains( bestHaplotype ) ) {
                    bestTwoHaplotypesPerSample.add( bestHaplotype );
                }
            }
            haplotypeLikehoodMatrixMap.put( sample, likelihoodCalculationEngine.haplotypeLikehoodMatrix );
        }

        final ArrayList<VariantContext> vcs = genotypingEngine.alignAndAssignGenotypeLikelihoods( getToolkit().getGenomeLocParser(), haplotypes, bestTwoHaplotypesPerSample, activeRegion.getReference(referenceReader, REFERENCE_PADDING), getPaddedLoc(activeRegion), activeRegion.getLocation(), haplotypeLikehoodMatrixMap );

        for( final VariantContext vc : vcs ) {
            if( activeRegion.getLocation().containsP(getToolkit().getGenomeLocParser().createGenomeLoc(vc).getStartLocation()) ) {
                final VariantCallContext vcOut = UG_engine.calculateGenotypes(vc, UG_engine.getUAC().GLmodel);
                if(vcOut != null) {
                    if( DEBUG ) { System.out.println(vcOut); }
                    vcfWriter.add(vcOut);
                }
            }
        }

        if( DEBUG ) { System.out.println("----------------------------------------------------------------------------------"); }

        return 1;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer cur, Integer sum) {
        return cur + sum;
    }

    @Override
    public void onTraversalDone(Integer result) {
        logger.info("Ran local assembly on " + result + " active regions");
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // private helper functions
    //
    //---------------------------------------------------------------------------------------------------------------

    private LocalAssemblyEngine makeAssembler(LocalAssemblyEngine.ASSEMBLER type, IndexedFastaSequenceFile referenceReader) {
        switch ( type ) {
            case SIMPLE_DE_BRUIJN:
                return new SimpleDeBruijnAssembler(graphWriter, referenceReader);
            default:
                throw new UserException.BadInput("Assembler type " + type + " is not valid/supported");
        }
    }

    private void finalizeActiveRegion( final ActiveRegion activeRegion ) {
        final ArrayList<GATKSAMRecord> finalizedReadList = new ArrayList<GATKSAMRecord>();
        final FragmentCollection<GATKSAMRecord> fragmentCollection = FragmentUtils.create( ReadUtils.sortReadsByCoordinate(activeRegion.getReads()) );
        activeRegion.clearReads();

        // Join overlapping paired reads to create a single longer read
        finalizedReadList.addAll( fragmentCollection.getSingletonReads() );
        for( final List<GATKSAMRecord> overlappingPair : fragmentCollection.getOverlappingPairs() ) {
            finalizedReadList.addAll( FragmentUtils.mergeOverlappingPairedFragments(overlappingPair) );
        }

        Collections.shuffle(finalizedReadList);

        // Loop through the reads hard clipping the adaptor and low quality tails
        for( final GATKSAMRecord myRead : finalizedReadList ) {
            final GATKSAMRecord postAdapterRead = ( myRead.getReadUnmappedFlag() ? myRead : ReadClipper.hardClipAdaptorSequence( myRead ) );
            if( postAdapterRead != null && !postAdapterRead.isEmpty() && postAdapterRead.getCigar().getReadLength() > 0 ) {
                final GATKSAMRecord clippedRead = ReadClipper.hardClipLowQualEnds(postAdapterRead, MIN_TAIL_QUALITY );

                // protect against INTERVALS with abnormally high coverage
                if( clippedRead.getReadLength() > 0 && activeRegion.size() < samplesList.size() * DOWNSAMPLE_PER_SAMPLE_PER_REGION ) {
                    activeRegion.add(clippedRead);
                }
            }
        }
    }

    private void filterNonPassingReads( final ActiveRegion activeRegion ) {
        final ArrayList<GATKSAMRecord> readsToRemove = new ArrayList<GATKSAMRecord>();
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {
            if( rec.getMappingQuality() <= 18 || BadMateFilter.hasBadMate(rec) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll( readsToRemove );
    }

    private GenomeLoc getPaddedLoc( final ActiveRegion activeRegion ) {
        final int padLeft = Math.max(activeRegion.getReferenceLoc().getStart()-REFERENCE_PADDING, 1);
        final int padRight = Math.min(activeRegion.getReferenceLoc().getStop()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(activeRegion.getReferenceLoc().getContig()).getSequenceLength());
        return getToolkit().getGenomeLocParser().createGenomeLoc(activeRegion.getReferenceLoc().getContig(), padLeft, padRight);
    }
    
    private HashMap<String, ArrayList<GATKSAMRecord>> splitReadsBySample( final ArrayList<GATKSAMRecord> reads ) {
        final HashMap<String, ArrayList<GATKSAMRecord>> returnMap = new HashMap<String, ArrayList<GATKSAMRecord>>();
        for( final String sample : samplesList) {
            ArrayList<GATKSAMRecord> readList = returnMap.get( sample );
            if( readList == null ) {
                readList = new ArrayList<GATKSAMRecord>();
                returnMap.put(sample, readList);
            }
        }
        for( final GATKSAMRecord read : reads ) {
            final String sample = read.getReadGroup().getSample();
            ArrayList<GATKSAMRecord> readList = returnMap.get( sample );
            readList.add(read);
        }

        return returnMap;
    }
}