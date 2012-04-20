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

import com.google.java.contract.Ensures;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
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
import org.broadinstitute.sting.utils.collections.Pair;
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
 * Call SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. Haplotypes are evaluated using an affine gap penalty Pair HMM.
 *
 * <h2>Input</h2>
 * <p>
 * Input bam file(s) from which to make calls
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * VCF file with raw, unrecalibrated SNP and indel calls.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *   java
 *     -jar GenomeAnalysisTK.jar
 *     -T HaplotypeCaller
 *     -R reference/human_g1k_v37.fasta
 *     -I input.bam
 *     -o output.raw.snps.indels.vcf
 * </pre>
 *
 * @author rpoplin
 * @since 8/22/11
 */

@PartitionBy(PartitionType.LOCUS)
@ActiveRegionExtension(extension=50,maxRegion=300)
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

    @Argument(fullName="keepRG", shortName="keepRG", doc="keepRG", required = false)
    protected String keepRG = null;
    
    @Argument(fullName="mnpLookAhead", shortName="mnpLookAhead", doc = "The number of bases to combine together to form MNPs out of nearby consecutive SNPs on the same haplotype", required = false)
    protected int MNP_LOOK_AHEAD = 0;

    @Argument(fullName="genotypeFullActiveRegion", shortName="genotypeFullActiveRegion", doc = "If specified, alternate alleles are considered to be the full active region for the purposes of genotyping", required = false)
    protected boolean GENOTYPE_FULL_ACTIVE_REGION = false;

    @Argument(fullName="fullHaplotype", shortName="fullHaplotype", doc = "If specified, output the full haplotype sequence instead of converting to individual variants w.r.t. the reference", required = false)
    protected boolean OUTPUT_FULL_HAPLOTYPE_SEQUENCE = false;

    @Argument(fullName="gcpHMM", shortName="gcpHMM", doc="gcpHMM", required = false)
    protected int gcpHMM = 10;

    @Argument(fullName="downsampleRegion", shortName="dr", doc="coverage per sample to downsample each region to", required = false)
    protected int DOWNSAMPLE_PER_SAMPLE_PER_REGION = 1000;

    @Argument(fullName="useExpandedTriggerSet", shortName="expandedTriggers", doc = "If specified, use additional, experimental triggers designed to capture larger indels but which may lead to an increase in the false positive rate", required=false)
    protected boolean USE_EXPANDED_TRIGGER_SET = false;

    @Argument(fullName="useAllelesTrigger", shortName="allelesTrigger", doc = "If specified, use additional, trigger on variants found in an external alleles file", required=false)
    protected boolean USE_ALLELES_TRIGGER = false;

    @ArgumentCollection
    private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // the calculation arguments
    private UnifiedGenotyperEngine UG_engine = null;
    private UnifiedGenotyperEngine UG_engine_simple_genotyper = null;
    
    @Argument(fullName="debug", shortName="debug", doc="If specified, print out very verbose debug information about each triggering active region", required = false)
    protected boolean DEBUG;

    @Argument(fullName="doBanded", shortName="doBanded", doc="If specified, use the banded option", required = false)
    protected boolean doBanded = false;

    // the assembly engine
    LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    LikelihoodCalculationEngine likelihoodCalculationEngine = null;

    // the genotyping engine
    GenotypingEngine genotypingEngine = null;

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    // fasta reference reader to supplement the edges of the reference sequence
    private IndexedFastaSequenceFile referenceReader;

    // reference base padding size
    private static final int REFERENCE_PADDING = 120;

    // bases with quality less than or equal to this value are trimmed off the tails of the reads
    private static final byte MIN_TAIL_QUALITY = 18;

    private ArrayList<String> samplesList = new ArrayList<String>();
    private final static double LOG_ONE_HALF = -Math.log10(2.0);
    private final static double LOG_ONE_THIRD = -Math.log10(3.0);
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
        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC.clone(), logger, null, null, samples, VariantContextUtils.DEFAULT_PLOIDY);
        UAC.OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY; // low values used for isActive determination only, default/user-specified values used for actual calling
        UAC.GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY; // low values used for isActive determination only, default/user-specified values used for actual calling
        UAC.STANDARD_CONFIDENCE_FOR_CALLING = (USE_EXPANDED_TRIGGER_SET ? 0.3 : 4.0); // low values used for isActive determination only, default/user-specified values used for actual calling
        UAC.STANDARD_CONFIDENCE_FOR_EMITTING = (USE_EXPANDED_TRIGGER_SET ? 0.3 : 4.0); // low values used for isActive determination only, default/user-specified values used for actual calling
        UG_engine_simple_genotyper = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, samples, VariantContextUtils.DEFAULT_PLOIDY);

        // initialize the output VCF header
        vcfWriter.writeHeader(new VCFHeader(new HashSet<VCFHeaderLine>(), samples));

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile, e);
        }

        assemblyEngine = new SimpleDeBruijnAssembler( DEBUG, graphWriter );
        likelihoodCalculationEngine = new LikelihoodCalculationEngine( (byte)gcpHMM, DEBUG, doBanded );
        genotypingEngine = new GenotypingEngine( DEBUG, MNP_LOOK_AHEAD, OUTPUT_FULL_HAPLOTYPE_SEQUENCE );
        annotationEngine = new VariantAnnotatorEngine(getToolkit());
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
    @Ensures({"result >= 0.0", "result <= 1.0"})
    public double isActive( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {

        /*
        if( UG_engine.getUAC().GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            for( final VariantContext vc : tracker.getValues(UG_engine.getUAC().alleles) ) {
                if(!allelesToGenotype.contains(vc)) {
                    allelesToGenotype.add(vc);
                }
            }
            if( tracker.getValues(UG_engine.getUAC().alleles).size() > 0 ) {
                return 1.0;
            }
        }
        */
        if( USE_ALLELES_TRIGGER ) {
            return ( tracker.getValues(UG_engine.getUAC().alleles).size() > 0 ? 1.0 : 0.0 );
        }

        if( context == null ) { return 0.0; }

        final List<Allele> noCall = new ArrayList<Allele>(); // used to noCall all genotypes until the exact model is applied
        noCall.add(Allele.NO_CALL);

        final Map<String, AlignmentContext> splitContexts = AlignmentContextUtils.splitContextBySampleName(context);
        final GenotypesContext genotypes = GenotypesContext.create(splitContexts.keySet().size());
        for( final String sample : splitContexts.keySet() ) {
            final double[] genotypeLikelihoods = new double[3]; // ref versus non-ref (any event)
            Arrays.fill(genotypeLikelihoods, 0.0);

            for( final PileupElement p : splitContexts.get(sample).getBasePileup() ) {
                final byte qual = ( USE_EXPANDED_TRIGGER_SET ?
                                        ( p.isNextToSoftClip() || p.isBeforeInsertion() || p.isAfterInsertion() ? ( p.getQual() > QualityUtils.MIN_USABLE_Q_SCORE ? p.getQual() : (byte) 20 ) : p.getQual() )
                                        : p.getQual() );
                if( p.isDeletion() || qual > (USE_EXPANDED_TRIGGER_SET ? QualityUtils.MIN_USABLE_Q_SCORE : (byte) 18) ) {
                    int AA = 0; final int AB = 1; int BB = 2;
                    if( USE_EXPANDED_TRIGGER_SET ) {
                        if( p.getBase() != ref.getBase() || p.isDeletion() || p.isBeforeDeletedBase() || p.isAfterDeletedBase() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip() ||
                                (!p.getRead().getNGSPlatform().equals(NGSPlatform.SOLID) && ((p.getRead().getReadPairedFlag() && p.getRead().getMateUnmappedFlag()) || BadMateFilter.hasBadMate(p.getRead()))) ) {
                            AA = 2;
                            BB = 0;
                        }
                    } else {
                        if( p.getBase() != ref.getBase() || p.isDeletion() || p.isBeforeDeletedBase() || p.isAfterDeletedBase() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip() ) {
                            AA = 2;
                            BB = 0;
                        }
                    }
                    genotypeLikelihoods[AA] += QualityUtils.qualToProbLog10(qual);
                    genotypeLikelihoods[AB] += MathUtils.approximateLog10SumLog10( QualityUtils.qualToProbLog10(qual) + LOG_ONE_HALF, QualityUtils.qualToErrorProbLog10(qual) + LOG_ONE_THIRD + LOG_ONE_HALF );
                    genotypeLikelihoods[BB] += QualityUtils.qualToErrorProbLog10(qual) + LOG_ONE_THIRD;
                }
            }

            final HashMap<String, Object> attributes = new HashMap<String, Object>();
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods(genotypeLikelihoods));
            genotypes.add(new Genotype(sample, noCall, Genotype.NO_LOG10_PERROR, null, attributes, false));
        }

        final ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add( Allele.create("A", true) ); // fake ref Allele
        alleles.add( Allele.create("G", false) ); // fake alt Allele
        final VariantCallContext vcOut = UG_engine_simple_genotyper.calculateGenotypes(new VariantContextBuilder("HCisActive!", context.getContig(), context.getLocation().getStart(), context.getLocation().getStop(), alleles).genotypes(genotypes).make(), GenotypeLikelihoodsCalculationModel.Model.INDEL);
        return ( vcOut == null ? 0.0 : QualityUtils.qualToProb( vcOut.getPhredScaledQual() ) );
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Integer map( final ActiveRegion activeRegion, final RefMetaDataTracker metaDataTracker ) {

        final ArrayList<VariantContext> activeAllelesToGenotype = new ArrayList<VariantContext>();
        /*
        if( UG_engine.getUAC().GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            for( final VariantContext vc : allelesToGenotype ) {
                if( activeRegion.getLocation().overlapsP( getToolkit().getGenomeLocParser().createGenomeLoc(vc) ) ) {
                    activeAllelesToGenotype.add(vc); // do something with these VCs during GGA mode
                }
            }
            allelesToGenotype.removeAll( activeAllelesToGenotype );
        }
        */

        if( !activeRegion.isActive ) { return 0; } // Not active so nothing to do!
        if( activeRegion.size() == 0 && UG_engine.getUAC().GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) { return 0; } // No reads here so nothing to do!
        if( UG_engine.getUAC().GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES && activeAllelesToGenotype.isEmpty() ) { return 0; } // No alleles found in this region so nothing to do!

        finalizeActiveRegion( activeRegion ); // merge overlapping fragments, clip adapter and low qual tails
        final Haplotype referenceHaplotype = new Haplotype(activeRegion.getActiveRegionReference(referenceReader)); // Create the reference haplotype which is the bases from the reference that make up the active region
        referenceHaplotype.setIsReference(true);
        final byte[] fullReferenceWithPadding = activeRegion.getFullReference(referenceReader, REFERENCE_PADDING);
        int PRUNE_FACTOR = determinePruneFactorFromCoverage( activeRegion );
        final ArrayList<Haplotype> haplotypes = assemblyEngine.runLocalAssembly( activeRegion.getReads(), referenceHaplotype, fullReferenceWithPadding, PRUNE_FACTOR );
        if( haplotypes.size() == 1 ) { return 1; } // only the reference haplotype remains so nothing else to do!

        activeRegion.hardClipToActiveRegion(); // only evaluate the parts of reads that are overlapping the active region
        final List<GATKSAMRecord> filteredReads = filterNonPassingReads( activeRegion ); // filter out reads from genotyping which fail mapping quality based criteria
        if( activeRegion.size() == 0 ) { return 1; } // no reads remain after filtering so nothing else to do!

        // evaluate each sample's reads against all haplotypes
        final HashMap<String, ArrayList<GATKSAMRecord>> perSampleReadList = splitReadsBySample( activeRegion.getReads() );
        final HashMap<String, ArrayList<GATKSAMRecord>> perSampleFilteredReadList = splitReadsBySample( filteredReads );
        likelihoodCalculationEngine.computeReadLikelihoods( haplotypes, perSampleReadList );

        // subset down to only the best haplotypes to be genotyped in all samples
        final ArrayList<Haplotype> bestHaplotypes = likelihoodCalculationEngine.selectBestHaplotypes( haplotypes );

        for( final Pair<VariantContext, ArrayList<ArrayList<Haplotype>>> callResult :
                ( GENOTYPE_FULL_ACTIVE_REGION
                  ? genotypingEngine.assignGenotypeLikelihoodsAndCallHaplotypeEvents( UG_engine, bestHaplotypes, fullReferenceWithPadding, getPaddedLoc(activeRegion), activeRegion.getLocation(), getToolkit().getGenomeLocParser() )
                  : genotypingEngine.assignGenotypeLikelihoodsAndCallIndependentEvents( UG_engine, bestHaplotypes, fullReferenceWithPadding, getPaddedLoc(activeRegion), activeRegion.getLocation(), getToolkit().getGenomeLocParser() ) ) ) {
            if( DEBUG && samplesList.size() <= 10 ) { System.out.println(callResult.getFirst()); }

            final Map<String, Map<Allele, List<GATKSAMRecord>>> stratifiedReadMap = LikelihoodCalculationEngine.partitionReadsBasedOnLikelihoods(perSampleReadList, perSampleFilteredReadList, callResult);
            final VariantContext annotatedCall = annotationEngine.annotateContext(stratifiedReadMap, callResult.getFirst());

            vcfWriter.add(annotatedCall);
        }

        if( DEBUG ) { System.out.println("----------------------------------------------------------------------------------"); }

        return 1; // One active region was processed during this map call
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

    private void finalizeActiveRegion( final ActiveRegion activeRegion ) {
        if( DEBUG ) { System.out.println("\nAssembling " + activeRegion.getLocation() + " with " + activeRegion.size() + " reads:"); }
        final ArrayList<GATKSAMRecord> finalizedReadList = new ArrayList<GATKSAMRecord>();
        final FragmentCollection<GATKSAMRecord> fragmentCollection = FragmentUtils.create( ReadUtils.sortReadsByCoordinate(activeRegion.getReads()) );
        activeRegion.clearReads();

        // Join overlapping paired reads to create a single longer read
        finalizedReadList.addAll( fragmentCollection.getSingletonReads() );
        for( final List<GATKSAMRecord> overlappingPair : fragmentCollection.getOverlappingPairs() ) {
            finalizedReadList.addAll( FragmentUtils.mergeOverlappingPairedFragments(overlappingPair) );
        }

        Collections.shuffle(finalizedReadList, GenomeAnalysisEngine.getRandomGenerator());

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

    private List<GATKSAMRecord> filterNonPassingReads( final ActiveRegion activeRegion ) {
        final ArrayList<GATKSAMRecord> readsToRemove = new ArrayList<GATKSAMRecord>();
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {
            if( rec.getReadLength() < 24 || rec.getMappingQuality() <= 20 || BadMateFilter.hasBadMate(rec) || (keepRG != null && !rec.getReadGroup().getId().equals(keepRG)) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll( readsToRemove );
        return readsToRemove;
    }

    private GenomeLoc getPaddedLoc( final ActiveRegion activeRegion ) {
        final int padLeft = Math.max(activeRegion.getReferenceLoc().getStart()-REFERENCE_PADDING, 1);
        final int padRight = Math.min(activeRegion.getReferenceLoc().getStop()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(activeRegion.getReferenceLoc().getContig()).getSequenceLength());
        return getToolkit().getGenomeLocParser().createGenomeLoc(activeRegion.getReferenceLoc().getContig(), padLeft, padRight);
    }

    private HashMap<String, ArrayList<GATKSAMRecord>> splitReadsBySample( final List<GATKSAMRecord> reads ) {
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
    
    private int determinePruneFactorFromCoverage( final ActiveRegion activeRegion ) {
        final ArrayList<Integer> readLengthDistribution = new ArrayList<Integer>();
        for( final GATKSAMRecord read : activeRegion.getReads() ) {
            readLengthDistribution.add(read.getReadLength());
        }
        final double meanReadLength = MathUtils.average(readLengthDistribution);
        final double meanCoveragePerSample = (double) activeRegion.getReads().size() / ((double) activeRegion.getExtendedLoc().size() / meanReadLength) / (double) samplesList.size();
        int PRUNE_FACTOR = 0;
        if( meanCoveragePerSample > 6.0 ) {
            PRUNE_FACTOR = (int) Math.floor( Math.sqrt( meanCoveragePerSample - 2.0 ) );
        } else if( meanCoveragePerSample > 2.5 ) {
            PRUNE_FACTOR = 1;
        }

        if( DEBUG ) { System.out.println(String.format("Mean coverage per sample = %.1f --> prune factor = %d", meanCoveragePerSample, PRUNE_FACTOR)); }
        return PRUNE_FACTOR;
    }
}