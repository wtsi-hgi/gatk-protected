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
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.gatk.walkers.recalibration.TableRecalibrationWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.fragments.FragmentCollection;
import org.broadinstitute.sting.utils.fragments.FragmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
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
@PartitionBy(PartitionType.INTERVAL)
@ReadFilters( {MappingQualityUnavailableFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class, FailsVendorQualityCheckFilter.class} )
public class HaplotypeCaller extends ReadWalker<GATKSAMRecord, Integer> implements TreeReducible<Integer> {

    /**
     * A raw, unfiltered, highly specific callset in VCF format.
     */
    @Output(doc="File to which variants should be written", required = true)
    protected VCFWriter vcfWriter = null;

    @Output(fullName="graphOutput", shortName="graph", doc="File to which debug assembly graph information should be written", required = false)
    protected PrintStream graphWriter = null;

    /**
     * After the header, data records occur one per line until the end of the file. The first several items on a line are the
     * values of the individual covariates and will change depending on which covariates were specified at runtime. The last
     * three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches,
     * and the raw empirical quality score calculated by phred-scaling the mismatch rate.
     */
    @Input(fullName="recal_file", shortName="recalFile", required=true, doc="Filename for the input indel covariates table recalibration .csv file")
    public File RECAL_FILE = null;

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

    @Argument(fullName="debug", shortName="debug", doc="If specified print out very verbose debug information about each triggering interval", required = false)
    protected boolean DEBUG;

    // the assembly engine
    LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    LikelihoodCalculationEngine likelihoodCalculationEngine = null;

    // the genotyping engine
    GenotypingEngine genotypingEngine = null;

    // the intervals input by the user
    private Iterator<GenomeLoc> intervals = null;

    // the current interval in the list
    private GenomeLoc currentInterval = null;

    // the reads that fall into the current interval
    private final ReadBin readsToAssemble = new ReadBin();

    // fasta reference reader to supplement the edges of the reference sequence
    private IndexedFastaSequenceFile referenceReader;

    // reference base padding size
    private static final int REFERENCE_PADDING = 120;

    // bases with quality less than or equal to this value are trimmed off the tails of the reads
    private static final byte MIN_TAIL_QUALITY = 8;

    private NestedHashMap kmerQualityTables = new NestedHashMap();
    private int contextSize = 0; // to be set when reading in the recal k-mer tables
    private ArrayList<String> samplesList = new ArrayList<String>();

    public void initialize() {

        // read in the input recal data from IndelCountCovariates and calculate the empirical gap open penalty for each k-mer
        parseInputRecalData();

        // get all of the unique sample names
        Set<String> samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        samplesList.addAll( samples );
        // initialize the UnifiedGenotyper Engine which is used to call into the exact model
        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, samples);
        // initialize the header
        vcfWriter.writeHeader(new VCFHeader(new HashSet<VCFHeaderLine>(), samples));

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile, e);
        }

        assemblyEngine = makeAssembler(ASSEMBLER_TO_USE, referenceReader);
        likelihoodCalculationEngine = new LikelihoodCalculationEngine(gopHMM, gcpHMM, DEBUG, true, false, kmerQualityTables, contextSize);
        genotypingEngine = new GenotypingEngine( DEBUG, gopSW, gcpSW );

        GenomeLocSortedSet intervalsToAssemble = getToolkit().getIntervals();
        if ( intervalsToAssemble == null || intervalsToAssemble.isEmpty() )
            throw new UserException.BadInput("Intervals must be provided with -L (preferably not larger than several hundred bp)");

        intervals = intervalsToAssemble.clone().iterator();
        currentInterval = intervals.hasNext() ? intervals.next() : null;
    }

    private void parseInputRecalData() {
        int lineNumber = 0;
        boolean foundAllCovariates = false;

        // Read in the data from the csv file and populate the data map and covariates list
        logger.info( "Reading in the data from input csv file..." );

        boolean sawEOF = false;
        try {
            for ( String line : new XReadLines(RECAL_FILE) ) {
                lineNumber++;
                if (TableRecalibrationWalker.EOF_MARKER.equals(line) ) {
                    sawEOF = true;
                } else if( TableRecalibrationWalker.COMMENT_PATTERN.matcher(line).matches() || TableRecalibrationWalker.OLD_RECALIBRATOR_HEADER.matcher(line).matches() ||
                        TableRecalibrationWalker.COVARIATE_PATTERN.matcher(line).matches() )  {
                    // Skip over the comment lines, (which start with '#')
                } else { // Found a line of data
                    final String[] vals = line.split(",");
                    final Object[] key = new Object[2];
                    key[0] = vals[0]; // read group
                    key[1] = vals[2]; // k-mer Context
                    final Double logGapOpenPenalty = Math.log10( (Double.parseDouble(vals[4]) + 1.0) / (Double.parseDouble(vals[3]) + 1.0) );
                    kmerQualityTables.put( logGapOpenPenalty, key );
                    if(contextSize == 0) { contextSize = key[1].toString().length(); }
                }
            }

        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(RECAL_FILE, "Can not find input file", e);
        } catch ( NumberFormatException e ) {
            throw new UserException.MalformedFile(RECAL_FILE, "Error parsing recalibration data at line " + lineNumber + ". Perhaps your table was generated by an older version of CovariateCounterWalker.");
        }
        logger.info( "...done!" );

        if ( !sawEOF ) {
            final String errorMessage = "No EOF marker was present in the recal covariates table; this could mean that the file is corrupted or was generated with an old version of the CountCovariates tool.";
            throw new UserException.MalformedFile(RECAL_FILE, errorMessage);
        }
    }

    private LocalAssemblyEngine makeAssembler(LocalAssemblyEngine.ASSEMBLER type, IndexedFastaSequenceFile referenceReader) {
        switch ( type ) {
            case SIMPLE_DE_BRUIJN:
                return new SimpleDeBruijnAssembler(graphWriter, referenceReader);
            default:
                throw new UserException.BadInput("Assembler type " + type + " is not valid/supported");
        }
    }

    public GATKSAMRecord map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {
        return currentInterval == null ? null : read;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(GATKSAMRecord read, Integer sum) {
        if ( read == null )
            return sum;

        GenomeLoc readLoc = getToolkit().getGenomeLocParser().createGenomeLoc(read);
        // hack to get around unmapped reads having screwy locations
        if ( readLoc.getStop() == 0 )
            readLoc = getToolkit().getGenomeLocParser().createGenomeLoc(readLoc.getContig(), readLoc.getStart(), readLoc.getStart());

        if ( readLoc.overlapsP(currentInterval) ) {
            readsToAssemble.add(read);
        } else {
            processReadBin( currentInterval );
            readsToAssemble.add(read); // don't want this triggering read which is past the interval to fall through the cracks
            sum++;

            do {
                currentInterval = intervals.hasNext() ? intervals.next() : null;
            } while ( currentInterval != null && currentInterval.isBefore(readLoc) );
        }

        return sum;
    }

    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone(Integer result) {
        processReadBin( currentInterval );
        result++;
        logger.info("Ran local assembly on " + result + " intervals");
    }

    private void processReadBin( final GenomeLoc curInterval ) {
        if ( readsToAssemble.size() == 0 ) { return; } // No reads here so nothing to do!

        if( DEBUG ) { System.out.println("Assembling " + curInterval.getLocation() + " with " + readsToAssemble.getReads().size() + " reads:"); }
        readsToAssemble.finalizeBin( curInterval );
        final ArrayList<Haplotype> haplotypes = assemblyEngine.runLocalAssembly( readsToAssemble.getReads(), new Haplotype(readsToAssemble.getReferenceNoPadding(referenceReader) ) );
        if( DEBUG ) { System.out.println("Found " + haplotypes.size() + " candidate haplotypes to evaluate"); }
        genotypingEngine.createEventDictionaryAndFilterBadHaplotypes( haplotypes, readsToAssemble.getReference(referenceReader), readsToAssemble.getLocation(), curInterval );
        if( DEBUG ) { System.out.println(haplotypes.size() + " candidate haplotypes remain after filtering"); }

        if( haplotypes.size() == 0 ) {
            if( DEBUG ) { System.out.println("WARNING! No haplotypes created during assembly!"); }
            return;
        }

        final HashMap<String, Double[][]> haplotypeLikehoodMatrixMap = new HashMap<String, Double[][]>();
        final ArrayList<Haplotype> bestTwoHaplotypesPerSample = new ArrayList<Haplotype>();
        final HashMap<String, ArrayList<GATKSAMRecord>> readListMap = splitReadsBySample( readsToAssemble.getPassingReads() );

        for( final String sample : readListMap.keySet() ) {
            if( DEBUG ) { System.out.println("Evaluating sample " + sample + " with " + readListMap.get( sample ).size() + " passing reads"); }
            likelihoodCalculationEngine.computeLikelihoods( haplotypes, readListMap.get( sample ) );
            bestTwoHaplotypesPerSample.addAll( likelihoodCalculationEngine.chooseBestHaplotypes(haplotypes) );
            haplotypeLikehoodMatrixMap.put( sample, likelihoodCalculationEngine.haplotypeLikehoodMatrix );
        }

        final ArrayList<VariantContext> vcs = genotypingEngine.alignAndAssignGenotypeLikelihoods( getToolkit().getGenomeLocParser(), haplotypes, bestTwoHaplotypesPerSample, readsToAssemble.getReference(referenceReader), readsToAssemble.getLocation(), curInterval, haplotypeLikehoodMatrixMap );

        for( final VariantContext vc : vcs ) {
            if( curInterval.containsP(getToolkit().getGenomeLocParser().createGenomeLoc(vc).getStartLocation()) ) {
                final VariantCallContext vcOut = UG_engine.calculateGenotypes(vc, UG_engine.getUAC().GLmodel);
                if(vcOut != null) {
                    if( DEBUG ) { System.out.println(vcOut); }
                    vcfWriter.add(vcOut);
                }
            }
        }

        if( DEBUG ) { System.out.println("----------------------------------------------------------------------------------"); }

        readsToAssemble.clear();
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

    // private class copied from IndelRealigner, used to bin together a bunch of reads and then retrieve the reference overlapping the full extent of the bin
    // the precursor to the Active Region Traversal
    private class ReadBin implements HasGenomeLocation {

        private final ArrayList<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
        private byte[] reference = null;
        private GenomeLoc loc = null;

        public ReadBin() { }

        // add each read to the bin and extend the reference genome loc if needed
        public void add( final GATKSAMRecord read ) {
            final GenomeLoc locForRead = getToolkit().getGenomeLocParser().createGenomeLoc( read );
            if ( loc == null )
                loc = locForRead;
            else if ( locForRead.getStop() > loc.getStop() )
                loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), loc.getStart(), locForRead.getStop());

            reads.add( read );
        }

        public void finalizeBin( final GenomeLoc curInterval ) {
            loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(),
                    Math.min(loc.getStart(), curInterval.getStart()), Math.max(loc.getStop(), curInterval.getStop()));
            final ArrayList<GATKSAMRecord> finalizedReadList = new ArrayList<GATKSAMRecord>();
            final FragmentCollection<GATKSAMRecord> fragmentCollection = FragmentUtils.create( reads );
            reads.clear();

            // Join overlapping paired reads to create a single longer read
            finalizedReadList.addAll( fragmentCollection.getSingletonReads() );
            for( final List<GATKSAMRecord> overlappingPair : fragmentCollection.getOverlappingPairs() ) {
                finalizedReadList.addAll(mergeOverlappingPairedReads(overlappingPair));
            }

            Collections.shuffle(finalizedReadList);

            // Loop through the reads hard clipping the adaptor and low quality tails
            for( final GATKSAMRecord myRead : finalizedReadList ) {
                final GATKSAMRecord postAdapterRead = ReadUtils.hardClipAdaptorSequence( myRead );
                if( postAdapterRead != null && postAdapterRead.getCigar().getReadLength() > 0 ) {
                    final GATKSAMRecord clippedRead = (new ReadClipper(postAdapterRead)).hardClipLowQualEnds( MIN_TAIL_QUALITY );

                    // protect against INTERVALS with abnormally high coverage
                    if( clippedRead.getReadLength() > 0 && reads.size() < samplesList.size() * DOWNSAMPLE_PER_SAMPLE_PER_REGION ) {
                        reads.add(clippedRead);
                    }
                }
            }
        }

        public ArrayList<GATKSAMRecord> getReads() { return reads; }

        public ArrayList<GATKSAMRecord> getPassingReads() {
            final ArrayList<GATKSAMRecord> passingReads = new ArrayList<GATKSAMRecord>();

            for( final GATKSAMRecord rec : reads ) {
                if( rec.getMappingQuality() > 18 && !BadMateFilter.hasBadMate(rec) ) {
                    passingReads.add(rec);
                }
            }

            return passingReads;
        }

        public byte[] getReference(IndexedFastaSequenceFile referenceReader) {
            // set up the reference if we haven't done so yet
            if ( reference == null ) {
                // first, pad the reference to handle deletions in narrow windows (e.g. those with only 1 read)
                int padLeft = Math.max(loc.getStart()-REFERENCE_PADDING, 1);
                int padRight = Math.min(loc.getStop()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(loc.getContig()).getSequenceLength());
                loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), padLeft, padRight);
                reference = referenceReader.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
            }

            return reference;
        }

        public byte[] getReferenceNoPadding(IndexedFastaSequenceFile referenceReader) {
            int padLeft = Math.max(loc.getStart(), 1);
            int padRight = Math.min(loc.getStop(), referenceReader.getSequenceDictionary().getSequence(loc.getContig()).getSequenceLength());
            final GenomeLoc thisLoc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), padLeft, padRight);
            return referenceReader.getSubsequenceAt(thisLoc.getContig(), thisLoc.getStart(), thisLoc.getStop()).getBases();
        }

        public GenomeLoc getLocation() { return loc; }

        public int size() { return reads.size(); }

        public void clear() {
            reads.clear();
            reference = null;
            loc = null;
        }
    }

    // BUGBUG: move this function to ReadUtils or FragmentUtils
    private List<GATKSAMRecord> mergeOverlappingPairedReads( List<GATKSAMRecord> overlappingPair ) {
        if( overlappingPair.size() != 2 ) { throw new ReviewedStingException("Found overlapping pair with " + overlappingPair.size() + " reads, but expecting exactly 2."); }

        GATKSAMRecord firstRead = overlappingPair.get(0);
        GATKSAMRecord secondRead = overlappingPair.get(1);
        if( !(secondRead.getUnclippedStart() <= firstRead.getUnclippedEnd() && secondRead.getUnclippedStart() >= firstRead.getUnclippedStart() && secondRead.getUnclippedEnd() >= firstRead.getUnclippedEnd()) ) {
            firstRead = overlappingPair.get(1);
            secondRead = overlappingPair.get(0);
        }
        if( !(secondRead.getUnclippedStart() <= firstRead.getUnclippedEnd() && secondRead.getUnclippedStart() >= firstRead.getUnclippedStart() && secondRead.getUnclippedEnd() >= firstRead.getUnclippedEnd()) ) {
            return overlappingPair; // can't merge them, yet:  AAAAAAAAAAA-BBBBBBBBBBB-AAAAAAAAAAAAAA
        }

        // DEBUG PRINTING
        //System.out.println("1: " + firstRead.getUnclippedStart() + " - " + firstRead.getUnclippedEnd() + "  > " + firstRead.getCigar() + " " + firstRead.getReadGroup().getId());
        //System.out.println("2: " + secondRead.getUnclippedStart() + " - " + secondRead.getUnclippedEnd() + "  > " + secondRead.getCigar()+ " " + firstRead.getReadGroup().getId());
        //System.out.println(firstRead.getReadString());
        //System.out.println(secondRead.getReadString());

        final int firstReadStop = ReadUtils.getReadCoordinateForReferenceCoordinate(firstRead, secondRead.getUnclippedStart(), ReadUtils.ClippingTail.RIGHT_TAIL);
        final int numBases = firstReadStop + secondRead.getReadLength();
        final byte[] bases = new byte[numBases];
        final byte[] quals = new byte[numBases];
        final byte[] firstReadBases = firstRead.getReadBases();
        final byte[] firstReadQuals = firstRead.getBaseQualities();
        final byte[] secondReadBases = secondRead.getReadBases();
        final byte[] secondReadQuals = secondRead.getBaseQualities();
        for(int iii = 0; iii < firstReadStop; iii++) {
            bases[iii] = firstReadBases[iii];
            quals[iii] = firstReadQuals[iii];
        }
        for(int iii = firstReadStop; iii < firstRead.getReadLength(); iii++) {
            bases[iii] = ( firstReadQuals[iii] > secondReadQuals[iii-firstReadStop] ? firstReadBases[iii] : secondReadBases[iii-firstReadStop] );
            quals[iii] = ( firstReadQuals[iii] > secondReadQuals[iii-firstReadStop] ? firstReadQuals[iii] : secondReadQuals[iii-firstReadStop] );
        }
        for(int iii = firstRead.getReadLength(); iii < numBases; iii++) {
            bases[iii] = secondReadBases[iii-firstReadStop];
            quals[iii] = secondReadQuals[iii-firstReadStop];
        }

        // DEBUG PRINTING
        /*
        if( DEBUG ) {
            System.out.println(firstRead.getReadString());
            for(int jjj = 0; jjj < firstReadStop; jjj++) {
                System.out.print(" ");
            }
            System.out.println(secondRead.getReadString());
            String displayString = "";
            for(int jjj = 0; jjj < bases.length; jjj++) {
                displayString += (char) bases[jjj];
            }
            System.out.println(displayString);
        }
        */
        
        final GATKSAMRecord returnRead = new GATKSAMRecord(firstRead.getHeader());
        returnRead.setAlignmentStart(firstRead.getUnclippedStart());
        returnRead.setReadBases( bases );
        returnRead.setBaseQualities( quals );
        returnRead.setReadGroup( firstRead.getReadGroup() );
        returnRead.setReferenceName( firstRead.getReferenceName() );
        final CigarElement c = new CigarElement(bases.length, CigarOperator.M);
        final ArrayList<CigarElement> cList = new ArrayList<CigarElement>();
        cList.add(c);
        returnRead.setCigar( new Cigar( cList ));
        returnRead.setMappingQuality( firstRead.getMappingQuality() );

        final ArrayList<GATKSAMRecord> returnList = new ArrayList<GATKSAMRecord>();
        returnList.add(returnRead);
        return returnList;
    }
}