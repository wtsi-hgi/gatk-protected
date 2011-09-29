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
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.indels.ConstrainedMateFixingManager;
import org.broadinstitute.sting.gatk.walkers.recalibration.Covariate;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDataManager;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum;
import org.broadinstitute.sting.gatk.walkers.recalibration.TableRecalibrationWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.codecs.samread.SAMReadCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
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
public class HaplotypeCaller extends ReadWalker<SAMRecord, Integer> implements TreeReducible<Integer> {

    /**
     * A raw, unfiltered, highly specific callset in VCF format.
     */
    @Output(doc="File to which variants should be written", required = true)
    protected VCFWriter vcfWriter = null;

    @Output(fullName="graphOutput", shortName="graph", doc="File to which debug assembly graph information should be written", required = false)
    protected PrintStream graphWriter = null;

    @Output(fullName="debugBam", shortName="debugBam", doc="File to which all possible haplotypes in bam format (aligned via SW) should be written", required = false)
    protected StingSAMFileWriter bamWriter = null;

    @Argument(fullName="realignReads", shortName="realignReads", doc="If provided, the debugBam will contain all reads in the interval realigned to the new haplotype", required = false)
    protected boolean realignReads = false;

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

    // should we print out verbose debug information about each triggering interval
    private final boolean DEBUG = true;

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
    private static final int REFERENCE_PADDING = 110;

    // bases with quality less than or equal to this value are trimmed off the tails of the reads
    private static final byte MIN_TAIL_QUALITY = 6;

    protected ConstrainedMateFixingManager manager = null;

    private NestedHashMap kmerQualityTables = new NestedHashMap();
    private int contextSize = 0;

    public void initialize() {

        // read in the input recal data from IndelCountCovariates and calculate the empirical gap open penalty for each k-mer
        parseInputRecalData();

        // get all of the unique sample names
        Set<String> samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        // initialize the header
        vcfWriter.writeHeader(new VCFHeader(new HashSet<VCFHeaderLine>(), samples));

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile,ex);
        }

        assemblyEngine = makeAssembler(ASSEMBLER_TO_USE, referenceReader);
        likelihoodCalculationEngine = new LikelihoodCalculationEngine(gopHMM, gcpHMM, false, true, false, kmerQualityTables, contextSize);
        genotypingEngine = new GenotypingEngine( DEBUG, gopSW, gcpSW );

        if( realignReads ) {
            manager = new ConstrainedMateFixingManager(bamWriter, getToolkit().getGenomeLocParser(), 3000, 200, 150000);
        }

        GenomeLocSortedSet intervalsToAssemble = getToolkit().getIntervals();
        if ( intervalsToAssemble == null || intervalsToAssemble.isEmpty() )
            throw new UserException.BadInput("Intervals must be provided with -L or -BTI (preferably not larger than several hundred bp)");

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

    public SAMRecord map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        return currentInterval == null ? null : read;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(SAMRecord read, Integer sum) {
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
            readsToAssemble.add(read); // don't want this triggering read which is past the interval to fall through the cracks?
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
        if ( readsToAssemble.size() > 0 ) {
            processReadBin( currentInterval );
            result++;
        }
        logger.info("Ran local assembly on " + result + " intervals");

        if( realignReads ) {
            manager.close();
        }
    }

    private void processReadBin( final GenomeLoc curInterval ) {

        if( DEBUG ) { System.out.println(curInterval.getLocation() + " with " + readsToAssemble.getReads().size() + " reads:"); }

        final List<Haplotype> haplotypes = assemblyEngine.runLocalAssembly( readsToAssemble.getReads() );

        // Add the full reference as a possible haplotype
        final Haplotype refHaplotype = new Haplotype(readsToAssemble.getReferenceNoPadding(referenceReader));
        haplotypes.add(refHaplotype);

        if( DEBUG ) { System.out.println("Found " + haplotypes.size() + " candidate haplotypes to evaluate"); }

        if( haplotypes.size() == 0 ) {
            if( DEBUG ) { System.out.println("WARNING! No haplotypes created during assembly!"); }
            return;
        }

        if( bamWriter != null && !realignReads ) {
            genotypingEngine.alignAllHaplotypes( haplotypes, readsToAssemble.getReference( referenceReader ), readsToAssemble.getLocation(), bamWriter, readsToAssemble.getReads().get(0) );
            return; // in assembly debug mode, so no need to run the rest of the procedure
        }

        final int pos = curInterval.getStart() + (curInterval.getStop() - curInterval.getStart()) / 2;
        final GenomeLoc evalWindow = getToolkit().getGenomeLocParser().createGenomeLoc(curInterval.getContig(), pos - 8, pos + 8);
        final GenomeLoc outputWindow = getToolkit().getGenomeLocParser().createGenomeLoc(curInterval.getContig(), pos - 1, pos + 1);

        final Pair<Haplotype, Haplotype> bestTwoHaplotypes = likelihoodCalculationEngine.computeLikelihoods( haplotypes, readsToAssemble.getReadsInWindow( evalWindow ) );
        final List<VariantContext> vcs = genotypingEngine.alignAndGenotype( bestTwoHaplotypes, readsToAssemble.getReference( referenceReader ), readsToAssemble.getLocation(), bestTwoHaplotypes.first.likelihood );

        if( bamWriter != null && realignReads ) {
            genotypingEngine.alignAllReads( bestTwoHaplotypes, readsToAssemble.getReference( referenceReader ), readsToAssemble.getLocation(), manager, readsToAssemble.getReadsInWindow( evalWindow ), likelihoodCalculationEngine.readLikelihoodsForBestHaplotypes );
        }

        for( final VariantContext vc : vcs ) {
            if( vc.getStart() >= outputWindow.getStart() && vc.getStart() <= outputWindow.getStop() ) {
                if( DEBUG ) { System.out.print("== "); }
                vcfWriter.add(vc);
            }

            if( DEBUG ) { System.out.println(vc); }
        }

        if( DEBUG ) { System.out.println("----------------------------------------------------------------------------------"); }

        readsToAssemble.clear();
    }

    // private class copied from IndelRealigner, used to bin together a bunch of reads and then retrieve the reference overlapping the full extent of the bin
    // the precursor to the Active Region Traversal
    private class ReadBin implements HasGenomeLocation {

        private final ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        private byte[] reference = null;
        private GenomeLoc loc = null;

        public ReadBin() { }

        // Return false if we can't process this read bin because the reads are not correctly overlapping.
        // This can happen if e.g. there's a large known indel with no overlapping reads.
        public void add( final SAMRecord read ) {

            if( reads.size() < 300 ) { // protection against pileups with abnormally deep coverage
                final GATKSAMRecord postAdapterRead = ReadUtils.hardClipAdaptorSequence(read);
                if( postAdapterRead != null ) {
                    final SAMRecord clippedRead = (new ReadClipper(postAdapterRead)).hardClipLowQualEnds( MIN_TAIL_QUALITY );

                    if( clippedRead.getReadLength() > 0 ) {
                        final GenomeLoc locForRead = getToolkit().getGenomeLocParser().createGenomeLoc(clippedRead);
                        if ( loc == null )
                            loc = locForRead;
                        else if ( locForRead.getStop() > loc.getStop() )
                            loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), loc.getStart(), locForRead.getStop());

                        reads.add(clippedRead);
                    }
                }
            }
        }

        public List<SAMRecord> getReads() { return reads; }

        public List<SAMRecord> getReadsInWindow( final GenomeLoc window ) {
            final ArrayList<SAMRecord> readsOverlappingVariant = new ArrayList<SAMRecord>();

            for( final SAMRecord rec : reads ) {
                if( rec.getMappingQuality() > 15 && !BadMateFilter.hasBadMate(rec) ) {
                    final GenomeLoc locForRead = getToolkit().getGenomeLocParser().createGenomeLoc(rec);
                    if( locForRead.overlapsP(window) ) {
                        readsOverlappingVariant.add(rec);
                    }
                }
            }

            return readsOverlappingVariant;
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
}