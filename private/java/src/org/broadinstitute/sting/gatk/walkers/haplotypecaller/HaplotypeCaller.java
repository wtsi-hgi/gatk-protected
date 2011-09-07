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
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.indels.ConstrainedMateFixingManager;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

@ReadFilters( {MappingQualityUnavailableFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class, FailsVendorQualityCheckFilter.class} )
public class HaplotypeCaller extends ReadWalker<SAMRecord, Integer> {

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

    @Argument(fullName = "assembler", shortName = "assembler", doc = "Assembler to use; currently only SIMPLE_DE_BRUIJN is available.", required = false)
    protected LocalAssemblyEngine.ASSEMBLER ASSEMBLER_TO_USE = LocalAssemblyEngine.ASSEMBLER.SIMPLE_DE_BRUIJN;

    // the assembly engine
    LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    LikelihoodCalculationEngine likelihoodCalculationEngine = new LikelihoodCalculationEngine(35.0, 10.0, false, true, false);

    // the genotyping engine
    GenotypingEngine genotypingEngine = new GenotypingEngine();

    // the intervals input by the user
    private Iterator<GenomeLoc> intervals = null;

    // the current interval in the list
    private GenomeLoc currentInterval = null;

    // the reads that fall into the current interval
    private final ReadBin readsToAssemble = new ReadBin();

    // fasta reference reader to supplement the edges of the reference sequence
    private IndexedFastaSequenceFile referenceReader;

    // reference base padding size
    private static final int REFERENCE_PADDING = 50;

    // bases with quality less than or equal to this value are trimmed off the tails of the reads
    private static final byte MIN_TAIL_QUALITY = 6;

    protected ConstrainedMateFixingManager manager = null;

    public void initialize() {

        // get all of the unique sample names
        // if we're supposed to assume a single sample, do so
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

        if( realignReads ) {
            manager = new ConstrainedMateFixingManager(bamWriter, getToolkit().getGenomeLocParser(), 3000, 200, 150000);
        }

        GenomeLocSortedSet intervalsToAssemble = getToolkit().getIntervals();
        if ( intervalsToAssemble == null || intervalsToAssemble.isEmpty() )
            throw new UserException.BadInput("Intervals must be provided with -L or -BTI (preferably not larger than several hundred bp)");

        intervals = intervalsToAssemble.clone().iterator();
        currentInterval = intervals.hasNext() ? intervals.next() : null;
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
            readsToAssemble.clear();
            readsToAssemble.add(read); // don't want this triggering read which is past the interval to fall through the cracks?
            sum++;

            do {
                currentInterval = intervals.hasNext() ? intervals.next() : null;
            } while ( currentInterval != null && currentInterval.isBefore(readLoc) );
        }

        return sum;
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

        System.out.println(curInterval.getLocation() + " with " + readsToAssemble.getReads().size() + " reads:");
        final List<Haplotype> haplotypes = assemblyEngine.runLocalAssembly( readsToAssemble.getReads() );
        System.out.println("Found " + haplotypes.size() + " potential haplotypes to evaluate");

        if( haplotypes.size() == 0 ) {
            System.out.println("WARNING! No haplotypes created during assembly!");
            return;
        }

        if( bamWriter != null && !realignReads ) {
            genotypingEngine.alignAllHaplotypes( haplotypes, readsToAssemble.getReference( referenceReader ), readsToAssemble.getLocation(), bamWriter, readsToAssemble.getReads().get(0) );
            return; // in assembly debug mode, so no need to run the rest of the procedure
        }

        final int pos = curInterval.getStart() + (curInterval.getStop() - curInterval.getStart()) / 2;
        final GenomeLoc window = getToolkit().getGenomeLocParser().createGenomeLoc(curInterval.getContig(), pos - 7, pos + 7);

        final Pair<Haplotype, Haplotype> bestTwoHaplotypes = likelihoodCalculationEngine.computeLikelihoods( haplotypes, readsToAssemble.getReadsInWindow( window ) );
        final List<VariantContext> vcs = genotypingEngine.alignAndGenotype( bestTwoHaplotypes, readsToAssemble.getReference( referenceReader ), readsToAssemble.getLocation() );

        if( bamWriter != null && realignReads ) {
            genotypingEngine.alignAllReads( bestTwoHaplotypes, readsToAssemble.getReference( referenceReader ), readsToAssemble.getLocation(), manager, readsToAssemble.getReadsInWindow( window ), likelihoodCalculationEngine.readLikelihoodsForBestHaplotypes );
        }

        for( final VariantContext vc : vcs ) {
            if( vc.getStart() >= window.getStart() && vc.getStart() <= window.getStop() ) {
                System.out.print("== ");
                vcfWriter.add(vc);
            }
            System.out.println(vc);
        }

        System.out.println("----------------------------------------------------------------------------------");

    }

    // private class copied from IndelRealigner, used to bin together a bunch of reads and then retrieve the reference overlapping the full extent of the bin
    private class ReadBin implements HasGenomeLocation {

        private final ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        private byte[] reference = null;
        private GenomeLoc loc = null;

        public ReadBin() { }

        // Return false if we can't process this read bin because the reads are not correctly overlapping.
        // This can happen if e.g. there's a large known indel with no overlapping reads.
        public void add( final SAMRecord read ) {

            final SAMRecord clippedRead = (new ReadClipper(read)).hardClipLowQualEnds( MIN_TAIL_QUALITY );
            
            if( clippedRead.getReadLength() > 0 ) {
                GenomeLoc locForRead = getToolkit().getGenomeLocParser().createGenomeLoc(clippedRead);
                if ( loc == null )
                    loc = locForRead;
                else if ( locForRead.getStop() > loc.getStop() )
                    loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), loc.getStart(), locForRead.getStop());

                reads.add(clippedRead);
            }
        }

        public List<SAMRecord> getReads() { return reads; }

        public List<SAMRecord> getReadsInWindow( final GenomeLoc window ) {
            final ArrayList<SAMRecord> readsOverlappingVariant = new ArrayList<SAMRecord>();

            for( final SAMRecord rec : reads ) {
                if( rec.getMappingQuality() > 18 && !BadMateFilter.hasBadMate(rec) ) {
                    GenomeLoc locForRead = getToolkit().getGenomeLocParser().createGenomeLoc(rec);
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
                StringUtil.toUpperCase(reference);
            }

            return reference;
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