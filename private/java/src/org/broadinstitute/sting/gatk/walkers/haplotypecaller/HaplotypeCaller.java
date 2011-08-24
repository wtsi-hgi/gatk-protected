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
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.InferredGeneticContext;
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
    protected VCFWriter writer = null;

    @Output(fullName="graphOutput", shortName="graph", doc="File to which debug assembly graph information should be written", required = true)
    protected PrintStream graphWriter = null;

    @Argument(fullName = "assembler", shortName = "assembler", doc = "Assembler to use; currently only SIMPLE_DE_BRUIJN is available.", required = false)
    protected LocalAssemblyEngine.ASSEMBLER ASSEMBLER_TO_USE = LocalAssemblyEngine.ASSEMBLER.SIMPLE_DE_BRUIJN;

    @Hidden
    @Argument(fullName = "readsToUse", shortName = "readsToUse", doc = "For debugging: how many reads to use", required = false)
    protected int numReadsToUse = -1;

    // the assembly engine
    LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    LikelihoodCalculationEngine likelihoodCalculationEngine = new LikelihoodCalculationEngine(45.0, 10.0, false, true, false);

    // the intervals input by the user
    private Iterator<GenomeLoc> intervals = null;

    // the current interval in the list
    private GenomeLoc currentInterval = null;

    // the reads that fall into the current interval
    private final ReadBin readsToAssemble = new ReadBin();

    // fasta reference reader to supplement the edges of the reference sequence
    private IndexedFastaSequenceFile referenceReader;

    // Smith-Waterman parameters copied from IndelRealigne
    private static final double SW_MATCH = 30.0;      // 1.0;
    private static final double SW_MISMATCH = -10.0;  //-1.0/3.0;
    private static final double SW_GAP = -10.0;       //-1.0-1.0/3.0;
    private static final double SW_GAP_EXTEND = -2.0; //-1.0/.0;

    // reference base padding size
    private static final int REFERENCE_PADDING = 50;

    public void initialize() {

        // get all of the unique sample names
        // if we're supposed to assume a single sample, do so
        Set<String> samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        // initialize the header
        writer.writeHeader(new VCFHeader(new HashSet<VCFHeaderLine>(), samples));

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile,ex);
        }

        assemblyEngine = makeAssembler(ASSEMBLER_TO_USE, referenceReader);

        GenomeLocSortedSet intervalsToAssemble = getToolkit().getIntervals();
        if ( intervalsToAssemble == null || intervalsToAssemble.isEmpty() )
            throw new UserException.BadInput("Intervals must be provided with -L or -BTI (preferably not larger than several hundred bp)");

        intervals = intervalsToAssemble.clone().iterator();
        currentInterval = intervals.hasNext() ? intervals.next() : null;
    }

    private LocalAssemblyEngine makeAssembler(LocalAssemblyEngine.ASSEMBLER type, IndexedFastaSequenceFile referenceReader) {
        switch ( type ) {
            case SIMPLE_DE_BRUIJN:
                return new SimpleDeBruijnAssembler(graphWriter, referenceReader, numReadsToUse);
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
            processReadBin();
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
            processReadBin();
            result++;
        }
        logger.info("Ran local assembly on " + result + " intervals");
    }

    private void processReadBin() {

        System.out.println(readsToAssemble.getLocation() + " with " + readsToAssemble.getReads().size() + " reads:");

        final List<Haplotype> haplotypes = assemblyEngine.runLocalAssembly( readsToAssemble.getReads() );
        final Pair<Haplotype, Haplotype> bestTwoHaplotypes = likelihoodCalculationEngine.computeLikelihoods( haplotypes, readsToAssemble.getReads() );

        final SWPairwiseAlignment swConsensus1 = new SWPairwiseAlignment( readsToAssemble.getReference( referenceReader ), bestTwoHaplotypes.first.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        final SWPairwiseAlignment swConsensus2 = new SWPairwiseAlignment( readsToAssemble.getReference( referenceReader ), bestTwoHaplotypes.second.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );

        System.out.println( bestTwoHaplotypes.first.toString() );
        System.out.println( "Cigar = " + swConsensus1.getCigar() );
        final List<VariantContext> vcs1 = generateVCsFromAlignment( swConsensus1, readsToAssemble.getReference( referenceReader ), bestTwoHaplotypes.first.bases, readsToAssemble.getLocation() );

        System.out.println( bestTwoHaplotypes.second.toString() );
        System.out.println( "Cigar = " + swConsensus2.getCigar() );
        final List<VariantContext> vcs2 = generateVCsFromAlignment( swConsensus2, readsToAssemble.getReference( referenceReader ), bestTwoHaplotypes.second.bases, readsToAssemble.getLocation() );

        final List<VariantContext> finalVCs = genotype( vcs1, vcs2 );
        for( final VariantContext vc : finalVCs ) {
            System.out.println(vc);
            writer.add(vc);
        }
        System.out.println("----------------------------------------------------------------------------------");

    }

    private static List<VariantContext> generateVCsFromAlignment( final SWPairwiseAlignment swConsensus, final byte[] ref, final byte[] read, final GenomeLoc loc ) {
        final ArrayList<VariantContext> vcs = new ArrayList<VariantContext>();

        int refPos = swConsensus.getAlignmentStart2wrt1();
        int readPos = 0;
        final int lookAhead = 5;

        for( final CigarElement ce : swConsensus.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
                case I:
                {
                    byte[] insertionBases = Arrays.copyOfRange( read, readPos, readPos + elementLength);
                    boolean allN = true;
                    for( byte b : insertionBases ) {
                        if( b != (byte) 'N' ) {
                            allN = false;
                            break;
                        }
                    }
                    if( !allN ) {
                        ArrayList<Allele> alleles = new ArrayList<Allele>();
                        alleles.add( Allele.create(Allele.NULL_ALLELE_STRING, true));
                        alleles.add( Allele.create(insertionBases, false));
                        System.out.println("Insertion: " + alleles);
                        vcs.add(new VariantContext("HaplotypeCaller", loc.getContig(), loc.getStart() + refPos - 1, loc.getStart() + refPos - 1, alleles, VariantContext.NO_GENOTYPES, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null, ref[refPos-1]));
                    }
                    readPos += elementLength;
                    break;
                }
                case S:
                {
                    readPos += elementLength;
                    refPos += elementLength;
                    break;
                }
                case D:
                {
                    byte[] deletionBases = Arrays.copyOfRange( ref, refPos, refPos + elementLength);
                    ArrayList<Allele> alleles = new ArrayList<Allele>();
                    alleles.add( Allele.create(deletionBases, true) );
                    alleles.add( Allele.create(Allele.NULL_ALLELE_STRING, false) );
                    System.out.println( "Deletion: " + alleles);
                    vcs.add( new VariantContext("HaplotypeCaller", loc.getContig(), loc.getStart() + refPos - 1, loc.getStart() + refPos + elementLength - 1, alleles, VariantContext.NO_GENOTYPES, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null, ref[refPos-1]) );
                    refPos += elementLength;
                    break;
                }
                case M:
                {
                    int numSinceMismatch = -1;
                    int stopOfMismatch = -1;
                    int startOfMismatch = -1;
                    int refPosStartOfMismatch = -1;
                    for( int iii = 0; iii < elementLength; iii++ ) {
                        if( ref[refPos] != read[readPos] ) {
                            // SNP or start of possible MNP
                            if( stopOfMismatch == -1 ) {
                                startOfMismatch = readPos;
                                stopOfMismatch = readPos;
                                numSinceMismatch = 0;
                                refPosStartOfMismatch = refPos;
                            } else {
                                stopOfMismatch = readPos;
                            }
                        }

                        if( stopOfMismatch != -1) {
                            numSinceMismatch++;
                        }

                        if( numSinceMismatch > lookAhead || (iii == elementLength - 1 && stopOfMismatch != -1) ) {
                            byte[] refBases = Arrays.copyOfRange( ref, refPosStartOfMismatch, refPosStartOfMismatch + (stopOfMismatch - startOfMismatch) + 1 );
                            byte[] mismatchBases = Arrays.copyOfRange( read, startOfMismatch, stopOfMismatch + 1 );
                            ArrayList<Allele> alleles = new ArrayList<Allele>();
                            alleles.add( Allele.create( refBases, true ) );
                            alleles.add( Allele.create( mismatchBases, false ) );
                            System.out.println( "SNP/MNP: " + alleles);
                            vcs.add( new VariantContext("HaplotypeCaller", loc.getContig(), loc.getStart() + refPosStartOfMismatch, loc.getStart() + refPosStartOfMismatch + (stopOfMismatch - startOfMismatch), alleles) );
                            numSinceMismatch = -1;
                            stopOfMismatch = -1;
                            startOfMismatch = -1;
                            refPosStartOfMismatch = -1;
                        }

                        refPos++;
                        readPos++;
                    }
                    break;
                }

                case N:
                case H:
                case P:
                default:
                    throw new ReviewedStingException( "Unsupported cigar operator: " + ce.getOperator() );
            }
        }

        return vcs;
    }

    private static List<VariantContext> genotype( final List<VariantContext> vcs1, final List<VariantContext> vcs2 ) {
        final ArrayList<VariantContext> vcs = new ArrayList<VariantContext>();

        final Iterator<VariantContext> vcs1Iter = vcs1.iterator();
        final Iterator<VariantContext> vcs2Iter = vcs2.iterator();

        VariantContext vc1Hold = null;
        VariantContext vc2Hold = null;

        do {
            final VariantContext vc1 = ( vc1Hold != null ? vc1Hold : (vcs1Iter.hasNext() ? vcs1Iter.next() : null) );
            final VariantContext vc2 = ( vc2Hold != null ? vc2Hold : (vcs2Iter.hasNext() ? vcs2Iter.next() : null) );

            vc1Hold = null;
            vc2Hold = null;


            if( vc1 == null && vc2 != null ) {
                ArrayList<Allele> alleles = new ArrayList<Allele>();
                alleles.addAll( vc2.getAlleles() );
                Genotype gt = new Genotype( "NA12878", alleles );
                HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                genotypeMap.put("NA12878", gt);
                vcs.add( VariantContext.modifyGenotypes( vc2, genotypeMap ) );
            } else if( vc1 != null && vc2 == null ) {
                ArrayList<Allele> alleles = new ArrayList<Allele>();
                alleles.addAll( vc1.getAlleles() );
                Genotype gt = new Genotype( "NA12878", alleles );
                HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                genotypeMap.put("NA12878", gt);
                vcs.add( VariantContext.modifyGenotypes( vc1, genotypeMap ) );
            } else if( vc1 != null ) { // && vc2 != null
                if( vc1.getStart() == vc2.getStart() ) {
                    ArrayList<Allele> alleles = new ArrayList<Allele>();
                    alleles.add( vc1.getAlternateAllele(0) );
                    alleles.add( vc2.getAlternateAllele(0) );
                    Genotype gt = new Genotype( "NA12878", alleles );
                    HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                    genotypeMap.put("NA12878", gt);
                    vcs.add( VariantContext.modifyGenotypes( vc1, genotypeMap ) );
                } else if( vc1.getStart() < vc2.getStart()) {
                    vc2Hold = vc2;
                    ArrayList<Allele> alleles = new ArrayList<Allele>();
                    alleles.addAll( vc1.getAlleles() );
                    Genotype gt = new Genotype( "NA12878", alleles );
                    HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                    genotypeMap.put("NA12878", gt);
                    vcs.add( VariantContext.modifyGenotypes( vc1, genotypeMap ) );
                } else {
                    vc1Hold = vc1;
                    ArrayList<Allele> alleles = new ArrayList<Allele>();
                    alleles.addAll( vc2.getAlleles() );
                    Genotype gt = new Genotype( "NA12878", alleles );
                    HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                    genotypeMap.put("NA12878", gt);
                    vcs.add( VariantContext.modifyGenotypes( vc2, genotypeMap ) );
                }
            }


        } while ( vcs1Iter.hasNext() || vcs2Iter.hasNext() );

        return vcs;
    }

    // private class copied from IndelRealigner
    private class ReadBin implements HasGenomeLocation {

        private final ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        private byte[] reference = null;
        private GenomeLoc loc = null;

        public ReadBin() { }

        // Return false if we can't process this read bin because the reads are not correctly overlapping.
        // This can happen if e.g. there's a large known indel with no overlapping reads.
        public void add(SAMRecord read) {

            GenomeLoc locForRead = getToolkit().getGenomeLocParser().createGenomeLoc(read);
            if ( loc == null )
                loc = locForRead;
            else if ( locForRead.getStop() > loc.getStop() )
                loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), loc.getStart(), locForRead.getStop());

            reads.add(read);
        }

        public List<SAMRecord> getReads() { return reads; }

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