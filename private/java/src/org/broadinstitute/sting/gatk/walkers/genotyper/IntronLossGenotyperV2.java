package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMUtils;
import net.sf.samtools.util.StringUtil;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportColumn;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.misc.intronloss.IntronLossGenotypeLikelihoodCalculationModel;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: loaner
 * Date: 8/20/11
 * Time: 10:35 PM
 * To change this template use File | Settings | File Templates.
 */
@ReadFilters({DuplicateReadFilter.class,FailsVendorQualityCheckFilter.class,MappingQualityZeroFilter.class})
public class IntronLossGenotyperV2 extends ReadWalker<SAMRecord,Integer> {

    /**
     * A raw, unfiltered, highly specific callset in VCF format.
     */
    @Output(doc="File to which variants should be written", required = true)
    protected VCFWriter vcfWriter = null;

    @Input(shortName="r",fullName="refSeq",required=true,doc="The RefSeq Gene definition track")
    public RodBinding<RefSeqFeature> refSeqRodBinding;

    @Input(shortName="H",fullName="insertHistogram",required=true,doc="The insert size histogram per read group, either flat file or GATK report formatted")
    public File readGroupInsertHistogram;

    @Hidden
    @Input(fullName="maxInsertSize",required=false)
    public int MAX_INSERT_SIZE = 1000;

    public final static byte MIN_TAIL_QUALITY = 6;

    private Map<String,byte[]> insertQualsByRG = new HashMap<String,byte[]>();

    private boolean initialized = false;

    // the likelihoods engine
    IntronLossGenotypeLikelihoodCalculationModel ilglcm= new IntronLossGenotypeLikelihoodCalculationModel(logger);

    // the reads that fall into the current interval
    private final Map<RefSeqFeature,PairedReadBin> readsToAssemble = new TreeMap<RefSeqFeature,PairedReadBin>(new Comparator<RefSeqFeature>() {
        @Override
            public int compare(RefSeqFeature refSeqFeature, RefSeqFeature refSeqFeature1) {
               GenomeLoc featureStop = getToolkit().getGenomeLocParser().createGenomeLoc(refSeqFeature).getStopLocation();
               GenomeLoc feature1Stop = getToolkit().getGenomeLocParser().createGenomeLoc(refSeqFeature1).getStopLocation();
                return featureStop.compareTo(feature1Stop);  //To change body of implemented methods use File | Settings | File Templates.
            }
    });

    // fasta reference reader to supplement the edges of the reference sequence
    private IndexedFastaSequenceFile referenceReader;

    // reference base padding size
    private static final int REFERENCE_PADDING = 10000;

    @Override
    public void initialize() {
        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile,ex);
        }

        Set<Sample> samples = getToolkit().getSAMFileSamples();
        Set<String> sampleStr = new HashSet<String>(samples.size());
        for ( Sample s : samples ) {
            sampleStr.add(s.getId());
        }

        ilglcm.setSamples(sampleStr);

        vcfWriter.writeHeader(new VCFHeader(new HashSet<VCFHeaderLine>(), sampleStr));

        try {

            XReadLines xrl = new XReadLines(readGroupInsertHistogram);
            if ( ! xrl.next().startsWith("##:") ) {
                xrl.close();
                for ( String entry : new XReadLines(readGroupInsertHistogram) ) {
                    String[] split1 = entry.split("\\t");
                    String id = split1[0];
                    String[] histogram = split1[1].split(";");
                    byte[] quals = new byte[histogram.length];
                    int idx = 0;
                    for ( String histEntry : histogram ) {
                        try {
                            quals[idx++] = Byte.parseByte(histEntry);
                        } catch( NumberFormatException e) {
                            quals[idx-1] = QualityUtils.probToQual(Double.parseDouble(histEntry));
                        }
                    }

                    insertQualsByRG.put(id,quals);
                }
            } else {
                xrl.close();
                GATKReport report = new GATKReport(readGroupInsertHistogram);
                GATKReportTable reportTable = report.getTable("InsertSizeDistribution");
                // rows are insert sizes, columns are read groups
                for (GATKReportColumn reportColumn : reportTable.getColumns() ) {
                    // annoyingly, the column has no knowledge of its own rows
                    int sum = 0;
                    for ( int row = 0; row < reportTable.getNumRows(); row++ ) {
                        sum += Integer.parseInt( (String) reportTable.get(row,reportColumn.getColumnName()));
                    }
                    byte[] rgHist = new byte[MAX_INSERT_SIZE];
                    for ( int row = 0; row < rgHist.length; row++) {
                        int val = 1;
                        if ( reportTable.containsKey(row) ) {
                            val = Integer.parseInt( (String) reportTable.get(row,reportColumn.getColumnName()));
                        }
                        rgHist[row] = QualityUtils.probToQual( 1.0-( ( (double) val )/sum ), Math.pow(10,-25.4) );
                    }

                    insertQualsByRG.put(reportColumn.getColumnName(),rgHist);
                }
            }

        } catch (FileNotFoundException e) {
            throw new UserException("Histogram file not found",e);
        } catch (IOException e) {
            throw new StingException("IO Exception",e);
        }

        StringBuffer debug = new StringBuffer();
        for ( Map.Entry<String,byte[]> etry : insertQualsByRG.entrySet() ) {
            debug.append("  ");
            debug.append(etry.getKey());
            debug.append("->");
            debug.append(Arrays.deepToString(ArrayUtils.toObject(etry.getValue())));
        }

        logger.debug(debug);

        ilglcm.initialize(refSeqRodBinding,getToolkit().getGenomeLocParser(),insertQualsByRG);
    }

    public SAMRecord map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        boolean overlapsExon = false;
        List<RefSeqFeature> refSeqFeatures = new ArrayList<RefSeqFeature>(16);
        for (GATKFeature feature : metaDataTracker.getAllCoveringRods() ) {
            if ( feature.getUnderlyingObject().getClass().isAssignableFrom(RefSeqFeature.class) ) {
                refSeqFeatures.add((RefSeqFeature) feature.getUnderlyingObject());
            }
        }

        for ( RefSeqFeature refSeqFeature : refSeqFeatures ) {
            if ( ! readsToAssemble.containsKey(refSeqFeature) ) {
                readsToAssemble.put(refSeqFeature,new PairedReadBin());
            }
            overlapsExon |= refSeqFeature.overlapsExonP(getToolkit().getGenomeLocParser().createGenomeLoc(read));
        }

        Set<RefSeqFeature> toRem = new TreeSet<RefSeqFeature>(new Comparator<RefSeqFeature>() {
            @Override
            public int compare(RefSeqFeature refSeqFeature, RefSeqFeature refSeqFeature1) {
               GenomeLoc featureStop = getToolkit().getGenomeLocParser().createGenomeLoc(refSeqFeature).getStopLocation();
               GenomeLoc feature1Stop = getToolkit().getGenomeLocParser().createGenomeLoc(refSeqFeature1).getStopLocation();
                return featureStop.compareTo(feature1Stop);  //To change body of implemented methods use File | Settings | File Templates.
            }
        });

        for ( RefSeqFeature refSeqFeature : readsToAssemble.keySet() ) {
            GenomeLoc featureStop = getToolkit().getGenomeLocParser().createGenomeLoc(refSeqFeature);
            if ( ref.getLocus().getStartLocation().isPast(featureStop) ) {
                toRem.add(refSeqFeature);
            }
        }

        for ( RefSeqFeature rsf : toRem ) {
            logger.debug(String.format("Processing %s at %s",rsf.getGeneName(),getToolkit().getGenomeLocParser().createGenomeLoc(rsf).toString()));
            processReadBin(rsf, readsToAssemble.get(rsf));
            readsToAssemble.get(rsf).clear();
            readsToAssemble.remove(rsf);
        }

        if ( ! overlapsExon ) {
            return null;
        }

        int start = read.getUnclippedStart();
        int stop= read.getUnclippedEnd();

        int leftHardClipAlready = read.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.H) ? read.getCigar().getCigarElement(0).getLength() : 0;
        byte[] refBases = getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(read.getReferenceName(),start, stop).getBases();
        SAMRecord clippedRead = ReadUtils.unclipSoftClippedBases((new ReadClipper(read)).hardClipLowQualEnds((byte) 12));

        if( clippedRead.getCigar().isEmpty()) {
            return null;
        }
        int offset = clippedRead.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.H) ? clippedRead.getCigar().getCigarElement(0).getLength() - leftHardClipAlready : 0;
        clippedRead.setAttribute("NM",AlignmentUtils.getMismatchCount(clippedRead,refBases,offset).numMismatches);

        return clippedRead;
    }

    public boolean filter(ReferenceContext ref, SAMRecord read) {
        // kill off unmapped reads
        boolean filter = ! read.getReadUnmappedFlag() && read.getMappingQuality() > 0 && ((( read.getAttribute("MQ") != null && ((Integer)read.getAttribute("MQ")) > 0)));
        return filter;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(SAMRecord read, Integer sum) {
        if ( read == null )
            return sum;
        for ( PairedReadBin rb : readsToAssemble.values() ) {
            rb.add(read);
        }
        return ++sum;
    }

    public void onTraversalDone(Integer result) {
        if ( readsToAssemble.size() > 0 ) {
            for ( Map.Entry<RefSeqFeature,PairedReadBin> refSeqEntry : readsToAssemble.entrySet() ) {
                processReadBin(refSeqEntry.getKey(),refSeqEntry.getValue());
                result++;
            }
        }
        logger.info("Ran ILG on " + result + " reads");
    }

    private void processReadBin(final RefSeqFeature geneFeature, PairedReadBin readBin) {
        logger.debug("Processing read bin");

       logger.info(readBin.getLocation() + " with " + readBin.getReadPairs().size() + " reads:");

        List<VariantContext> likVC = ilglcm.getLikelihoods(readBin,geneFeature,referenceReader);

        for ( VariantContext vc : likVC ) {
            vcfWriter.add(vc);
        }
    }

    // private class copied from IndelRealigner, used to bin together a bunch of reads and then retrieve the reference overlapping the full extent of the bin
    public class PairedReadBin implements HasGenomeLocation, Iterable<Pair<SAMRecord,SAMRecord>> {

        private final HashMap<String,Pair<SAMRecord,SAMRecord>> reads = new HashMap<String, Pair<SAMRecord, SAMRecord>>(3200);
        private byte[] reference = null;
        private GenomeLoc loc = null;

        public PairedReadBin() { }

        // Return false if we can't process this read bin because the reads are not correctly overlapping.
        // This can happen if e.g. there's a large known indel with no overlapping reads.
        public void add(SAMRecord read) {
            if ( read.getCigarString().equals("*") ) {
                return;
            }
            GenomeLoc locForRead = getToolkit().getGenomeLocParser().createGenomeLoc(read);
            if ( loc == null )
                loc = locForRead;
            else if ( locForRead.getStop() > loc.getStop() )
                loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), loc.getStart(), locForRead.getStop());

            String name = read.getReadName();
            if ( reads.containsKey(name) ) {
                // mate is already in the map
                reads.get(name).second = read;
            } else {
                Pair<SAMRecord,SAMRecord> samPair = new Pair<SAMRecord,SAMRecord>(read,null);
                reads.put(name,samPair);
            }
        }

        public Collection<Pair<SAMRecord,SAMRecord>> getReadPairs() { return reads.values(); }

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

        public Iterator<Pair<SAMRecord,SAMRecord>> iterator() {
            return  reads.values().iterator();
        }
    }
}
