package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 8/3/11
 * Time: 2:24 PM
 * To change this template use File | Settings | File Templates.
 */
public class SlidingWindow {
    protected final static String RG_POSTFIX = ".ReducedReads";

    public LinkedList<SlidingRead> SlidingReads = new LinkedList<SlidingRead>();
    public LinkedList<CountWithBase> countsWithBases = new LinkedList<CountWithBase>();
    private String contig = null;
    private int start; //of the first read
    SAMFileHeader header;
    SAMRecord runningConsensus;
    final SAMReadGroupRecord reducedReadGroup;
    int consensusCounter = 0;
    private int contigIndex;

    public SlidingWindow(final String sampleName, String Contig, SAMFileHeader Header) {
        this.contig = Contig;
        //this.contigIndex = ContigIndex;
        this.header = Header;
        runningConsensus = null;
        this.reducedReadGroup = createReducedReadGroup(sampleName);
    }

     /**
     * Helper function to create a read group for these reduced reads
     * @param sampleName
     * @return
     */
    private static final SAMReadGroupRecord createReducedReadGroup(final String sampleName) {
        SAMReadGroupRecord rg = new SAMReadGroupRecord(sampleName + RG_POSTFIX);
        rg.setSample(sampleName);
        return rg;
    }

    public void slide(int position) {

        //position is the start of the new window


        LinkedList<SlidingRead> slidReads = new LinkedList<SlidingRead>();
        for ( SlidingRead read: SlidingReads ) {
            if ( read.getAlignmentStart() < position )
                read = read.clipStart(position);
            if ( !read.getBasesAndQuals().isEmpty())
                slidReads.add(read);
                // TODO there might be a better way
        }
        if (slidReads.isEmpty())
            countsWithBases.clear();
        else {
            while ( countsWithBases.getFirst().location < position)
                countsWithBases.pop();
        }
        SlidingReads = slidReads;
    }

    public int getEnd() {
        if( !(SlidingReads.isEmpty()) )
            return SlidingReads.getLast().getAlignmentStop();
        else
            return -1;

    }

    protected class CountWithBase {
        BaseCounts counts;
        int location;
        boolean isVariant;

        public CountWithBase(byte base, byte qual, int Location) {
            this.counts = new BaseCounts();
            if (qual >= 20)
                this.counts.incr(base);
            this.location = Location;
            this.isVariant = false;
        }

        public CountWithBase(int Location) {
            this.counts = new BaseCounts();
            this.location = Location;
            this.isVariant = false;
        }
        // TODO minMapQual filters and minBaseQual filters
        public boolean addBase(byte base, byte qual) {
            // return true if a variant site was CREATED
            boolean result = false;
            if ( qual >= 20 ) {
                if (!this.isVariant) {
                    if (base != this.counts.baseWithMostCounts() && this.counts.totalCount() != 0 ) {
                        this.isVariant = true;
                        result = true;
                    }
                }
                this.counts.incr(base);
            }
            return result;

        }

    }

    public boolean addRead( SAMRecord read ) {
        // Assuming reads are ordered, should return variant location if created variance
        boolean result = addToCounts(read);
        SlidingReads.add(new SlidingRead(read));
        return result;
    }

    public int getStart() {
        if ( !countsWithBases.isEmpty() )
            start = countsWithBases.getFirst().location;
        else
            start = 0;
        return start;
    }

    // look through the read and slidingReads
    // return true if a variant site was created variance
    // add read to counts index
    private boolean addToCounts(SAMRecord read) {
        boolean result = false;
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();
        if (SlidingReads.isEmpty()) {
            // update contigs
            contig = read.getReferenceName();
            contigIndex = read.getReferenceIndex();
            header = read.getHeader();
            // TODO add Quality check for low qual scores here

            for (int i = 0 ; i < read.getReadLength(); i++) {
                countsWithBases.add( new CountWithBase( bases[i], quals[i], i + read.getAlignmentStart()) );
            }
            return result;
        }
        else {
            // do this while there is a slidingWindow position in the linked list
            // we have to make new elements for the rest.

            int i = 0;

            int index = read.getAlignmentStart() - countsWithBases.getFirst().location;

            Iterator<CountWithBase> I = countsWithBases.iterator();
            CountWithBase cBase = I.next();
            // move by the indexed amount
            for ( int j = 0; (j < index) && ( I.hasNext() ); j++ )
                cBase = I.next();
            /*
            while ( j < index  ) {
                countsWithBases.add( new CountWithBase(countsWithBases.getFirst().location + j));
                cBase = I.next();
                j++;
            }
            */
            // increment elements while they exist in the window, AND while you have reads
            // TODO FIX: If I does not have next, cBase is last element, but it never gets added
            while (  i < read.getReadLength()) {

                // If variant was created
                if (cBase.addBase(bases[i], quals[i]))
                    result = true;
                i++;
                if (I.hasNext())
                    cBase = I.next();
                else
                    break;
            }
            // create new elements
            while (i < read.getReadLength() ) {
                //add new  element to sliding window
                countsWithBases.add( new CountWithBase( bases[i], quals[i], i + read.getAlignmentStart()) );
                i++;
            }
        }
        return result;
    }

    // gets the region(s) of variable sites present in the window
    public List<VariableRegion> getVariableRegions(int contextSize) {
        // TODO require contextSize be non-Negative

        List<VariableRegion> rawRegions = new LinkedList<VariableRegion>();

        // Create every variant site
        for (CountWithBase cBase: countsWithBases) {
            if ( cBase.isVariant ) {
                rawRegions.add(createVariableRegion(cBase, contextSize));
            }
        }

        if (rawRegions.isEmpty())
            return rawRegions;

        // Merge the proximate sites
        List<VariableRegion> regions = new LinkedList<VariableRegion>();
        Iterator<VariableRegion> i = rawRegions.iterator();
        VariableRegion p = i.next();
        while ( i.hasNext() ) {
            VariableRegion q = i.next();
            while ( p.end >= q.start ) {
                p = p.merge(q);
                if (i.hasNext() )
                    q = i.next();
                else
                    break;
            }
            regions.add(p);
            p = q;
        }

        if (regions.isEmpty())
            regions.add(p);
        //for ( VariableRegion vr: regions )
            //System.out.println(String.format("INFO -- #### -- VARIANT FOUND Creating Variable Region, %d - %d", vr.start, vr.end));

        return regions;
    }

    private VariableRegion createVariableRegion(CountWithBase cBase, int ContextSize) {
        return new VariableRegion(cBase.location - ContextSize, cBase.location + ContextSize);
    }

    public List<SAMRecord> finalizeVariableRegion(VariableRegion variableRegion) {
        //System.out.println(String.format(String.format("INFO ######################### \n \n Finalizing Variable Region: %d-%d ", variableRegion.start, variableRegion.end)));
        List<SAMRecord> output = new LinkedList<SAMRecord>();

        for ( SlidingRead read: SlidingReads ) {
            SAMRecord SAM = read.trimToVariableRegion(variableRegion);
            SAM.setReadName(SAM.getReadName()+".trim");
            if ( SAM.getReadLength() > 0 ) {
                output.add(SAM);
                //System.out.println(String.format("Output Variable Read: %d-%d", SAM.getAlignmentStart(), SAM.getAlignmentEnd()));
            }
        }
        slide(variableRegion.end+1);
        createRunningConsensus();
        return output;
    }

    private void createRunningConsensus() {

        //header stuff
        runningConsensus = new SAMRecord(header);
        //runningConsensus.setAttribute("RG", reducedReadGroup.getId());
        //runningConsensus.setAttribute(ReadUtils.REDUCED_READ_QUALITY_TAG, Integer.valueOf(REDUCED_READ_BASE_QUALITY)); //Also problems: Should it be avg
        //runningConsensus.setAttribute("QE", QualityEquivalent);    // Qual equivs
        runningConsensus.setReferenceName(contig);
        runningConsensus.setReferenceIndex(contigIndex);
        //runningConsensus.setReadName(String.format("%s.read.%d", reducedReadGroup.getId(), consensusCounter++));
        runningConsensus.setReadPairedFlag(false);
        runningConsensus.setReadUnmappedFlag(false);
        //runningConsensus.setCigarString(String.format("%dM", size));    // TODO fix cigar
        runningConsensus.setAlignmentStart(getStart());

    }

    public List<SAMRecord> finalizeConsensusRead(VariableRegion variableRegion) {
        List<SAMRecord> result = new LinkedList<SAMRecord>();
        addToConsensus(variableRegion.start);
        SAMRecord consensus = runningConsensus;
        // This determines the end of read
        if ( runningConsensus != null ){
            if ( consensus.getReadLength() != 0 ) {
                consensus.setMappingQuality(60); // TODO set up rms
                consensus.setCigarString(String.format("%dM", consensus.getReadLength()));
                consensus.setReadName(String.format("%s.read.%d", reducedReadGroup.getId(), consensusCounter++));
                result.add(consensus);
            }
        }
        runningConsensus = null;

        return result;
    }

    public void addToConsensus(int position) {
        // compresses the sliding reads to a streaming(incomplete) SamRecord
        if ( position != getStart() ) {
            SAMRecord consensus;
            /*
            if ( position == -1 )
                position = (getEnd() + 1);
            */
            int size = position - getStart();

            if ( size >= 0 ) {

                //System.out.println(String.format("INFO -- Compressing running Consensus from %d to %d  ", getStart(), position));

                byte[] newBases = new byte[size];
                byte[] newQuals = new byte[size];
                if ( !countsWithBases.isEmpty() ) {
                    Iterator<SlidingWindow.CountWithBase> I = countsWithBases.iterator();
                    SlidingWindow.CountWithBase cBase = I.next();
                    for ( int i = 0; (i < size) && I.hasNext(); i++ ) {
                        BaseCounts Counts= cBase.counts;
                        newBases[i] = Counts.baseWithMostCounts();
                        newQuals[i] = QualityUtils.boundQual(Counts.countOfMostCommonBase(), (byte) 64);
                        cBase = I.next();
                    }
                }
                if ( runningConsensus == null ) {
                    createRunningConsensus();
                    //header stuff
                    consensus = runningConsensus;

                    consensus.setReadBases(newBases);
                    consensus.setBaseQualities(newQuals);
                }
                else {
                    consensus = runningConsensus;
                    //size += consensus.getReadLength();

                    byte[] bases = consensus.getReadBases();
                    byte[] quals = consensus.getBaseQualities();

                    bases = ArrayUtils.addAll(bases, newBases);
                    quals = ArrayUtils.addAll(quals, newQuals);

                    consensus.setReadBases(bases);
                    consensus.setBaseQualities(quals);
                }
                // TODO MQ RMS calc
                // Once the new consensus is generated, we need to slide the window
                slide(position);


                runningConsensus = consensus;
            }
        }
    }

}

