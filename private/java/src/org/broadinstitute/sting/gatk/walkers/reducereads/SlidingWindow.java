package org.broadinstitute.sting.gatk.walkers.reducereads;

import com.google.java.contract.Requires;
import net.sf.samtools.*;

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

    private LinkedList<SlidingRead> SlidingReads = new LinkedList<SlidingRead>();
    private LinkedList<CountWithBase> countsWithBases = new LinkedList<CountWithBase>();
    private int contextSize;
    private String contig;
    private int startLocation;
    private int stopLocation;

    private RunningConsensus runningConsensus;

    private int contigIndex;
    private final double MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT = 0;   // proportion has to be greater than this value
    private final int MIN_BASE_QUAL_TO_COUNT = 20;                         // qual has to be greater than or equal to this value
    private final int MAX_QUAL_COUNT = 64;


    protected class CountWithBase {
        private BaseCounts counts;
        private int insertionsToTheRight;
        private int location;

        public CountWithBase () {
            this.counts = new BaseCounts();
            this.insertionsToTheRight = 0;
            this.location = 0;
        }

        public CountWithBase (int location) {
            this();
            this.location = location;
        }

        public boolean isVariant() {
            return counts.totalCount() > 1 && ( isVariantFromInsertions() || isVariantFromMismatches() );
        }

        public void addBase(byte base, byte qual, boolean hasInsertionToTheRight) {
            if ( qual >= MIN_BASE_QUAL_TO_COUNT ) {
                this.counts.incr(base);
                if (hasInsertionToTheRight)
                    insertionsToTheRight++;
            }
        }

        private boolean isVariantFromInsertions() {
            return ((double) insertionsToTheRight / counts.totalCount()) > MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT;
        }

        private boolean isVariantFromMismatches() {
            return counts.baseCountProportion(counts.baseWithMostCounts()) < (1 - MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT);
        }
    }

    public SlidingWindow(String contig, int contextSize) {
        this.contig = contig;
        this.runningConsensus = null;
        this.startLocation = -1;
        this.stopLocation = -1;
        this.contextSize = contextSize;
    }

    public int getStartLocation() {
        return startLocation;
    }

    public int getStopLocation() {
        return stopLocation;
    }

    public List<SAMRecord> addRead( SAMRecord read ) {
        List<SAMRecord> finalizedReads = slideIfPossible(read.getUnclippedStart());

        addToCounts(read);

        SlidingReads.add(new SlidingRead(read));

        return finalizedReads;
    }

    /**
     * Determines if the window can be slid given the new incoming read.
     *
     * We check from the start of the window to the (unclipped) start of the new incoming read if there
     * is any variant.
     *
     * If there are no variant sites, we slide the left side of the window to the unclipped start of the
     * new incoming read minus the context size (we can only guarantee that the reads before this point
     * will not be in any variant region ever)
     *
     * If there are variant sites, we check if it's time to close the variant region.
     *
     * @param incomingReadUnclippedStart the incoming read's start position. Must be the unclipped start!
     * @return all reads that have fallen to the left of the sliding window after the slide
     */
    private List<SAMRecord> slideIfPossible (int incomingReadUnclippedStart) {

        List<SAMRecord> finalizedReads = new LinkedList<SAMRecord>();

        // No point doing any calculation if the read's unclipped start is behind
        // the start of the window
        if (incomingReadUnclippedStart > startLocation) {
            boolean [] variantSite = markSites(incomingReadUnclippedStart);

            // this is the limit of what we can close/send to consensus (non-inclusive)
            int breakpoint = incomingReadUnclippedStart - contextSize;
            int start = 0;
            int i = 0;

            while (i < breakpoint) {
                while (i<breakpoint && !variantSite[i]) i++;
                addToConsensus(start, i);

                start = i;
                while(i<breakpoint && variantSite[i]) i++;
                if (i < breakpoint || i == breakpoint && !variantSite[i])
                   finalizedReads.addAll(closeVariantRegion(start, i));
            }

            slide(start);
        }


        return finalizedReads;
    }

    /**
     * returns an array marked with variant and non-variant regions
     *
     * @param upto check the window from start to 'upto'
     * @return a boolean array with 'true' marking variant regions and false marking consensus sites
     */
    private boolean [] markSites(int upto) {

        boolean [] markedSites = new boolean [upto - startLocation];

        Iterator<CountWithBase> countWithBaseIterator = countsWithBases.iterator();
        for (int i = startLocation; i < upto; i++) {
            if (countWithBaseIterator.hasNext()) {
                CountWithBase cb = countWithBaseIterator.next();
                if (cb.isVariant())
                    markVariantRegion(markedSites, i - startLocation);
            }
            else
                break;
        }
        return markedSites;
    }

    /**
     * performs the variant region marks (true) around a variant site
     *
     * @param markedSites the boolean array to bear the marks
     * @param variantSiteLocation the location where a variant site was found
     */
    private void markVariantRegion (boolean [] markedSites, int variantSiteLocation) {
        int from = (variantSiteLocation < contextSize) ? 0 : variantSiteLocation - contextSize;
        int to = (variantSiteLocation + contextSize + 1 > markedSites.length) ? markedSites.length : variantSiteLocation + contextSize + 1;
        for (int i=from; i<to; i++)
            markedSites[i] = true;
    }

    private void addToConsensus (int end) {
        addToConsensus(0, end);
    }

    private void addToConsensus (int start, int end) {
        if (start == end)
            return;

        int i = 0;
        for (CountWithBase cwb : countsWithBases) {
            if (i == end)
                break;
            if (i >= start) {
                runningConsensus.add(cwb.counts.baseWithMostCounts(), (byte) Math.min(cwb.counts.countOfMostCommonBase(), MAX_QUAL_COUNT));
            }
        }
    }

    private List<SAMRecord> closeVariantRegion(int start, int end) {
        List<SAMRecord> finalizedReads = new LinkedList<SAMRecord>();

        if (start < end) {
            if (runningConsensus != null)
                finalizedReads.add(runningConsensus.close());

            // Clipping operations are reference based, not read based
            int refStart = countsWithBases.get(start).location;
            int refEnd = countsWithBases.get(end).location;

            for ( SlidingRead read: SlidingReads ) {
                SAMRecord SAM = read.trimToVariableRegion(refStart, refEnd);
                SAM.setReadName(SAM.getReadName()+".trim");
                if ( SAM.getReadLength() > 0 ) {
                    finalizedReads.add(SAM);
                    //System.out.println(String.format("Output Variable Read: %d-%d", SAM.getAlignmentStart(), SAM.getAlignmentEnd()));
                }
            }
        }
        return finalizedReads;
    }


    // TODO -- CONTINUE REFACTOR FROM HERE



    /**
     * Slides the window to the new position
     *
     * @param leftPosition the new starting location of the window
     */
    private void slide (int leftPosition) {
        //position is the start of the new window
        LinkedList<SlidingRead> slidReads = new LinkedList<SlidingRead>();
        for ( SlidingRead read: SlidingReads ) {
            if ( read.getAlignmentStart() < leftPosition )
                read = read.clipStart(leftPosition);
            if ( !read.getBasesAndQuals().isEmpty())
                slidReads.add(read);
                // TODO there might be a better way
        }
        if (slidReads.isEmpty())
            countsWithBases.clear();
        else {
            while ( !countsWithBases.isEmpty() && (countsWithBases.getFirst().location < leftPosition) )
                countsWithBases.pop();
        }
        SlidingReads = slidReads;
    }

    public List<SAMRecord> close() {
//        // TODO fix ending at variable region, it must complete the consensus
//
//        LinkedList<SAMRecord> result = new LinkedList<SAMRecord>();
//        for ( VariableRegion variableRegion : getVariableRegions(contextSize) ) {
//            result.addAll(finalizeVariableRegion(variableRegion));
//
//        }
//        result.addAll(addToConsensus(stopLocation + 1));
//        result.addAll(finalizeConsensusRead());
//        return result;
        return null;
    }

    // look through the read and slidingReads
    // return true if a variant site was created variance
    // add read to counts index
    @Requires("read.getAlignmentStart >= startLocation")
    private void addToCounts(SAMRecord read) {
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();
        Cigar cigar = read.getCigar();

        // todo -- this is definitely out of place. Sliding Window should be initialized before we try to add things to it.
        // set up the Sliding Window parameters
        if (SlidingReads.isEmpty()) {
            // update contigs
            contig = read.getReferenceName();
            contigIndex = read.getReferenceIndex();
            header = read.getHeader();
            startLocation = read.getAlignmentStart();
            stopLocation = read.getAlignmentEnd();
            for (int i = startLocation; i<=stopLocation; i++)
                countsWithBases.add(new CountWithBase(read.getAlignmentStart() + i));
        }

        int readBaseIndex = 0;
        int locationIndex = read.getAlignmentStart() - startLocation;

        // Do we need to add extra elements to the end of the list?
        if (stopLocation < read.getAlignmentEnd()) {
            int elementsToAdd = read.getAlignmentEnd() - stopLocation;
            while (elementsToAdd-- > 0)
                countsWithBases.addLast(new CountWithBase(read.getAlignmentEnd() - elementsToAdd));
        }


        for (CigarElement cigarElement : cigar.getCigarElements()) {
            switch (cigarElement.getOperator()) {
                case H:
                case S:
                    // nothing to add to the window
                    break;
                case I:
                    // insertions are added to the base to the left (previous element) with the quality score of the first inserted base
                    countsWithBases.get(locationIndex - 1).addBase((byte) 'I', quals[readBaseIndex], true);
                    readBaseIndex += cigarElement.getLength();
                    break;
                case D:
                    // deletions are added to the counts with the read mapping quality as it's quality score
                    int nDeletionsToAdd = cigarElement.getLength();
                    while(nDeletionsToAdd-- > 0) {
                        countsWithBases.get(locationIndex).addBase((byte) 'D', (byte) read.getMappingQuality(), false);
                        locationIndex++;
                    }
                    break;
                case M:
                case P:
                case EQ:
                case X:
                    int nBasesToAdd = cigarElement.getLength();
                    while(nBasesToAdd-- > 0) {
                        countsWithBases.get(locationIndex).addBase(bases[readBaseIndex], quals[readBaseIndex], false);
                        readBaseIndex++;
                        locationIndex++;
                    }
                    break;
            }
        }
    }

//    // gets the region(s) of variable sites present in the window
//    private List<VariableRegion> getVariableRegions(int contextSize) {
//        // TODO require contextSize be non-Negative
//
//        List<VariableRegion> rawRegions = new LinkedList<VariableRegion>();
//
//        // Create every variant site
//        for (CountWithBase cBase: countsWithBases) {
//            if ( cBase.isVariant() ) {
//                rawRegions.add(createVariableRegion(cBase, contextSize));
//            }
//        }
//
//        if (rawRegions.isEmpty())
//            return rawRegions;
//
//        // Merge the proximate sites
//        List<VariableRegion> regions = new LinkedList<VariableRegion>();
//        Iterator<VariableRegion> i = rawRegions.iterator();
//        VariableRegion p = i.next();
//        while ( i.hasNext() ) {
//            VariableRegion q = i.next();
//            while ( p.end >= q.start ) {
//                p = p.merge(q);
//                if (i.hasNext() )
//                    q = i.next();
//                else
//                    break;
//            }
//            regions.add(p);
//            p = q;
//        }
//
//        if (regions.isEmpty())
//            regions.add(p);
//        //for ( VariableRegion vr: regions )
//            //System.out.println(String.format("INFO -- #### -- VARIANT FOUND Creating Variable Region, %d - %d", vr.start, vr.end));
//
//        return regions;
//    }
//
//    private VariableRegion createVariableRegion(CountWithBase cBase, int ContextSize) {
//        return new VariableRegion(cBase.location - ContextSize, cBase.location + ContextSize);
//    }
//
//    private List<SAMRecord> finalizeVariableRegion(VariableRegion variableRegion) {
//        //System.out.println(String.format(String.format("INFO ######################### \n \n Finalizing Variable Region: %d-%d ", variableRegion.start, variableRegion.end)));
//        List<SAMRecord> output = new LinkedList<SAMRecord>();
//
//        for ( SlidingRead read: SlidingReads ) {
//            SAMRecord SAM = read.trimToVariableRegion(variableRegion);
//            SAM.setReadName(SAM.getReadName()+".trim");
//            if ( SAM.getReadLength() > 0 ) {
//                output.add(SAM);
//                //System.out.println(String.format("Output Variable Read: %d-%d", SAM.getAlignmentStart(), SAM.getAlignmentEnd()));
//            }
//        }
//        slide(variableRegion.end+1);
//        createRunningConsensus();
//        return output;
//    }
//
//    private void createRunningConsensus() {
//
//        //header stuff
//        runningConsensus = new SAMRecord(header);
//        //runningConsensus.setAttribute("RG", reducedReadGroup.getId());
//        //runningConsensus.setAttribute(ReadUtils.REDUCED_READ_QUALITY_TAG, Integer.valueOf(REDUCED_READ_BASE_QUALITY)); //Also problems: Should it be avg
//        //runningConsensus.setAttribute("QE", QualityEquivalent);    // Qual equivs
//        runningConsensus.setReferenceName(contig);
//        runningConsensus.setReferenceIndex(contigIndex);
//        //runningConsensus.setReadName(String.format("%s.read.%d", reducedReadGroup.getId(), consensusCounter++));
//        runningConsensus.setReadPairedFlag(false);
//        runningConsensus.setReadUnmappedFlag(false);
//        //runningConsensus.setCigarString(String.format("%dM", size));    // TODO fix cigar
//        runningConsensus.setAlignmentStart(startLocation);
//
//    }
//
//    private List<SAMRecord> finalizeConsensusRead() {
//        List<SAMRecord> result = new LinkedList<SAMRecord>();
//        SAMRecord consensus = runningConsensus;
//        // This determines the end of read
//        if ( runningConsensus != null ){
//            if ( consensus.getReadLength() != 0 ) {
//                consensus.setMappingQuality(60); // TODO set up rms
//                consensus.setCigarString(String.format("%dM", consensus.getReadLength()));
//                consensus.setReadName(String.format("%s.read.%d", reducedReadGroup.getId(), consensusCounter++));
//                // Add only if no errors are present
//                if ( consensus.getReadLength() > 0)
//                    result.add(consensus);
//            }
//        }
//        runningConsensus = null;
//
//        return result;
//    }
//
//    private List<SAMRecord> finalizeConsensusRead(VariableRegion variableRegion) {
//        List<SAMRecord> result = new LinkedList<SAMRecord>();
//        addToConsensus(variableRegion.start);
//        SAMRecord consensus = runningConsensus;
//        // This determines the end of read
//        if ( runningConsensus != null ){
//            if ( consensus.getReadLength() != 0 ) {
//                consensus.setMappingQuality(60); // TODO set up rms
//                consensus.setCigarString(String.format("%dM", consensus.getReadLength()));
//                consensus.setReadName(String.format("%s.read.%d", reducedReadGroup.getId(), consensusCounter++));
//                result.add(consensus);
//            }
//        }
//        runningConsensus = null;
//
//        return result;
//    }
//
//    /**
//     * compresses the sliding reads to a streaming(incomplete) SamRecord
//     */
//    private List<SAMRecord> addToConsensus(int position) {
//
//
//        List<SAMRecord> result = new LinkedList<SAMRecord>();
//        if ( position > startLocation ) {
//            if ( !countsWithBases.isEmpty() ) {
//
//                if ( runningConsensus == null ) {
//                    createRunningConsensus();
//                }
//                if ( runningConsensus.getReadLength() == 0 ) {
//                    createRunningConsensus();
//                }
//                SAMRecord consensus;
//                /*
//                if ( position == -1 )
//                    position = (getEnd() + 1);
//                */
//
//                int size = position - startLocation;
//
//                //System.out.println(String.format("INFO -- Compressing running Consensus from %d to %d  ", getStart(), position));
//                int i = 0;
//                byte[] newBases = new byte[size];
//                byte[] newQuals = new byte[size];
//
//                int currentPosition = startLocation;
//                Iterator<SlidingWindow.CountWithBase> I = countsWithBases.iterator();
//                SlidingWindow.CountWithBase cBase = I.next();
//                while ( currentPosition < position ) {
//                    BaseCounts Counts = cBase.counts;
//                    if ( Counts.totalCount() == 0 ) {
//
//                        int total = Counts.totalCount();
//                        // This counts how many null bases are present
//                        int gap = 0;
//                        while ( total == 0) {
//                            gap++;
//                            if (I.hasNext()) {
//                                cBase = I.next();
//                                Counts = cBase.counts;
//                                total = Counts.totalCount();
//                            }
//                            else
//                                break;
//                        }
//
//                        result.addAll(finalizeConsensusRead());
//                        slide(cBase.location + 1);
//                        createRunningConsensus();
//                        break;
//                        // The recursive call acts as a retry to make the second consensus
//                        //result.addAll(addToConsensus(position));
//                    }
//                    else {
//                        newBases[i] = Counts.baseWithMostCounts();
//                        newQuals[i] = QualityUtils.boundQual(Counts.countOfMostCommonBase(), (byte) 64);
//                        //cBase = I.next();
//                        i++;
//                    }
//
//                    if (I.hasNext()){
//                        cBase = I.next();
//                        currentPosition = cBase.location;
//                    }
//                    else
//                        break;
//                }
//                byte[] tempBases = new byte[i];
//                byte[] tempQuals = new byte[i];
//                System.arraycopy(newBases, 0, tempBases, 0, i);
//                System.arraycopy(newQuals, 0, tempQuals, 0, i);
//                newBases = tempBases;
//                newQuals = tempQuals;
//
//            if ( runningConsensus.getReadLength() == 0 ) {
//                //createRunningConsensus();
//                //header stuff
//                consensus = runningConsensus;
//
//                consensus.setReadBases(newBases);
//                consensus.setBaseQualities(newQuals);
//            }
//            else {
//                consensus = runningConsensus;
//                //size += consensus.getReadLength();
//
//                byte[] bases = consensus.getReadBases();
//                byte[] quals = consensus.getBaseQualities();
//
//                bases = ArrayUtils.addAll(bases, newBases);
//                quals = ArrayUtils.addAll(quals, newQuals);
//
//                consensus.setReadBases(bases);
//                consensus.setBaseQualities(quals);
//            }
//            // TODO MQ RMS calc
//            // Once the new consensus is generated, we need to slide the window
//
//            slide(position);
//
//            runningConsensus = consensus;
//
//            }
//        }
//
//        return result;
//    }
//

}

