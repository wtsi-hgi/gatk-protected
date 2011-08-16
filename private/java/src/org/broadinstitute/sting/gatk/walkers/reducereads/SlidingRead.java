package org.broadinstitute.sting.gatk.walkers.reducereads;

import net.sf.samtools.*;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 8/3/11
 * Time: 2:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class SlidingRead {


    public SlidingRead(SAMRecord read) {
        alignmentStart = read.getAlignmentStart();
        this.read = read;
        // todo We should have a SAMRecordMetadata class
        // TODO add more field for cigar implementation, (unaligned start and end)
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();


        for (int i = 0; i < read.getReadBases().length; i++) {
            byte base = bases[i];
            byte qual = quals[i];

            BasesAndQuals.add(new BaseAndQual(base, qual));
        }
    }

    public LinkedList<BaseAndQual> getBasesAndQuals() {
        return BasesAndQuals;
    }

    public void setBasesAndQuals(LinkedList<BaseAndQual> basesAndQuals) {
        BasesAndQuals = basesAndQuals;
    }

    public int getAlignmentStart() {
        return alignmentStart;
    }

    public void setAlignmentStart(int alignmentStart) {
        this.alignmentStart = alignmentStart;
    }

    public SlidingRead trimToVariableRegion(VariableRegion variableRegion) {
        SlidingRead read = new SlidingRead(this.read);

        int start = getAlignmentStart();
        int stop = getAlignmentStop();

        // check to see if read is contained in region
        if ( start < variableRegion.end || stop > variableRegion.start ) {
            if ( start < variableRegion.start )
                read = read.clipStart(variableRegion.start);
            if ( stop > variableRegion.end )
                read = read.clipEnd(variableRegion.end);
        }
        return read;
    }
    //makes position the last element in LL
    private SlidingRead clipEnd(int position) {
        // Like hard clip but for sliding reads
        if ( position < getAlignmentStop()) {
            int toRemove = getAlignmentStop() - position ;
            // TODO Cigar handling here
            while ( toRemove > 0 && !BasesAndQuals.isEmpty() ) {
                BasesAndQuals.removeLast();
                toRemove--;
            }
        }
        return this;
    }

    public int getAlignmentStop() {
        return (getAlignmentStart() + BasesAndQuals.size() - 1 );
    }

    public SAMRecord toSAMRecord() {
        if (this == null)
            return null;
        SAMRecord output = read;
        output.setReadBases(getBaseArray());
        output.setBaseQualities(getQualArray());
        output.setCigarString(String.format("%dM", getQualArray().length)); // TODO fix cigar string
        output.setAlignmentStart(getAlignmentStart());

        return output;  //To change body of created methods use File | Settings | File Templates.
    }

    public byte[] getBaseArray() {
        byte[] output = new byte[BasesAndQuals.size()];
        Iterator<BaseAndQual> I = BasesAndQuals.iterator();

        for (int i = 0; i < output.length; i++ ) {
            BaseAndQual bQual = I.next();
            output[i] = bQual.base;
        }
        return output;
    }

   public byte[] getQualArray() {
        byte[] output = new byte[BasesAndQuals.size()];
        Iterator<BaseAndQual> I = BasesAndQuals.iterator();

        for (int i = 0; i < output.length; i++ ) {
            BaseAndQual bQual = I.next();
            output[i] = bQual.qual;
        }
        return output;
    }

    private class BaseAndQual {
        public byte base;
        public byte qual;

        public BaseAndQual( byte Base, byte Qual) {
            base = Base;
            qual = Qual;
        }

    }

    private LinkedList<BaseAndQual> BasesAndQuals = new LinkedList<BaseAndQual>();
    private int alignmentStart;
    private SAMRecord read;




    public SlidingRead clipStart(int position) {
        // Like hard clip but for sliding reads
        // position becomes start of the new read

        if ( position > this.alignmentStart) {
            int toRemove = position - this.alignmentStart;
            // TODO Cigar handling here
            while ( toRemove > 0  && !BasesAndQuals.isEmpty() ) {
                this.BasesAndQuals.pop();
                toRemove--;
            }
            setAlignmentStart(position);
        }
        return this;
    }


}

