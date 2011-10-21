package org.broadinstitute.sting.gatk.walkers.newassociation.hierarchicalbinning;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 9/28/11
 * Time: 6:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class PairedReadBin implements HasGenomeLocation, Iterable<Pair<SAMRecord,SAMRecord>> {

    private final HashMap<String,Pair<SAMRecord,SAMRecord>> reads = new HashMap<String, Pair<SAMRecord, SAMRecord>>(3200);
    private GenomeLoc loc = null;
    private GenomeLocParser parser;

    public PairedReadBin(GenomeLocParser p) { parser = p; }

    // Return false if we can't process this read bin because the reads are not correctly overlapping.
    // This can happen if e.g. there's a large known indel with no overlapping reads.
    public void add(SAMRecord read) {
        if ( read.getCigarString().equals("*") ) {
            return;
        }
        GenomeLoc locForRead = parser.createGenomeLoc(read);
        if ( loc == null )
            loc = locForRead;
        else if ( locForRead.getStop() > loc.getStop() )
            loc = parser.createGenomeLoc(loc.getContig(), loc.getStart(), locForRead.getStop());

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

    public GenomeLoc getLocation() { return loc; }

    public int size() { return reads.size(); }

    public void clear() {
        reads.clear();
        loc = null;
    }

    public Iterator<Pair<SAMRecord,SAMRecord>> iterator() {
        return  reads.values().iterator();
    }
}
