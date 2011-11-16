package org.broadinstitute.sting.gatk.walkers.newassociation.hierarchicalbinning;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
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
public class PairedReadBin implements HasGenomeLocation, Iterable<Pair<GATKSAMRecord,GATKSAMRecord>> {

    private final HashMap<String,Pair<GATKSAMRecord,GATKSAMRecord>> reads = new HashMap<String, Pair<GATKSAMRecord, GATKSAMRecord>>(3200);
    private GenomeLoc loc = null;
    private GenomeLocParser parser;

    public PairedReadBin(GenomeLocParser p) { parser = p; }

    // Return false if we can't process this read bin because the reads are not correctly overlapping.
    // This can happen if e.g. there's a large known indel with no overlapping reads.
    public void add(GATKSAMRecord read) {
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
            Pair<GATKSAMRecord,GATKSAMRecord> samPair = new Pair<GATKSAMRecord,GATKSAMRecord>(read,null);
            reads.put(name,samPair);
        }
    }

    public Collection<Pair<GATKSAMRecord,GATKSAMRecord>> getReadPairs() { return reads.values(); }

    public GenomeLoc getLocation() { return loc; }

    public int size() { return reads.size(); }

    public void clear() {
        reads.clear();
        loc = null;
    }

    public Iterator<Pair<GATKSAMRecord,GATKSAMRecord>> iterator() {
        return  reads.values().iterator();
    }
}
