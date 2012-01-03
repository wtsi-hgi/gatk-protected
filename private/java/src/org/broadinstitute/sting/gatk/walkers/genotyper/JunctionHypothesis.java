package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.sting.utils.codecs.table.TableFeature;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 11/16/11
 * Time: 9:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class JunctionHypothesis implements Comparable, HasGenomeLocation {

    List<Pair<Integer,Integer>> exonJunctionNumbers;
    List<Pair<GenomeLoc,GenomeLoc>> exonJunctionLocs;
    private Map<GenomeLoc,Integer> baseBeforeCache;
    List<GenomeLoc> exons;
    boolean isForward;
    String sequence;
    String name;

    public JunctionHypothesis(String transcript, RefSeqFeature feature, List<Pair<Integer,Integer>> junctions, IndexedFastaSequenceFile refReader) {
        if ( ! transcript.equals(feature.getTranscriptUniqueGeneName()) ) {
            throw new ReviewedStingException("Attempting to form a junction sequence for a RefSeq record whose feature name does not match");
        }
        exons = new ArrayList<GenomeLoc>(feature.getExons()); // ensures same ordering as hypothesis
        Collections.sort(exons);
        exonJunctionNumbers = junctions;
        exonJunctionLocs = new ArrayList<Pair<GenomeLoc,GenomeLoc>>(exonJunctionNumbers.size());
        Integer last = null;
        StringBuilder seqBuilder = new StringBuilder();
        exons = new ArrayList<GenomeLoc>(exonJunctionLocs.size());
        for ( Pair<Integer,Integer> junction : exonJunctionNumbers ) {
            GenomeLoc first = feature.getSortedExonLoc(junction.first);
            GenomeLoc second = feature.getSortedExonLoc(junction.second);
            exonJunctionLocs.add(new Pair<GenomeLoc,GenomeLoc>(first,second));
            if ( last == null || last.equals(junction.first) ) {
                String exonSeq = new String(refReader.getSubsequenceAt(first.getContig(),first.getStart(),first.getStop()).getBases());
                seqBuilder.append(exonSeq);
            } else {
                throw new ReviewedStingException("Attempting to form a junction sequence consisting of discontinous ending and starting exons");
            }
            last = junction.second;
        }
        for ( Pair<GenomeLoc,GenomeLoc> locPair : exonJunctionLocs ) {
            if ( exons.size() == 0 ) {
                exons.add(locPair.first);
            }
            exons.add(locPair.second);
        }
        GenomeLoc lastLoc = exonJunctionLocs.get(exonJunctionLocs.size()-1).second;
        String lastString = new String(refReader.getSubsequenceAt(lastLoc.getContig(),lastLoc.getStart(),lastLoc.getStop()).getBases());
        seqBuilder.append(lastString);
        sequence = seqBuilder.toString();
        isForward = feature.getStrand() == 1;
        name = transcript;
        baseBeforeCache = new HashMap<GenomeLoc,Integer>(exonJunctionLocs.size());
    }

    public String getSequence() {
        return sequence;
    }

    public int size() {
        return sequence.length();
    }

    public int getBaseOffset(GenomeLoc position) {
        // count up the number of bases in exons prior to this position (involved in the hypothesis)
        if ( baseBeforeCache.size() == 0 ) {
            // fill in the cache
            boolean first = true;
            for ( Pair<GenomeLoc,GenomeLoc> jun : exonJunctionLocs ) {
                if ( first ) {
                    baseBeforeCache.put(jun.first,0);
                    first = false;
                }
                baseBeforeCache.put(jun.second,baseBeforeCache.get(jun.first)+ (int)jun.first.size());
            }
        }

        if ( getExonLoc(position) == null ) {
            System.out.println("foo");
        }
            GenomeLoc exon = getExonLoc(position);
        int nBefore = baseBeforeCache.get(exon);

        return nBefore + position.getStart()-exon.getStart();
    }

    public GenomeLoc getLocation() {
        GenomeLoc exonStart = exons.get(0).getStartLocation();
        GenomeLoc exonStop = exons.get(exons.size()-1).getStopLocation();
        return exonStart.endpointSpan(exonStop);
    }

    public int compareTo(Object other) {
        if ( other instanceof JunctionHypothesis ) {
            int locComp = getLocation().compareTo(((JunctionHypothesis) other).getLocation());
            if ( locComp == 0 ) {
                return this.name.compareTo(((JunctionHypothesis) other).name);
            }
        }
        return Integer.MIN_VALUE;
    }

    private double bestUnmappedScore = Double.POSITIVE_INFINITY;
    private double bestMappedScore = Double.POSITIVE_INFINITY;

    public void updateBestScore(double score, boolean unmapped) {
        if ( unmapped ) {
            bestUnmappedScore = score < bestUnmappedScore ? score : bestUnmappedScore;
        } else {
            bestMappedScore = score < bestMappedScore ? score : bestMappedScore;
        }
    }

    public double getScoreDifference() {
        return bestUnmappedScore - bestMappedScore;
    }

    public int getInsertAdjustment(GenomeLoc read, GenomeLoc mate) {
        int dist = -1;
        GenomeLoc ren = getExonLoc(read);
        if ( ren == null ) {
            return dist;
        }
        GenomeLoc men = getExonLoc(mate);
        if ( men == null || men.equals(ren) ) {
            return dist;
        }
        dist = men.distance(ren);
        return dist;
    }

    private GenomeLoc getExonLoc(GenomeLoc loc) {
        for ( GenomeLoc eLoc : exons ) {
            if ( eLoc.overlapsP(loc)) {
                return eLoc;
            } else if ( eLoc.isPast(loc) ) {
                return null;
            }
        }

        return null;
    }

    public boolean overlapsExonIntron(GenomeLoc loc) {
        for ( GenomeLoc eLoc : exons ) {
            if ( eLoc.overlapsP(loc) ) {
                // see if it starts before eloc or ends after eloc
                if ( loc.startsBefore(eLoc) ) {
                    // make sure that eloc isn't hte first
                    if ( ! eLoc.equals(exons.get(0)) ) {
                        return true;
                    } else {
                        return false;
                    }
                } else if ( eLoc.startsBefore(loc) ) {
                    // make sure that eloc isn't hte last
                    if ( ! eLoc.equals(exons.get(exons.size()-1))) {
                        return true;
                    } else {
                        return false;
                    }
                }
            }
        }

        return false;
    }

    public String toString() {
        return String.format("<%s:%s>",name,joinExonNums());
    }

    private String joinExonNums() {
        StringBuilder builder = new StringBuilder();
        if ( exonJunctionNumbers.size() <= 0 ) {
            return "";
        }
        for ( Pair<Integer,Integer> e : exonJunctionNumbers ) {
            builder.append(e.first);
            builder.append('-');
        }
        builder.append(exonJunctionNumbers.get(exonJunctionNumbers.size()-1).second);
        return builder.toString();
    }

    public static List<JunctionHypothesis> generateHypotheses(List<TableFeature> rawHypos, Map<String,RefSeqFeature> refSeqFeaturesByGene, IndexedFastaSequenceFile reader) {
        List<JunctionHypothesis> hypotheses = new ArrayList<JunctionHypothesis>(rawHypos.size());
        for ( TableFeature rhy : rawHypos ) {
            RefSeqFeature feature = refSeqFeaturesByGene.get(rhy.getValue(1));
            Set<List<Pair<Integer,Integer>>> pairs = unwrapPairs(rhy.getValue(3));
            for ( List<Pair<Integer,Integer>> hypo : pairs ) {
                hypotheses.add(new JunctionHypothesis(rhy.getValue(2),feature,hypo,reader));
            }
        }

        return hypotheses;
    }

    public static Set<List<Pair<Integer,Integer>>> unwrapPairs(String rawString) {
        Set<List<Pair<Integer,Integer>>> pairSet = new HashSet<List<Pair<Integer,Integer>>>();
        for ( String hypo : rawString.split("\\|") ) {
            List<Pair<Integer,Integer>> pairs = new ArrayList<Pair<Integer,Integer>>(hypo.split(";").length);
            for ( String rawPair : hypo.split(";")) {
                String[] exons = rawPair.split(",");
                Integer e1 = Integer.parseInt(exons[0]);
                Integer e2 = Integer.parseInt(exons[1]);
                Pair<Integer,Integer> pair = new Pair<Integer,Integer>(e1,e2);
                pairs.add(pair);
            }

            pairSet.add(pairs);
        }

        return pairSet;
    }
}
