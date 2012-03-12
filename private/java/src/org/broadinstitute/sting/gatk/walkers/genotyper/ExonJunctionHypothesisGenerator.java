package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 11/16/11
 * Time: 2:06 PM
 * To change this template use File | Settings | File Templates.
 */
@ReadFilters({DuplicateReadFilter.class,FailsVendorQualityCheckFilter.class,MappingQualityZeroFilter.class})
public class ExonJunctionHypothesisGenerator extends ReadWalker<TreeSet<ExonJunctionHypothesisGenerator.IntronLossJunctions>,TreeSet<ExonJunctionHypothesisGenerator.IntronLossJunctions>> {

    /**
     * A raw, unfiltered, highly specific callset in VCF format.
     */
    @Output
    protected PrintStream out = null;

    @Input(shortName="r",fullName="refSeq",required=true,doc="The RefSeq Gene definition track")
    public RodBinding<RefSeqFeature> refSeqRodBinding;

    @Input(shortName="mis",fullName="minInsertSize",required=false,doc="The minimum insert size for a read pair to consider as evidence of joined exons")
    public int minInferredInsert = 600;


    private IndexedFastaSequenceFile referenceReader;

    private TreeMap<GenomeLoc,byte[]> exonSequences;

    public TreeSet<IntronLossJunctions> reduceInit() { return new TreeSet<IntronLossJunctions>(); }

    public class IntronLossJunctions implements Comparable{

        private Set<Pair<Integer,Integer>> junctions;
        private String uniqueName;
        private GenomeLoc geneLocus;

        private IntronLossJunctions() {
            junctions = new TreeSet<Pair<Integer,Integer>>(new Comparator<Pair<Integer, Integer>>() {
                @Override
                public int compare(Pair<Integer, Integer> integerIntegerPair, Pair<Integer, Integer> integerIntegerPair1) {
                    int first = integerIntegerPair.first - integerIntegerPair1.first;
                    int second = integerIntegerPair.second - integerIntegerPair1.second;
                    return first != 0 ? first : second;
                }
            });
        }

        public IntronLossJunctions(RefSeqFeature feature) {
            uniqueName = feature.getTranscriptUniqueGeneName();
            geneLocus = feature.getLocation();
            junctions = new TreeSet<Pair<Integer,Integer>>(new Comparator<Pair<Integer, Integer>>() {
                @Override
                public int compare(Pair<Integer, Integer> integerIntegerPair, Pair<Integer, Integer> integerIntegerPair1) {
                    int first = integerIntegerPair.first - integerIntegerPair1.first;
                    int second = integerIntegerPair.second - integerIntegerPair1.second;
                    return first != 0 ? first : second;
                }
            });
        }

        public void addJunction(Pair<Integer,Integer> junction) {
            junctions.add(junction);
        }

        public int compareTo(Object other) {
            if ( other instanceof IntronLossJunctions ) {
                int locComp = this.geneLocus.compareTo(((IntronLossJunctions)other).geneLocus);
                int nameComp = this.uniqueName.compareTo(((IntronLossJunctions)other).uniqueName);
                if ( locComp == 0 && nameComp == 0 ) {
                    return 0;
                } else if ( locComp == 0 ) {
                    return nameComp;
                } else {
                    return locComp;
                }
            } else {
                return Integer.MIN_VALUE;
            }
        }

        public void add(IntronLossJunctions other) {
            if ( this.compareTo(other) != 0 ) {
                throw new ReviewedStingException("Attempting to merge junctions across different genes or gene transcripts");
            }
            this.junctions.addAll(other.junctions);
        }

        public String toString() {
            return String.format("%s\t%s\t%s\t%s", geneLocus, uniqueName, joinJunctions(), reconcile().joinJunctions());
        }

        private String joinJunctions() {
            StringBuilder buf = new StringBuilder();
            for ( Pair<Integer,Integer> j : junctions ) {
                buf.append(j.first);
                buf.append(',');
                buf.append(j.second);
                buf.append(';');
            }
            if ( buf.length() > 0 && buf.charAt(buf.length()-1) == ';') {
                buf.deleteCharAt(buf.lastIndexOf(";"));
            }
            return buf.toString();
        }

        private IntronLossJunctions reconcile() {
            IntronLossJunctions reconciled = new IntronLossJunctions();
            reconciled.geneLocus = this.geneLocus;
            reconciled.uniqueName = this.uniqueName;
            // build up an adjacency map, where each entry contains a list of
            // all entries it connects to (before and after)
            HashMap<Integer,TreeSet<Integer>> adjMap = new HashMap<Integer,TreeSet<Integer>>();
            for ( Pair<Integer,Integer> j : this.junctions ) {
                if ( ! adjMap.containsKey(j.first) ) {
                    adjMap.put(j.first,new TreeSet<Integer>());
                }
                adjMap.get(j.first).add(j.second);
                if ( ! adjMap.containsKey(j.second) ) {
                    adjMap.put(j.second, new TreeSet<Integer>());
                }
                adjMap.get(j.second).add(j.first);
            }

            Integer prevStart = null;
            Integer prevStop = null;
            HashSet<Integer> usedPositions = new HashSet<Integer>();
            for ( Pair<Integer,Integer> etry : this.junctions ) {
                Integer potentialStart = etry.first;
                // rule 1: always take the shortest (e.g. first) connection
                boolean isFirstConnection = !usedPositions.contains(etry.first);
                Integer potentialStop = etry.second;
                // rule 2: check to see if the stop position has a shorter (prior) connection
                TreeSet<Integer> stopConnections = adjMap.get(potentialStop);
                boolean targetHasShorter = false;
                for ( Integer p : stopConnections ) {
                    if ( p >= potentialStop ) {
                        break;
                    }
                    targetHasShorter |= ( p > potentialStart && p < potentialStop );
                }
                // rule 3: check to see if the junction has subordinated junctions (hooray n^2)
                boolean targetHasSubordinated = false;
                for ( int i = potentialStart+1; i < potentialStop; i++ ) {
                    if ( adjMap.containsKey(new Integer(i)) ){
                        if ( adjMap.get(new Integer(i)).first() <= potentialStop ) {
                            targetHasSubordinated = true;
                            break;
                        }
                    }
                }
                // if there are no shorter connections for the target and no subordinated junctions
                if ( ! targetHasShorter && ! targetHasSubordinated && isFirstConnection ) {
                    // do we need to fill anything in?
                    if ( prevStart == null || potentialStart.equals(prevStop) ) {
                        // no
                        reconciled.addJunction(new Pair<Integer,Integer>(potentialStart,potentialStop));
                        usedPositions.add(potentialStart);
                    } else {
                        if ( prevStop < potentialStart ) {
                            for ( int j = prevStop; j < potentialStart; j++ ) {
                                reconciled.addJunction(new Pair<Integer,Integer>(j,j+1));
                                usedPositions.add(j);
                            }
                        } else {
                            for ( int j = prevStart; j < potentialStart; j++ ) {
                                reconciled.addJunction(new Pair<Integer,Integer>(j,j+1));
                                usedPositions.add(j);
                            }
                        }
                        reconciled.addJunction(new Pair<Integer,Integer>(potentialStart,potentialStop));
                        usedPositions.add(potentialStart);
                    }

                    prevStart = potentialStart;
                    prevStop = potentialStop;
                }
            }

            return reconciled;
        }
    }

    public TreeSet<IntronLossJunctions> map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {
        List<RefSeqFeature> refSeqFeatures = new ArrayList<RefSeqFeature>(16);
        for (GATKFeature feature : metaDataTracker.getAllCoveringRods() ) {
            if ( feature.getUnderlyingObject().getClass().isAssignableFrom(RefSeqFeature.class) ) {
                refSeqFeatures.add((RefSeqFeature) feature.getUnderlyingObject());
            }
        }

        TreeSet<IntronLossJunctions> junctionsSet = new TreeSet<IntronLossJunctions>();
        GenomeLoc readLoc = getToolkit().getGenomeLocParser().createGenomeLoc(read);
        GenomeLoc mateLoc = ( (! read.getReadPairedFlag()) || read.getMateUnmappedFlag() ) ? null : getToolkit().getGenomeLocParser().createGenomeLoc(read.getMateReferenceName(),read.getMateAlignmentStart(),read.getMateAlignmentStart()+read.getReadLength());
        for ( RefSeqFeature refSeqFeature : refSeqFeatures ) {
            int readOverlapInt = refSeqFeature.getSortedOverlapInteger(readLoc);
            if ( mateLoc != null && readOverlapInt > -1 ) {
                int inferredDistance = Math.abs(read.getInferredInsertSize());
                int mateOverlapInt = refSeqFeature.getSortedOverlapInteger(mateLoc);
                Integer mateMapQ = (Integer) read.getAttribute("MQ");
                if ( mateMapQ == null ) {
                    throw new UserException.MalformedFile("Paired reads in the BAM file must have a mate mapping quality tag (MQ) filled in.");
                }
                if ( mateOverlapInt > -1 && mateOverlapInt != readOverlapInt && inferredDistance >= minInferredInsert && mateMapQ > 23 ) {
                    IntronLossJunctions j = new IntronLossJunctions(refSeqFeature);
                    if ( readOverlapInt < mateOverlapInt ) {
                        j.addJunction(new Pair<Integer,Integer>(readOverlapInt,mateOverlapInt));
                    } else {
                        j.addJunction(new Pair<Integer,Integer>(mateOverlapInt,readOverlapInt));
                    }
                    junctionsSet.add(j);
                }
            }

            int unclStart = read.getUnclippedStart();
            int unclEnd = read.getUnclippedEnd();
            if ( readOverlapInt > -1 && readLoc.getStop() > refSeqFeature.getSortedExonLoc(readOverlapInt).getStop() ) {
                for ( int k = readOverlapInt+1; k < refSeqFeature.getNumExons() ; k++ ) {
                    if ( readMatches(read,readLoc,refSeqFeature.getSortedExonLoc(k), unclStart, unclEnd) ) {
                        IntronLossJunctions j = new IntronLossJunctions(refSeqFeature);
                        j.addJunction(new Pair<Integer,Integer>(readOverlapInt,k));
                        junctionsSet.add(j);
                        break;
                    }
                }
            } else if ( readOverlapInt > -1 && readLoc.getStart() < refSeqFeature.getSortedExonLoc(readOverlapInt).getStart() ) {
                for ( int k = readOverlapInt-1; k > -1; k--) {
                    if ( readMatches(read,readLoc,refSeqFeature.getSortedExonLoc(k),unclStart,unclEnd) ) {
                        IntronLossJunctions j = new IntronLossJunctions(refSeqFeature);
                        j.addJunction(new Pair<Integer,Integer>(k,readOverlapInt));
                        junctionsSet.add(j);
                        break;
                    }
                }
            }
        }

        if ( exonSequences != null && exonSequences.size() > 50 ) {
            TreeMap<GenomeLoc,byte[]> keep = new TreeMap<GenomeLoc, byte[]>();
            for ( GenomeLoc l : exonSequences.keySet() ) {
                if ( ! readLoc.isPast(l) || readLoc.distance(l) < 10000 ) {
                    keep.put(l,exonSequences.get(l));
                }
            }
            exonSequences = keep;
        }

        return junctionsSet;
    }


    private boolean readMatches(GATKSAMRecord read, GenomeLoc readLoc, GenomeLoc exonPos, int readUnStart, int readUnEnd) {
        int match = 0;
        if ( readLoc.isPast(exonPos) ) {
            int nBases = read.getAlignmentStart()-readUnStart;
            if (nBases < 10) {
                return false;
            }
            // first bases of read should match the previous exon
            byte[] refBases = getExonBases(exonPos, nBases, false);
            //referenceReader.getSubsequenceAt(exonPos.getContig(),exonPos.getStop()-nBases,exonPos.getStop()).getBases();
            byte[] readBases = read.getReadBases();
            int max = Math.min(refBases.length,10);
            for ( int i = 0; i < max; i++ ) {
                match += (refBases[i] == readBases[i]) ? 1 : 0;
            }
        } else {
            // last bases of read should match the next exon
            int nBases = readUnEnd-read.getAlignmentEnd();
            if ( nBases < 10 ) {
                return false;
            }
            byte[] refBases = getExonBases(exonPos, nBases, true);
            //referenceReader.getSubsequenceAt(exonPos.getContig(),exonPos.getStart(),exonPos.getStart()+nBases).getBases();
            byte[] readBases = read.getReadBases();
            ArrayUtils.reverse(readBases);
            for ( int i = 0; i < 10; i++) {
                match += (refBases[i] == readBases[i]) ? 1 : 0;
            }
        }

        return match >= 8;
    }


    private byte[] getExonBases(GenomeLoc exonPosition, int numBases, boolean startOfExon) {
        if (this.exonSequences == null || ! this.exonSequences.containsKey(exonPosition) ) {
            if ( this.exonSequences == null ) {
                exonSequences = new TreeMap<GenomeLoc,byte[]>();
            }
            exonSequences.put(exonPosition,referenceReader.getSubsequenceAt(exonPosition.getContig(),
                    exonPosition.getStart(),exonPosition.getStop()).getBases());
        }

        byte[] sq = exonSequences.get(exonPosition);
        //logger.debug(String.format("%d: %d %d",sq.length,sq.length,numBases));
        return startOfExon ? Arrays.copyOfRange(sq,0,numBases) : Arrays.copyOfRange(sq,Math.max(0,sq.length-numBases),sq.length);
    }

    public TreeSet<IntronLossJunctions> reduce(TreeSet<IntronLossJunctions> pMap, TreeSet<IntronLossJunctions> pRed) {
        for ( IntronLossJunctions ilj : pMap ) {
            if ( pRed.contains(ilj) ) {
                pRed.tailSet(ilj).first().add(ilj);
            } else {
                pRed.add(ilj);
            }
        }

        return pRed;
    }

    public void onTraversalDone(TreeSet<IntronLossJunctions> junctions) {
        out.printf("HEADER%s\t%s\t%s\t%s%n","Loc","Transcript","Junctions","Potential Event");
        for ( IntronLossJunctions ilj : junctions ) {
            out.printf("%s%n",ilj.toString());
        }
    }

    public void initialize() {
        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile,ex);
        }
    }
}
