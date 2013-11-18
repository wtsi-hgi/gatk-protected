/*
*  By downloading the PROGRAM you agree to the following terms of use:
*
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
*
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/
package org.broadinstitute.sting.gatk.walkers.haplotypecaller.graphs;

import org.broadinstitute.sting.gatk.walkers.haplotypecaller.Kmer;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.haplotype.Haplotype;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: valentin
 * Date: 11/8/13
 * Time: 4:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class KmerSearchableSeqGraph extends SeqGraph {

    private final Map<Kmer,KmerMap> uniqueKmers;
    private final Map<Kmer,Set<KmerMap>> nonUniqueKmers;
    private final Map<SeqVertex,Pair<SeqVertex,Integer>> vertexReplacement;

    private boolean needToCreateSearchTables;

    public KmerSearchableSeqGraph(final int kmerSize) {
        super(kmerSize);
        uniqueKmers = new HashMap<>(100);
        nonUniqueKmers = new HashMap<>(100);
        vertexReplacement = new HashMap<>(50);
    }

    @Override
    public boolean addVertex(SeqVertex v) {
        needToCreateSearchTables = true;
        return super.addVertex(v);
    }

    @Override
    public void addVertices(SeqVertex ... v) {
        needToCreateSearchTables = true;
        super.addVertices(v);
    }

    @Override
    public BaseEdge addEdge(SeqVertex v1, SeqVertex v2) {
        needToCreateSearchTables = true;
        return super.addEdge(v1,v2);
    }

    /**
     * Update the kmer -> vertex search maps if needed.
     */
    private void createSearchTablesIfNeeded() {
        if (!needToCreateSearchTables) {
            return;
        }
        uniqueKmers.clear();
        nonUniqueKmers.clear();
        final Map<SeqVertex,List<Set<KmerMap>>> doneSet = new HashMap<>(this.vertexSet().size());
        for (final SeqVertex sv : this.vertexSet()) {
            calculateKmerMap(sv,doneSet);
        }
        needToCreateSearchTables = false;
    }


    /**
     * Given a SeqVertex return all the possible kmer maps that overlap it sequences with a end base within its sequence.
     * @param sv the target SeqVertex
     * @param done reusable results cache.
     * @return a List of sets where each set contain possible kmer maps finishing at that offset in the SeqVertex.
     */
    private List<Set<KmerMap>> calculateKmerMap(final SeqVertex sv, final Map<SeqVertex, List<Set<KmerMap>>> done) {
        if (done.containsKey(sv)) {
            return done.get(sv);
        } else {
            final List<Set<KmerMap>> kmerMapList;
            if (isSource(sv)) {
                final byte[] sequence = sv.sequence;
                final int count = sequence.length - kmerSize + 1;
                kmerMapList = new ArrayList<Set<KmerMap>>(count);
                for (int i = 0; i < count; i++) {
                    kmerMapList.add(Collections.singleton(new KmerMap(this,new Kmer(sequence, i, kmerSize), sv, i, sv, i + kmerSize)));
                }
            } else {
                final Set<SeqVertex> parents = (Set<SeqVertex>) (Set) this.incomingVerticesOf(sv);
                final byte[] sequence = sv.sequence;
                final int count = sequence.length;
                kmerMapList = new ArrayList<>(count);
                for (int i = 0; i < count; i++) {
                    kmerMapList.add(new HashSet<KmerMap>(parents.size() * 5));
                }

                for (final SeqVertex parent : parents) {
                    final List<Set<KmerMap>> parentKmerMapList = calculateKmerMap(parent, done);
                    Set<KmerMap> lastKmerMapSet = parentKmerMapList.get(parentKmerMapList.size() - 1);
                    for (int i = 0; i < count; i++) {
                        final byte nextChar = sequence[i];
                        final Set<KmerMap> currentKmerMapSet = kmerMapList.get(i);
                        for (final KmerMap km : lastKmerMapSet) {
                            final KmerMap nextKmerMap = shiftKmerMap(km, nextChar);
                            currentKmerMapSet.add(nextKmerMap);

                        }
                        lastKmerMapSet = currentKmerMapSet;
                    }
                }
            }
            done.put(sv, kmerMapList);
            addKmerMaps(kmerMapList);
            return kmerMapList;
        }
    }

    /**
     * Produces the kmer map results of moving along the graph (sliding kmer) by one base.
     * @param km
     * @param nextChar
     * @return
     */
    private KmerMap shiftKmerMap(final KmerMap km, final byte nextChar) {
        SeqVertex startVertex = km.getStartVertex();
        int startOffset = km.getStartOffset();
        if (startOffset == startVertex.length() - 1) {
            startVertex = outgoingVertexOf(startVertex, km.kmer.base(1));
            startOffset = 0;
        }
        if (startVertex == null) {
            throw new StingException("");
        }
        SeqVertex endVertex = km.getEndVertex();
        int endOffset = km.getEndOffset();
        if (endOffset == endVertex.length()) {
            endVertex = outgoingVertexOf(startVertex, nextChar);
            endOffset = 1;
        }
        if (endVertex == null) {
            throw new StingException("");
        }
        final Kmer newKmer = km.kmer.shift(nextChar);
        return new KmerMap(this,newKmer, startVertex, startOffset, endVertex, endOffset);
    }


    /**
     * Return the SeqVertex that follows this one based on the next base (kmer-sliding).
     * @param startVertex
     * @param nextChar
     * @return
     */
    private SeqVertex outgoingVertexOf(final SeqVertex startVertex, final byte nextChar) {
        for (final SeqVertex nextStart : this.outgoingVerticesOf(startVertex)) {
            if (nextStart.sequence[0] == nextChar) {
                return nextStart;
            }
        }
        return null;
    }


    /**
     * Helping method when updating the kmer search maps, it adds the set of kmer maps into the unique and non-unique kmers table
     * depending.
     *
     * @param kmerMapList
     */
    private void addKmerMaps(final List<Set<KmerMap>> kmerMapList) {
        for (final Set<KmerMap> kms : kmerMapList) {
            for (final KmerMap km : kms) {
                if (uniqueKmers.containsKey(km.kmer)) {
                    final KmerMap otherKm = uniqueKmers.get(km.kmer);
                    if (!otherKm.equals(km)) {
                        final Set<KmerMap> maps = new HashSet<>(5);
                        nonUniqueKmers.put(km.kmer,maps);
                        maps.add(km);
                        uniqueKmers.remove(km.kmer);
                    }
                } else if (nonUniqueKmers.containsKey(km.kmer)) {
                    nonUniqueKmers.get(km.kmer).add(km);
                } else {
                    uniqueKmers.put(km.kmer,km);
                }
            }
        }
    }

    /**
     * Look up the mapping position of a kmer given it is unique.
     * @param kmer the target kmer.
     * @return never {@code null} but potentially unmapped kmerMap.
     */
    public KmerMap findUniqueKmerMap(final Kmer kmer) {
        createSearchTablesIfNeeded();
        KmerMap result = uniqueKmers.get(kmer);
        if (result == null) {
            result = KmerMap.unmapped(this,kmer);
        }
        return result;
    }

    /**
     * Find all unique kmer maps for a give sequence.
     *
     * <p>
     *     There will be one map in the result for each kmer in the sequence.
     * </p>
     * @param sequence target sequence.
     * @return never {@code null}, a list with as many elements as kmer in the input sequence.
     */
    public List<KmerMap> findKmerMaps(final byte[] sequence) {
        final int ks = this.kmerSize;
        final int count = sequence.length - ks + 1;
        final List<KmerMap> result = new ArrayList<>(count);
        for (int i = 0; i < count; i++) {
            final Kmer kmer = new Kmer(sequence,i,ks);
            final KmerMap kmerMap = findUniqueKmerMap(kmer);
            result.add(kmerMap);
        }
        return result;
    }

    /**
     * Find all unique kmer maps for all kmers in a haplotype.
     * @param h target haplotype
     * @return never {@code null} but a list with with as many elements are kmers in the haplotype.
     */
    public List<KmerMap> findKmerMaps(final Haplotype h) {
        return findKmerMaps(h.getBases());
    }

    /**
     * Find all unique kmer maps for all kmers in a read
     *
     * @param r target read.
     * @return never {@code null } but a list with as many elements are kmers in the read.
     */
    public List<KmerMap> findKmerMaps(final GATKSAMRecord r) {
        return findKmerMaps(r.getReadBases());
    }


    /**
     * Reference to the vertex sorted for this graph.
     */
    private final VertexSorter<SeqVertex,BaseEdge> vertexSorter = new VertexSorter(this);

    /**
     * Reference to the kmer map comparator based on the sorter {@link #vertexSorter}
     */
    private final Comparator<KmerMap> kmerMapComparator = new Comparator<KmerMap>() {
        @Override
        public int compare(final KmerMap o1, final KmerMap o2) {
            final SeqVertex v1 = o1.getEndVertex();
            final SeqVertex v2 = o2.getEndVertex();

            VertexOrder vo = vertexSorter.vertexOrder(o1.getEndVertex(), o2.getEndVertex());
            switch (vo) {
                case SAME:
                    int o1offset = o1.getEndOffset();
                    int o2offset = o2.getEndOffset();
                    return (o1offset == o2offset) ? 0 : o1offset < o2offset ? -1 : 1;
                case BEFORE: return -1;
                case AFTER: return 1;
                case PARALLEL:
                    return (v1.getId() < v2.getId()) ? -1 : 1;
                default:
                    throw new RuntimeException("unexpected vertex-order: " + vo);
            }
        }
    };

    /**
     * Returns any replacement registered on a vertex.
     * @param vertex target vertex.
     * @return never {@code return}.
     */
    public Pair<SeqVertex, Integer> vertexReplacement(final SeqVertex vertex) {
        return vertexReplacement.get(vertex);
    }


    private class SequenceMap {
        protected final byte[] sequence;
        protected final List<KmerMap> kmerMaps;
        protected final int[] order;
        protected final int length;
        protected int sortCost = 0;

        private SequenceMap(final byte[] sequence) {
            this.sequence = sequence;
            kmerMaps = findKmerMaps(sequence);
            length = sequence.length;
            order = calculateOrder(kmerMaps);
        }

        /**
         * Calculates an array indicating the order of each kmer maps onto the array.
         *
         * This method uses bubble sort with the conviction that in most cases the list
         * is already sorted.
         *
         *
         * @param kmerMaps list of kmerMaps.
         * @return
         */
        private int[] calculateOrder(final List<KmerMap> kmerMaps) {
            final int[] result = new int[length];
            // First we move any unmapped kmerMaps to the end and mapped ones too the beginning of
            // result index array.
            int nextMapped = 0;
            int nextUnmapped = length - 1;
            for (int i = 0; i < length; i++) {
                final KmerMap km = kmerMaps.get(i);
                if (km.isUnmapped()) {
                    result[nextUnmapped--] = i;
                } else {
                    result[nextMapped++] = i;
                }
            }
            // nextMapped points the first unmapped at this point, and is the effective size of the
            // mapped kmers.
            // Now be proceed bubble-sorting.
            // The worst scenario log10MLE is O(n^2), how ever it is much faster if we know that the
            // array is mostly sorted.

            for (int next = 1; next < nextMapped; next++) {
                final KmerMap nextKm = kmerMaps.get(result[next]);
                int prev = next - 1;
                final KmerMap prevKm = kmerMaps.get(result[prev]);
                while (prev > 0) {
                    if (kmerMapComparator.compare(nextKm,prevKm) >= 0)
                        continue;
                    result[prev + 1] = result[prev];
                    sortCost++;
                    prev--;
                }
                result[prev + 1] = next;
            }
            return result;

        }

        private int sortCost() {
            return sortCost;
        }

        private int[][] alignTo(final SequenceMap h) {
            final int maxLength = this.length < h.length ? length : h.length;
            final int[][] result = new int[2][];
            result[0] = new int[maxLength];
            result[1] = new int[maxLength];

            int j = 0;
            int next = 0;
            for (int i = 0; i < h.length; i++) {
                final KmerMap hKm = h.kmerMaps.get(i);
                while (j < length) {
                    final int comp = kmerMapComparator.compare(kmerMaps.get(order[j]),hKm);
                    if (comp < 0) {
                        j++;
                        continue;
                    } else {
                        if (comp == 0) {
                            result[0][next] = j;
                            result[1][next] = i;
                            next++;
                        }
                        break;
                    }
                }
                if (j >= length) {
                    break;
                }
            }
            return result;
        }


        public class ReadMap extends SequenceMap {
            private final GATKSAMRecord read;

            public ReadMap(final GATKSAMRecord r) {
                super(r.getReadBases());
                read = r;
            }


            private ReadHaplotypeMapStats costStats(final Haplotype h) {
                final SequenceMap haplotypeMap = new SequenceMap(h.getBases());
                final int[][] alignment = alignTo(haplotypeMap);
                final ReadHaplotypeMapStats stats = new ReadHaplotypeMapStats();
                stats.length = alignment[0].length;
                if (stats.length > 0) {
                    return stats;
                }
                stats.prefixMissing = Math.max(alignment[0][0],alignment[1][0]);
                stats.suffixMissing = Math.min(length - alignment[0][stats.length - 1] - 1,
                        haplotypeMap.length - alignment[1][stats.length - 1] - 1);
                stats.complexity = stats.length;

                int bubbleArea = 0;
                int rIdx = alignment[0][0];
                int hIdx = alignment[1][0];
                for (int i = 1; i < stats.length; i++) {
                    final int newRIdx = alignment[0][i];
                    final int newHIdx = alignment[1][i];
                    if (newRIdx != rIdx + 1 && newHIdx != hIdx + 1) {
                        bubbleArea += (newRIdx - rIdx) * (newHIdx - hIdx);
                    }
                    rIdx = newRIdx;
                    hIdx = newHIdx;
                }
                stats.complexity += bubbleArea;
                return stats;
            }

        }

        public class ReadHaplotypeMapStats {
            public int length;
            public int prefixMissing;
            public int suffixMissing;
            public int complexity;
        }

    }

    protected SeqVertex mergeLinearChainVertices(final List<SeqVertex> chain) {

        if (chain.size() == 0) {
            throw new IllegalArgumentException("cannot merge an empty chain");
        } else if (chain.size() == 1) {
            return chain.get(0);
        }

        final SeqVertex first = chain.get(0);
        int sequenceLength = first.length();
        for (final SeqVertex v : chain)
            sequenceLength += v.length();

        final byte[] sequence = new byte[sequenceLength];
        int offset = 0;
        for (final SeqVertex sv : chain) {
            int svLength = sv.length();
            System.arraycopy(sv.getSequence(),0,sequence,offset,svLength);
            offset += svLength;
        }
        final SeqVertex result =  new SeqVertex(sequence);
        offset = 0;
        for (final SeqVertex sv : chain) {
            int svLength = sv.length();
            setReplacement(sv,result, offset);
            offset += svLength;
        }
        return result;
    }

    protected void setReplacement(final SeqVertex oldVertex, final SeqVertex newVertex, final int offset) {
        vertexReplacement.put(oldVertex,new Pair<>(newVertex,offset));
    }


}
