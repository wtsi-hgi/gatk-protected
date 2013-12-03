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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.pairhmm.FlexibleHMM;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * Maximum likelihood calculation engine, a single 1-D array PairHMM implementation.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.com&gt;
 *
 */
public class MLECalculationEngine implements FlexibleHMM {


    private long computationSteps = 0;
    private int computationStepsReused;

    /**
     * Resets the computation step counts.
     */
    public void resetComputationSteps() {
        computationSteps = 0;
    }

    /**
     * Returns the number of computation step reused since the last reset.
     *
     * <p>
     * A computation step is a reused one when its CPU cost was avoided by using cached values.
     * </p>
     *
     * @return 0 or greater.
     */
    public long getComputationStepsReused(){
        return computationStepsReused;
    }

    /**
     * Total number of computation steps either those consumed CPU time and those that were spared by using a cached value.
     *
     * @return 0 or greater.
     */
    public long getComputationSteps() {
        return computationSteps;
    }

    /**
     * Returns the bases of the current read.
     *
     * @return null if no read has yet been loaded.
     */
    public byte[] getReadBases() {
        return base;
    }


    /**
     * Perform the partial PairHMM maximum-likelihood calculation of a section of a read and haplotype.
     *
     * <p>It makes that assumption that the section of the table is preceded and followed by a match
     * between the read and the haplotype unless the segment correspond to the beginning or end of the
     * either haplotype or read.</p>
     *
     * @param readStart 0-based inclusive start on the read.
     * @param readEnd 0-based exclusive end on the read.
     * @param hapStart 0-based inclusive start on the haplotype.
     * @param hapEnd 0-based exclusive end on the haplotype.
     * @return log-scaled MLE for that pair-hmm section. Thus the value should be 0 or less.
     */
    public double fullPairHMM(final int readStart, final int readEnd, final int hapStart, final int hapEnd) {
        return calculateLocalLikelihood(readStart, readEnd, hapStart, hapEnd, false);
    }

    /**
     * Approximates the partial PairHMM maximum-likelihood calculation of a section of a read and haplotype.
     *
     * <p>It uses matches between unique kmer shared by the read and haplotype to constraint the maximum-likelihood
     * path to pass through </p>
     *
     * @param readStart 0-based inclusive start on the read.
     * @param readEnd 0-based exclusive end on the read.
     * @param hapStart 0-based inclusive start on the haplotype.
     * @param hapEnd 0-based exclusive end on the haplotype.
     * @param alignment match between read kmers and haplotype kmers. This array should have length at least
     *                  {@code readEnd - readStart} and its {@code i}th position indicates what kmer (0-based) in the
     *                  haplotype correspond to the {@code (i + readStart)}th kmer on the read.
     *
     * @return log-scaled MLE for that pair-hmm section. Thus the value should be 0 or less.
     */
    public double fastPairHMM(final int readStart, final int readEnd, final int hapStart,
                              final int hapEnd, final int kmerSize, final int[] alignment) {
        int hapSegStart = hapStart;
        int readSegStart = readStart;
        double cost = 0;
        while (readSegStart < readEnd) {
            if (alignment[readSegStart - readStart] == -1) {
                int readSegEnd = readSegStart;
                while (readSegEnd < readEnd && alignment[readSegEnd - readStart] == -1) {
                    readSegEnd++;
                }
                final int hapSegEnd = readSegEnd < readEnd ? alignment[readSegEnd - readStart]  : hapEnd;
                cost += calculateLocalLikelihood(readSegStart, readSegEnd, hapSegStart, hapSegEnd, false);
                hapSegStart = hapSegEnd;
                readSegStart = readSegEnd;
            } else {
                int readSeqEnd = readSegStart;
                hapSegStart = alignment[readSegStart - readStart];
                while (readSeqEnd < readEnd - 1 && alignment[readSeqEnd + 1 - readStart]  == alignment[readSeqEnd - readStart] + 1) {
                    readSeqEnd++;
                }
                int hapSegEnd = alignment[readSeqEnd - readStart] + 1;
                readSeqEnd++;
                hapSegEnd += kmerSize - 1;   // be greedy in the case of a kmer match stretch.
                readSeqEnd += kmerSize - 1;   // be greedy in the case of a kmer match stretch.
                cost += calculateLocalLikelihood(readSegStart, readSeqEnd, hapSegStart, hapSegEnd, true);
                hapSegStart = hapSegEnd;
                readSegStart = readSeqEnd;
                if (readSegStart < readEnd && alignment[readSegStart - readStart] != -1) {
                    hapSegEnd = Math.max(alignment[readSegStart - readStart], hapSegStart);
                    cost += calculateLocalLikelihood(readSegStart, readSegStart, hapSegStart, hapSegEnd, false);
                    hapSegStart = hapSegEnd;
                }
            }
        }
        return cost;
    }



    public MLECalculationEngine(final byte[] bases, final byte[] bq, final byte[] iq, final byte[] dq, int mq) {
        initCostCalculationStructures(bases, bq, iq, dq, mq);
    }

    public MLECalculationEngine(final GATKSAMRecord read) {
        initCostCalculationStructures(read.getReadBases(), read.getBaseQualities(), read.getBaseInsertionQualities(), read.getBaseDeletionQualities(), read.getMappingQuality());
    }

    @Override
    public void loadRead(final GATKSAMRecord read) {
        initCostCalculationStructures(read.getReadBases(), read.getBaseQualities(), read.getBaseInsertionQualities(), read.getBaseDeletionQualities(), read.getMappingQuality());
    }


    @Override
    public void loadRead(final byte[] bases, final byte[] bq, final byte[] iq, final byte[] dq, final int mq) {
        initCostCalculationStructures(bases, bq, iq, dq, mq);
    }



    protected Map<Problem, Double> cachedLikelihoods;
    private static final int MtoM = 0;
    private static final int ItoM = 1;
    private static final int DtoM = 1;
    private static final int DtoD = 2;
    private static final int MtoD = 3;
    private static final int ItoI = 4;
    private static final int MtoI = 5;
    private static final int TRANSITION_CELL_SIZE = 6;


    private static final int MATRIX_CELL_SIZE = 3;
    private static final int M = 0;
    private static final int I = 1;
    private static final int D = 2;

    private byte[] haplotype;
    private byte[] base;
    private byte[] bq;
    private double[] bMatchProb;
    private double[] bMissProb;
    private byte[] iq;
    private byte[] dq;
    private double[] trans;
    private int length;
    private byte gapExtensionPenalty = 10;
    private double minGapOpenCost;

    private void initCostCalculationStructures(final byte[] bases, final byte[] bq, final byte[] iq, final byte[] dq, final int mq) {
        if (base != null) {
            return;
        }
        minGapOpenCost = 0;
        base = bases;
        this.bq = bq.clone();
        for (int i = 0; i < this.bq.length; i++) {
            this.bq[i] = (byte) Math.min(0xff & bq[i], mq);
            this.bq[i] = ( this.bq[i] < (byte) 18 ? QualityUtils.MIN_USABLE_Q_SCORE : this.bq[i]);
        }

        final double log10_3 = StrictMath.log10(3.0);
        bMatchProb = new double[bq.length];
        bMissProb = new double[bq.length];
        this.iq = iq;
        this.dq = dq;
        length = bases.length;
        trans = new double[length * TRANSITION_CELL_SIZE];
        for (int i = 0, offset = 0; i < length; i++, offset += TRANSITION_CELL_SIZE) {
            final int qualIndexGOP = Math.min(iq[i] + dq[i], Byte.MAX_VALUE);
            if (qualIndexGOP < minGapOpenCost) {
                minGapOpenCost = qualIndexGOP;
            }
            bMatchProb[i] = QualityUtils.qualToProbLog10(this.bq[i]);
            bMissProb[i] = QualityUtils.qualToErrorProbLog10(this.bq[i]) - log10_3;
            trans[offset + MtoM] = QualityUtils.qualToProbLog10((byte) qualIndexGOP);
            trans[offset + ItoM] = QualityUtils.qualToProbLog10(gapExtensionPenalty);
            trans[offset + MtoI] = QualityUtils.qualToErrorProbLog10(iq[i]);
            trans[offset + ItoI] = QualityUtils.qualToErrorProbLog10(gapExtensionPenalty);
            trans[offset + MtoD] = QualityUtils.qualToErrorProbLog10(dq[i]);
            trans[offset + DtoD] = QualityUtils.qualToErrorProbLog10(gapExtensionPenalty);
        }

        cachedLikelihoods = new HashMap<>();
    }

    @Override
    public void loadHaplotypeBases(final byte[] haplotypeSequence) {
        haplotype = haplotypeSequence;
    }


    /**
     * Calculates the likelihood of a local alignment between a segment of a read vs a segment of the haplotype.
     *
     * <p>
     *     The read and haplotypes should have been loaded previously using {@link #loadRead}
     *     and {@link #loadHaplotypeBases} respectively.
     * </p>
     *
     * <p>
     *     When {@code assumeIsExactMatch} is true two things happen: firstly, ats the parameter name imply the
     *     algorithm is free to assume that indeed both segments have the same nucleotides and secondly there is no
     *     further penalties added before and after (e.g. match-to-match, match-to-indel etc).
     *     In contrast when {@code assumeIsExactMatch} is false the algorithm won't make such assumption thus probably
     *     resulting a regular PairHMM calculation (but not necessarily so) and if the segments are not located to
     *     the beginning or end of the read
     *     transition penalties will be added assuming that what came before or after is a match.
     * </p>
     *
     *
     * @param readStart 0-based inclusive start position on the read.
     * @param readEnd 0-based exclusive end position on the read.
     * @param hapStart 0-based inclusive start position on the haplotype.
     * @param hapEnd 0-based exclusive end position on the haplotype.
     * @param assumeIsExactMatch whether we can assume that there is a perfect match between both read and haplotype
     *                           segments.
     * @return 1 or less.
     */
    public double calculateLocalLikelihood(final int readStart, final int readEnd, final int hapStart, final int hapEnd,
                                           final boolean assumeIsExactMatch) {
        final int hapSegmentLength = hapEnd - hapStart;
        final int readSegmentLength = readEnd - readStart;

        if (assumeIsExactMatch) {
            return calculateLocalLikelihoodWithExactMatch(readStart, readEnd);
        } else if (hapSegmentLength == readSegmentLength) {
            if (hapSegmentLength == 0)
                return calculateLocalLikelihoodEmptySquare(readStart, readEnd);
            else if (hapSegmentLength == 1)
                return calculateLocalLikelihoodSingleBase(readStart, hapStart);
            else if (BaseUtils.basesAreEqualIgnoreAmbiguous(haplotype, hapStart, base, readStart, hapSegmentLength))
                return calculateLocalLikelihoodSquareEqualBases(readStart, readEnd);
            else
                return calculateLocalLikelihoodsGeneralized(readStart, readEnd, hapStart, hapEnd);
        } else if (hapSegmentLength == 0)
            return calculateLocalLikelihoodsInsertion(readStart, readEnd);
        else if (readSegmentLength == 0)
            return readStart <= 0 ? 0 :calculateLocalLikelihoodsDeletion(readStart, readEnd, hapStart, hapEnd);
        else
            return calculateLocalLikelihoodsGeneralized(readStart, readEnd, hapStart, hapEnd);
    }

    private double calculateLocalLikelihoodsDeletion(final int readStart, final int readEnd,
                                                     final int hapStart, final int hapEnd) {
        double total = 0;
        int transOffset = (readStart - 1) * TRANSITION_CELL_SIZE;
        total += trans[transOffset + MtoD];
        for (int i = hapStart + 1; i < hapEnd; i++) {
            total += trans[transOffset + DtoD];
        }
        if (readEnd < base.length) {
            total += trans[transOffset + DtoM];
        }
        computationSteps += hapStart - hapEnd;
        return total;
    }

    private double calculateLocalLikelihoodsInsertion(final int readStart, final int readEnd) {
        double total = 0;
        int transOffset = readStart * TRANSITION_CELL_SIZE;
        total += trans[transOffset + MtoI];
        transOffset += TRANSITION_CELL_SIZE;
        for (int i = readStart + 1; i < readEnd; i++, transOffset += TRANSITION_CELL_SIZE) {
            total += trans[transOffset + ItoI];
        }
        if (readEnd < base.length) {
            total += trans[transOffset + ItoM];
        }
        computationSteps += readEnd - readStart;
        return total;
    }

    /**
     * Resolve a general problem using a PairHMM.
     */
    private double calculateLocalLikelihoodsGeneralized(final int readStart, final int readEnd, final int hapStart,
                                                        final int hapEnd) {

        final int hapSegmentLength = hapEnd - hapStart;
        final int readSegmentLength = readEnd - readStart;
        computationSteps += hapSegmentLength * readSegmentLength;
        final Problem p = new Problem(readStart, readEnd, hapStart, hapEnd);
        final Double cachedCost = cachedLikelihoods.get(p);
        if (cachedCost != null) {
            computationStepsReused += hapSegmentLength * readSegmentLength;
            return cachedCost;
        }
        final double cost = solve(p);
        cachedLikelihoods.put(p, cost);
        return cost;
    }

    /**
     * Calculates the likelihood of a local squared alignment (same length for read and hap segements) assuming
     * is surrounded by base matches (unless at the beginning or end of the read).
     */
    private double calculateLocalLikelihoodSquareEqualBases(final int readStart, final int readEnd) {
        double total = 0;
        for (int i = readStart; i < readEnd; i++)
            total += bMatchProb[readStart];
        final int from = Math.max(1, readStart);
        final int to = Math.min(readEnd + 1, base.length);
        int transOffsetMtoM = from * TRANSITION_CELL_SIZE + MtoM;
        for (int i = from; i < to; i++, transOffsetMtoM += TRANSITION_CELL_SIZE) {
            total += trans[transOffsetMtoM];
        }
        computationSteps += readEnd - readStart;
        return total;
    }

    /**
     * Calculates the likelihood of a single base pairwise alignment assuming it is surrounded by base matches
     * (unless at the end or beginning of the read).
     */
    private double calculateLocalLikelihoodSingleBase(final int readIndex, final int hapIndex) {
        double total;
        final byte readBase = base[readIndex];
        final byte hapBase = haplotype[hapIndex];
        final boolean match = (readBase == hapBase || readBase == 'N' || hapBase == 'N');
        total = match ? bMatchProb[readIndex] : bMissProb[readIndex];
        if (readIndex > 0) {
            total += trans[readIndex * TRANSITION_CELL_SIZE + MtoM];
        }
        if (readIndex + 1 < base.length) {
            total += trans[(readIndex + 1) * TRANSITION_CELL_SIZE + MtoM];
        }
        computationSteps += 1;
        return total;
    }

    /**
     * Local alignment likelihood of an empty square assuming is surrounded by base matches.
     *
     *
     * @param readStart
     * @param readEnd
     * @return 0 or less.
     */
    private double calculateLocalLikelihoodEmptySquare(final int readStart, final int readEnd) {
        computationSteps += 1;
        return readStart > 0 && readEnd < base.length ? trans[readStart * TRANSITION_CELL_SIZE + MtoM] : 0;
    }

    /**
     * Local alignment likelihood of an exact match section between read and haplotype.
     *
     * @param readStart
     * @param readEnd
     * @return
     */
    private double calculateLocalLikelihoodWithExactMatch(final int readStart, final int readEnd) {
        final int length = readEnd - readStart;
        computationSteps += length;
        if (length == 1)
            return bMatchProb[readStart];
        else {
            int transOffset = readStart * TRANSITION_CELL_SIZE + MtoM;
            int total = 0;
            for (int i = readStart; i < readEnd; i++, transOffset += TRANSITION_CELL_SIZE) {
                total += bMatchProb[i];
                if (i > readStart) {
                    total += trans[transOffset];
                }
            }
            return total;
        }
    }

    /**
     * Calculate the likelihood of a local alignment.
     *
     * <p>
     *     If the same problem has been solved before then it will return the catched result value.
     * </p>
     *
     * @param p
     * @return 0 or less.
     */
    private double solve(final Problem p) {

        p.allocateMatrix();

        int offset = p.width; // skip first column.
        int baseOffset = p.seqStart;
        int transOffset = p.startTransOffset;
        for (int i = 1; i < p.rowCount; i++, baseOffset++, transOffset += TRANSITION_CELL_SIZE) {
            offset += MATRIX_CELL_SIZE; // skip first column
            for (int j = 1; j < p.columnCount; j++, offset += MATRIX_CELL_SIZE) {
                solve(p, offset, i, j, base[baseOffset], p.path[j - 1], bMatchProb[baseOffset], bMissProb[baseOffset], transOffset);
            }
        }
        offset -= MATRIX_CELL_SIZE; // offset: from just beyond the matrix last position, go back one cell to top right cell offset.
        // transOffset : points to the next base transition probs that must always be a match by construction.


        // if trailing we only allow for the last state to be a match or deletion costless.
        // if non-trailing it can be anything but need to penalize as if it is followed my a match.

        double cost = p.trailing ? Math.max(Math.max(p.matrix[offset + M], p.matrix[offset + D]), p.matrix[offset + I])
                : Math.max(Math.max(p.matrix[offset + M] + trans[transOffset + MtoM],
                p.matrix[offset + I] + trans[transOffset + ItoM]),
                p.matrix[offset + D] + trans[transOffset + DtoM]);
        p.freeMatrix();
        return cost;
    }

    private void solve(final Problem p, final int offset, final int r, final int c,
                       final byte seqBase, final byte haplotypeBase, final double matchProb,
                       final double missProb, final int transOffset) {
        final int offsetM = offset + M;
        final int offsetI = offset + I;
        final int offsetD = offset + D;

        // in the last row, deletions are free for trailing sequence problems:
        if (p.trailing && r == p.rowCount - 1 && c >= 2) {
            p.matrix[offsetD] = Math.max(p.matrix[offsetM - MATRIX_CELL_SIZE], p.matrix[offsetD - MATRIX_CELL_SIZE]);
        } else {
            p.matrix[offsetD] = Math.max(p.matrix[offsetM - MATRIX_CELL_SIZE] + trans[transOffset + MtoD],
                    p.matrix[offsetD - MATRIX_CELL_SIZE] + trans[transOffset + DtoD]);
        }

        final double prior = seqBase == haplotypeBase || seqBase == 'N' || haplotypeBase == 'N' ? matchProb : missProb;

        final int matchOffsetDiff = p.width + MATRIX_CELL_SIZE;

        p.matrix[offsetM] = prior +
                Math.max(Math.max(p.matrix[offsetM - matchOffsetDiff] + (p.leading && r == 1 ? 0 : trans[transOffset + MtoM]),
                        p.matrix[offsetD - matchOffsetDiff] + (p.leading && r == 1 ? 0 : trans[transOffset + DtoM])),
                        p.matrix[offsetI - matchOffsetDiff] + trans[transOffset + ItoM]);

        final int insertOffsetDiff = p.width;

        p.matrix[offsetI] = Math.max(p.matrix[offsetM - insertOffsetDiff] + trans[transOffset + MtoI],
                p.matrix[offsetI - insertOffsetDiff] + trans[transOffset + ItoI]);
    }


    /**
     * Represents the information that characterize a computational problem to be solved such as the haplotype and read
     * segment to pair-hmm.
     */
    private class Problem {
        private final byte[] path;
        private final int seqStart;
        private final int seqEnd;
        private transient final int width;
        private transient final int columnCount;
        private transient final int rowCount;
        private transient final int height;
        private transient final int hashCode;
        private transient final boolean trailing;
        private transient final boolean leading;
        private transient final int startTransOffset;

        /**
         * the matrix is organized rows in the <MID> order so
         * <p/>
         * I state max at row r and column c is located at
         * r * (3* columnCount) + 3 * c + 1;
         * <p/>
         * heigth = r * 3;
         * height = rowCount * 3;
         * width = columnCount * 3;
         * MATRIX_CELL_SIZE = 3
         * <p/>
         * r * width + MATRIX_CELL_SIZE * c + 1;
         */
        private transient double[] matrix;

        /**
         * Construct a new problem given the read and haplotype segment ranges.
         */
        private Problem(final int readStart, final int readEnd, final int hapStart, final int hapEnd) {
            if (readStart < 0 || readStart > base.length) {
                throw new IllegalArgumentException("bad readStart index " + readStart);
            }
            if (readEnd < readStart || readEnd > base.length) {
                throw new IllegalArgumentException("bad readEnd index " + readEnd);
            }
            if (hapStart < 0 || hapStart > haplotype.length) {
                throw new IllegalArgumentException("bad hapStart index " + hapStart);
            }
            if (hapEnd < hapStart || hapEnd > haplotype.length) {
                throw new IllegalArgumentException("bad hapEnd index " + hapEnd + " outside [" + hapStart + "," + haplotype.length + "]");
            }

            path = Arrays.copyOfRange(haplotype, hapStart, hapEnd);
            seqStart = readStart;
            seqEnd = readEnd;
            columnCount = path.length + 1;
            rowCount = readEnd - readStart + 1;
            width = columnCount * 3;
            height = rowCount * 1;

            hashCode = (readStart + 1) + (readEnd + 1) + Arrays.hashCode(path);
            leading = seqStart == 0;
            trailing = seqEnd == base.length;
            startTransOffset = seqStart * TRANSITION_CELL_SIZE;
        }

        /**
         * Initializes the matrix's first row and column of 3-states.
         */
        protected void allocateMatrix() {
            // create the calculation matrix.
            final int area = width * height;
            matrix = new double[area];


            // fill first row with -Inf fot M and I but not for Deletion if leading
            // to allow for free deletions at the beginning.
            if (leading) {
                // First row initialization:
                for (int offset = 0; offset < width; offset += MATRIX_CELL_SIZE) {
                    matrix[offset + M] = matrix[offset + I] = Double.NEGATIVE_INFINITY;
                    matrix[offset + D] = 0;
                }
                // First column initialization:
                matrix[M] = 0;
                int transOffset = 0;
                for (int offset = width; offset < area; offset += width, transOffset += TRANSITION_CELL_SIZE) {
                    matrix[offset + I] = Math.max(trans[transOffset + MtoI] + matrix[offset + M - width], trans[transOffset + ItoI] + matrix[offset + I - width]);
                    matrix[offset + M] = matrix[offset + D] = Double.NEGATIVE_INFINITY;
                }


            } else { // If not leading set the appropiate matching 1.0 prob and deletion + extension.

                // trans probabilities for the preceeding character.
                // only make sense (and is ussed if leading if leading
                final int prevMBaseOffset = (seqStart - 1) * TRANSITION_CELL_SIZE;

                // Row initialization:
                Arrays.fill(matrix, 0, width, Double.NEGATIVE_INFINITY);
                matrix[M] = 0;
                matrix[MATRIX_CELL_SIZE + D] = trans[prevMBaseOffset + MtoD];
                for (int offset = (MATRIX_CELL_SIZE << 1) + D; offset < width; offset += MATRIX_CELL_SIZE) {
                    matrix[offset] = matrix[offset - MATRIX_CELL_SIZE] + trans[prevMBaseOffset + DtoD];
                }
                //Column initialization.
                //matrix[M] = 0; done before.
                int offset = width;
                matrix[offset + D] = matrix[offset + M] = Double.NEGATIVE_INFINITY;
                matrix[offset + I] = /* matrix[M] + // is zero so... */ trans[prevMBaseOffset + MtoI];
                for (offset += width; offset < area; offset += width) {
                    matrix[offset + M] = matrix[offset + D] = Double.NEGATIVE_INFINITY;
                    matrix[offset + I] = matrix[offset - width + I] + trans[prevMBaseOffset + ItoI];
                }

            }

        }

        // Free the space require to compute the result to this problem.
        private void freeMatrix() {
            matrix = null;
        }

        @Override
        public int hashCode() {
            return hashCode;
        }

        @Override
        public boolean equals(Object o) {
            if (o == null) {
                return false;
            } else if (o == this) {
                return true;
            } else if (o.hashCode() != this.hashCode) {
                return false;
            } else if (o.getClass() != this.getClass()) {
                return false;
            } else {
                final Problem other = (Problem) o;
                if (this.seqStart != other.seqStart) {
                    return false;
                } else if (this.seqEnd != other.seqEnd) {
                    return false;
                } else if (Arrays.equals(path, other.path)) {
                    return true;
                } else {
                    return false;
                }
            }
        }
    }

}
