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

package org.broadinstitute.sting.utils.pairhmm;

import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: valentin
 * Date: 7/25/13
 * Time: 4:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class MLLog10PairHMM extends Log10PairHMM  implements FlexibleHMM {


    private int readCapacity = 200;
    private int haplotypeCapacity = 400;

    private boolean captureActions = true;
    private Action[][][] actions;
    private double[][][] actionCost;
    private Action bestFinalAction;
    private int bestFinalColumn;
    private boolean usingPairHMMInterface = false;

    enum Action {
        MATCH, INSERTION, DELETION;
    }

    /**
     *
     */
    public MLLog10PairHMM(final byte gcp) {
        super(true);
        constantGCP = gcp;
        initialize(readCapacity,haplotypeCapacity);
        doNotUseTristateCorrection = false;
    }


    @Override
    public double subComputeReadLikelihoodGivenHaplotypeLog10(final byte[] haplotypeBases,
                   final byte[] readBases,
                   final byte[] readQuals,
                   final byte[] insertionGOP,
                   final byte[] deletionGOP,
                   final byte[] overallGCP,
                   final int hapStartIndex,
                   final boolean recacheReadValues, final int nextHapStartIndex ) {
        this.readBases = readBases;
        this.haplotypeBases = haplotypeBases;
        this.usingPairHMMInterface = true;
        final double result = super.subComputeReadLikelihoodGivenHaplotypeLog10(haplotypeBases,readBases,readQuals,insertionGOP,deletionGOP,overallGCP,hapStartIndex,recacheReadValues,nextHapStartIndex);
        this.usingPairHMMInterface = false;
        return result;
    }

    @Override
    public void initialize(final int rc, final int hc) {
        super.initialize(rc,hc);
        if (captureActions && actions == null || actions.length  < rc + 1 || (actions[0].length < hc + 1)) {
            actions = new Action[rc+1][hc+1][3];
            actionCost = new double[rc+1][hc+1][3];
        }
    }


    private boolean graphLikelihoodTest = false;
    public void setupForGraphLikelihoodTest() {
        doNotUseTristateCorrection = false;
        graphLikelihoodTest = true;
    }

    /**
     * Updates a cell in the HMM matrix
     *
     * The read and haplotype indices are offset by one because the state arrays have an extra column to hold the
     * initial conditions

     * @param indI             row index in the matrices to update
     * @param indJ             column index in the matrices to update
     * @param prior            the likelihood editing distance matrix for the read x haplotype
     * @param transition        an array with the six transition relevant to this location
     */
    @Override
    protected void updateCell(final int indI, final int indJ, final double prior, final double[] transition) {
        if (captureActions) {
            updateActionCell(indI,indJ,prior,transition);
        }
        if (!graphLikelihoodTest) {
            matchMatrix[indI][indJ] = prior + Math.max(Math.max(matchMatrix[indI - 1][indJ - 1] + transition[matchToMatch],
                    insertionMatrix[indI - 1][indJ - 1] + transition[indelToMatch]),
                    deletionMatrix[indI - 1][indJ - 1] + transition[indelToMatch]);
            insertionMatrix[indI][indJ] = Math.max(matchMatrix[indI - 1][indJ] + transition[matchToInsertion],
                    insertionMatrix[indI - 1][indJ] + transition[insertionToInsertion]);
            deletionMatrix[indI][indJ] = Math.max(matchMatrix[indI][indJ - 1] + transition[matchToDeletion],
                    deletionMatrix[indI][indJ - 1] + transition[deletionToDeletion]);
            return;
        }
        matchMatrix[indI][indJ] = prior + Math.max(Math.max(matchMatrix[indI - 1][indJ - 1] + (indI == 1 ? 0 : transition[matchToMatch]),
                insertionMatrix[indI - 1][indJ - 1] + transition[indelToMatch]),
                deletionMatrix[indI - 1][indJ - 1] + (indI == 1 ? 0 : transition[indelToMatch]));
        insertionMatrix[indI][indJ] = Math.max(matchMatrix[indI - 1][indJ] + transition[matchToInsertion],
                insertionMatrix[indI - 1][indJ] + transition[insertionToInsertion]);
        deletionMatrix[indI][indJ] = Math.max(matchMatrix[indI][indJ - 1] + transition[matchToDeletion],
                deletionMatrix[indI][indJ - 1] + transition[deletionToDeletion]);
    }

    private void updateActionCell(final int indI, final int indJ, final double prior, final double[] transition) {
        double[] matchCosts = new double[] {
                prior + ((graphLikelihoodTest && indI == 1) ? 0 : transition[matchToMatch]),
                prior + transition[indelToMatch],
                prior + ((graphLikelihoodTest && indI == 1) ? 0 : transition[indelToMatch])
        };

        double[] insertCosts = new double[] {
                transition[matchToInsertion],
                transition[insertionToInsertion],
                Double.NEGATIVE_INFINITY
        };

        double[] deletionCosts = new double[] {
                transition[matchToDeletion],
                Double.NEGATIVE_INFINITY,
                transition[deletionToDeletion]
        };

        double[] matchPredecesors = new double[] {
                matchMatrix[indI - 1][indJ - 1] + matchCosts[0],
                insertionMatrix[indI - 1][indJ - 1] + matchCosts[1],
                deletionMatrix[indI - 1][indJ - 1] + matchCosts[2] };

        double[] insertPredecesors = new double[] {
                matchMatrix[indI - 1][indJ] + insertCosts[0],
                insertionMatrix[indI - 1][indJ] + insertCosts[1],
                deletionMatrix[indI - 1][indJ] + insertCosts[2] };
        double[] deletionPredecesors = new double[] {
                matchMatrix[indI][indJ - 1] + deletionCosts[0],
                insertionMatrix[indI][indJ - 1] + deletionCosts[1],
                deletionMatrix[indI][indJ - 1] + deletionCosts[2]};



        int matchActionOrdinal = maxIndex(matchPredecesors);
        int insertionActionOrdinal = maxIndex(insertPredecesors);
        int deletionActionOrdinal = maxIndex(deletionPredecesors);
        actions[indI][indJ][Action.MATCH.ordinal()] = Action.values()[matchActionOrdinal];
        actions[indI][indJ][Action.INSERTION.ordinal()] = Action.values()[insertionActionOrdinal];
        actions[indI][indJ][Action.DELETION.ordinal()] = Action.values()[deletionActionOrdinal];

        actionCost[indI][indJ][Action.MATCH.ordinal()] = matchCosts[matchActionOrdinal];
        actionCost[indI][indJ][Action.INSERTION.ordinal()] = insertCosts[insertionActionOrdinal];
        actionCost[indI][indJ][Action.DELETION.ordinal()] = deletionCosts[deletionActionOrdinal];

    }

    private int maxIndex(final double[] v) {
        int maxIndex = 0;
        for (int i = 1; i < v.length; i++) {
            if (v[i] > v[maxIndex]) {
                maxIndex = i;
            }
        }
        return maxIndex;
    }

    @Override
    protected double finalLikelihoodCalculation() {
       // if (!graphLikelihoodTest) {
       //     return super.finalLikelihoodCalculation();
       // }
        return finalLikelihoodCalculation(paddedReadLength - 1, 2,
                paddedHaplotypeLength);
    }

    protected double finalLikelihoodCalculation(final int endI,
                                                final int fromCol, final int toCol) {
        double finalLikelihood = matchMatrix[endI][1]; // Math.max(matchMatrix[endI][1],
        bestFinalAction = Action.MATCH;
        bestFinalColumn = 1;
        // insertionMatrix[endI][1]);
        for (int j = fromCol; j < toCol; j++) {
            if (matchMatrix[endI][j] > finalLikelihood) {
                bestFinalAction = Action.MATCH;
                bestFinalColumn = j;
            }
            finalLikelihood = Math.max(finalLikelihood,
                    matchMatrix[endI][j]);// Math.max(matchMatrix[endI][j],
            if (!graphLikelihoodTest) {
                if (insertionMatrix[endI][j] > finalLikelihood) {
                    bestFinalAction = Action.INSERTION;
                    bestFinalColumn = j;
                }
                finalLikelihood = Math.max(finalLikelihood,insertionMatrix[endI][j]);
            }
        }
        if (captureActions && usingPairHMMInterface) {
            calculateAlignment();
        }
        // insertionMatrix[endI][j]));
        return finalLikelihood;
    }

    private String readAlignmentString;
    private String haplotypeAlignmentString;
    private double[] alignmentCosts;

    private void calculateAlignment() {
        final StringBuffer readAlignmentBuffer = new StringBuffer(paddedReadLength + paddedHaplotypeLength);
        final StringBuffer haplotypeAlignmentBuffer = new StringBuffer(paddedReadLength + paddedHaplotypeLength);
        final java.util.List<Double> costs = new ArrayList<>(paddedHaplotypeLength + paddedReadLength);
        for (int i = haplotypeBases.length; i > bestFinalColumn; i--) {
            readAlignmentBuffer.append('-');
            haplotypeAlignmentBuffer.append((char)haplotypeBases[i-1]);
            costs.add(!graphLikelihoodTest ? 0 : transition[readBases.length][deletionToDeletion]);
        }

        int row = readBases.length;
        int col = bestFinalColumn;
        Action lastAction = bestFinalAction;
        while (row > 0 && col > 0) {

            final Action nextAction;
            try {
              nextAction = actions[row][col][lastAction.ordinal()];
            }catch (NullPointerException ex) {
                throw new NullPointerException(" " + actions + " " + row + " " + col + " " + lastAction);
            }
            final double nextCost = actionCost[row][col][lastAction.ordinal()];
            switch (lastAction) {
                case MATCH:
                    readAlignmentBuffer.append((char)readBases[row - 1]);
                    haplotypeAlignmentBuffer.append((char)haplotypeBases[col - 1]);
                    row--; col--;
                    break;
                case INSERTION:
                    readAlignmentBuffer.append((char)readBases[row - 1]);
                    haplotypeAlignmentBuffer.append('-');
                    row--;
                    break;
                case DELETION:
                    readAlignmentBuffer.append('-');
                    haplotypeAlignmentBuffer.append((char)haplotypeBases[col - 1]);
                    col--;
            }
            costs.add(nextCost);
            lastAction = nextAction;
        }
        while (col > 0) {
            readAlignmentBuffer.append('-');
            haplotypeAlignmentBuffer.append((char)haplotypeBases[col - 1]);
            costs.add(0.0);
            col--;
        }
        while (row > 0) {
            readAlignmentBuffer.append((char)readBases[row - 1]);
            haplotypeAlignmentBuffer.append('-');
            costs.add(0.0);
            row--;
        }
        readAlignmentString = readAlignmentBuffer.reverse().toString();
        haplotypeAlignmentString = haplotypeAlignmentBuffer.reverse().toString();
        alignmentCosts = new double[costs.size()];
        for (int i = 0; i < alignmentCosts.length; i++) {
            alignmentCosts[i] = costs.get(costs.size() - i - 1);
        }

    }

    @Override
    protected void initializeMatrixValues(final byte[] haplotypeBases) {
        this.haplotypeBases = haplotypeBases;
        if (!graphLikelihoodTest) {
            super.initializeMatrixValues(haplotypeBases);
            return;
        }
        final double initialValue = 0.0;
        for( int j = 0; j < paddedHaplotypeLength; j++ ) {
            deletionMatrix[0][j] = initialValue;
        }
        matchMatrix[0][0] = 0.0;
        insertionMatrix[0][0] = Double.NEGATIVE_INFINITY;
        for ( int i = 1; i < paddedReadLength; i++ ) {
            insertionMatrix[i][0] = Math.max(insertionMatrix[i-1][0] + transition[i][insertionToInsertion],
                    matchMatrix[i-1][0] + transition[i][matchToInsertion]);
        }
    }

    protected void initializeMatrixValues(final Problem p) {
        final int zeroRow = p.readStart;
        final int zeroCol = p.hapStart;
        final int toRow = p.readEnd;
        final int toCol = p.hapEnd;

        // fill first row with -Inf fot M and I but not for Deletion if leading
        // to allow for free deletions at the beginning.
        if (p.leading) {
            // First row initialization:
            for (int i = zeroCol; i < toCol; i++) {
                matchMatrix[zeroRow][i] = insertionMatrix[zeroRow][i] = Double.NEGATIVE_INFINITY;
                deletionMatrix[zeroRow][i] = 0;
            }
            // First column initialization:
            matchMatrix[zeroRow][zeroCol] = 0;
            for (int i = zeroRow + 1; i < toRow; i++) {
                insertionMatrix[i][zeroCol] = Math.max(
                        transition[i][matchToInsertion]
                                + matchMatrix[i - 1][zeroCol],
                        transition[i][insertionToInsertion]
                                + insertionMatrix[i - 1][zeroCol]);
                matchMatrix[i][zeroCol] = deletionMatrix[i][zeroCol] = Double.NEGATIVE_INFINITY;
            }
        } else { // If not leading set the appropiate matching 1.0 prob and
            // deletion + extension.

            Arrays.fill(matchMatrix[zeroRow], zeroCol + 1, toCol,
                    Double.NEGATIVE_INFINITY);
            Arrays.fill(insertionMatrix[zeroRow], zeroCol, toCol,
                    Double.NEGATIVE_INFINITY);
            matchMatrix[zeroRow][zeroCol] = 0;
            deletionMatrix[zeroRow][zeroCol + 1] = transition[zeroRow][matchToDeletion];
            for (int i = zeroCol + 2; i < toCol; i++) {
                deletionMatrix[zeroRow][i] = deletionMatrix[zeroRow][i - 1]
                        + transition[zeroRow][deletionToDeletion];
            }

            matchMatrix[zeroRow + 1][zeroCol] = deletionMatrix[zeroRow + 1][zeroCol] = Double.NEGATIVE_INFINITY;
            insertionMatrix[zeroRow + 1][zeroCol] = transition[zeroRow + 1][matchToInsertion];

            for (int i = zeroRow + 2; i < toRow; i++) {
                matchMatrix[i][zeroCol] = deletionMatrix[i][zeroCol] = Double.NEGATIVE_INFINITY;
                insertionMatrix[i][zeroCol] = insertionMatrix[i - 1][zeroCol]
                        + transition[i][insertionToInsertion];
            }
        }
    }

    private byte constantGCP;


    private byte[] haplotypeBases;
    private byte[] readBases;
    private byte[] readQuals;
    private Map<Problem, Double> cachedCosts = new HashMap<Problem, Double>();

    @Override
    public void loadRead(final GATKSAMRecord read) {
        loadRead(read.getReadBases(),read.getBaseQualities(),read.getBaseInsertionQualities(),read.getBaseDeletionQualities(),read.getMappingQuality());
    }

    public void loadRead(final byte[] readBases, final byte[] readQuals, final byte[] readInsQuals, final byte[] readDelQuals, int mq) {
        // TODO Copy & paste, kind of akward to share from PairHMM*Engine.

        if (readBases.length > readCapacity) {
            readCapacity = readBases.length;
            initialize(readCapacity,haplotypeCapacity);
        }
        paddedReadLength = readBases.length + 1;
        final byte[] overallGCP = new byte[readBases.length];
        Arrays.fill(overallGCP, constantGCP); // Is there a way to derive
        // empirical estimates for this
        // from the data?
        // NOTE -- must clone anything that gets modified here so we don't screw
        // up future uses of the read

        for (int kkk = 0; kkk < readQuals.length; kkk++) {
            readQuals[kkk] = (byte) Math.min(0xff & readQuals[kkk],
                    mq); // cap base quality by mapping
            // quality, as in UG
            // readQuals[kkk] = ( readQuals[kkk] > readInsQuals[kkk] ?
            // readInsQuals[kkk] : readQuals[kkk] ); // cap base quality by base
            // insertion quality, needs to be evaluated
            // readQuals[kkk] = ( readQuals[kkk] > readDelQuals[kkk] ?
            // readDelQuals[kkk] : readQuals[kkk] ); // cap base quality by base
            // deletion quality, needs to be evaluated
            // TODO -- why is Q18 hard-coded here???
            readQuals[kkk] = (readQuals[kkk] < (byte) 18 ? QualityUtils.MIN_USABLE_Q_SCORE
                    : readQuals[kkk]);
        }
        this.readBases = readBases;
        this.readQuals = readQuals;
        super.initializeProbabilities(readInsQuals, readDelQuals, overallGCP);
        cachedCosts.clear();
    }

    @Override
    public void loadHaplotypeBases(final byte[] hap) {
        if (readBases == null)
            throw new IllegalStateException(
                    "no read was loaded before the haplotype");

        haplotypeBases = hap;
        paddedHaplotypeLength = haplotypeBases.length + 1;
        if (haplotypeCapacity < haplotypeBases.length) {
            haplotypeCapacity = haplotypeBases.length;
            initialize(readCapacity,haplotypeCapacity);
        }
        initializePriors(haplotypeBases, readBases, readQuals, 0);
    }

    public void initializePriors(final byte[] hapBases, final byte[] readBases, final byte[] baseQuals, final int idx) {
        haplotypeBases = hapBases;
        this.readBases = readBases;
        super.initializePriors(hapBases,readBases,baseQuals,idx);
    }

    private long computationSteps = 0;
    private int computationStepsReused;

    public void resetComputationSteps() {
        computationSteps = 0;
    }

    public long getComputationStepsReused() {
        return computationStepsReused;
    }

    public long getComputationSteps() {
        return computationSteps;
    }

    @Override
    public double calculateLocalLikelihood(final int readStart, final int readEnd,
                                           final int hapStart, final int hapEnd, final boolean kmerMatch) {
        if (readBases == null || haplotypeBases == null)
            throw new IllegalStateException("read or haplotype was not loaded");
        final int hapSegmentLength = hapEnd - hapStart;
        final int readSegmentLength = readEnd - readStart;
        double total = 0;
        // trivial case when there is a single base match.
        if (kmerMatch) {
            if (hapSegmentLength == 1) {
                total += prior[readStart + 1][hapStart + 1];
            } else {
                for (int i = 0; i < readSegmentLength; i++) {
                    total += prior[readStart + i + 1][hapStart + i + 1];
                    if (i > 0) {
                        total += transition[readStart + i + 1][matchToMatch];
                    }
                }
            }
            computationSteps += hapSegmentLength;
        } else if (hapSegmentLength == readSegmentLength) {
            if (hapSegmentLength == 0) {
                if (readStart > 0 && readEnd < readBases.length) {
                    total = transition[readStart + 1][matchToMatch];
                }
                computationSteps += 1;
            } else if (hapSegmentLength == 1) {
                total = prior[readStart + 1][hapStart + 1];
                if (readStart > 0) {
                    total += transition[readStart + 1][matchToMatch];
                }
                if (readEnd < readBases.length) {
                    total += transition[readEnd + 1 + 1][matchToMatch];
                }
                computationSteps += 1;
                // another easy.. same length and all are matches:
            } else if (byteEquals(haplotypeBases, hapStart, readBases,
                    readStart, hapSegmentLength)) {

                for (int i = 0; i < readSegmentLength; i++) {
                    total += prior[readStart + i + 1][hapStart + i + 1];
                }
                final int from = Math.max(1, readStart);
                final int to = Math.min(readEnd + 1, readBases.length);
                for (int i = from; i < to; i++) {
                    total += transition[i + 1][matchToMatch];
                }
                computationSteps += hapSegmentLength;
            } else { // general (slower) solution.
                final Problem p = new Problem(readStart, readEnd, hapStart,
                        hapEnd);
                final Double cachedCost = cachedCosts.get(p);
                if (cachedCost != null) {
                    computationStepsReused += hapSegmentLength
                            * readSegmentLength;
                    return cachedCost;
                }
                double cost = doSolve(p);
                cachedCosts.put(p, cost);
                computationSteps += hapSegmentLength * readSegmentLength;
                return cost;
            }
        } else if (hapSegmentLength == 0) { // must be full insertion we
            total += transition[readStart + 1][matchToInsertion];
            for (int i = readStart + 1; i < readEnd; i++) {
                total += transition[i + 1][insertionToInsertion];
            }
            if (readEnd < readBases.length) {
                total += transition[readEnd + 1 + 1][indelToMatch];
            }
            computationSteps += readSegmentLength;
        } else if (readSegmentLength == 0) { // full deletion.
            if (readStart > 0) { // no penalty if at the beginning.
                total += transition[readStart][matchToDeletion];
                total += (hapEnd - hapStart - 1)
                        * transition[readStart][deletionToDeletion];
                total += transition[readStart][indelToMatch];
                computationSteps += hapSegmentLength;
            }
        } else { // general (slower) solution.
            final Problem p = new Problem(readStart, readEnd, hapStart, hapEnd);
            final Double cachedCost = cachedCosts.get(p);
            if (cachedCost != null) {
                computationStepsReused += hapSegmentLength * readSegmentLength;
                return cachedCost;
            }
            double cost = doSolve(p);
            cachedCosts.put(p, cost);
            computationSteps += hapSegmentLength * readSegmentLength;
            return cost;
        }
        return total;
    }

    private boolean byteEquals(final byte[] a, final int aStart,
                               final byte[] b, final int bStart, final int length) {
        if (aStart < 0) {
            throw new IllegalArgumentException("aStart must be greater than 0");
        }
        if (bStart < 0) {
            throw new IllegalArgumentException("bStart must be greater than 0");
        }
        final int aEnd = aStart + length;
        int i, j;
        for (i = aStart, j = bStart; i < aEnd; i++, j++) {
            if (a[i] == b[j]) { // we check this first as it makes it faster.
                continue;
            }
            if (a[i] != 'N' && b[j] != 'N') {
                return false;
            }
        }
        return true;
    }



    private double doSolve(final Problem p) {
        initializeMatrixValues(p);

        updateTableForFullProblem(p.readStart + 1, p.readEnd + 1,
                p.hapStart + 1, p.hapEnd + 1);

        if (p.trailing) {
            return finalLikelihoodCalculation(p.readEnd,2,p.hapEnd + 1);
        } else {
            return Math.max(Math.max(matchMatrix[p.readEnd][p.hapEnd]
                    + transition[p.readEnd][matchToMatch],
                    insertionMatrix[p.readEnd][p.hapEnd]
                            + transition[p.readEnd][indelToMatch]),
                    deletionMatrix[p.readEnd][p.hapEnd]
                            + transition[p.readEnd][indelToMatch]);
        }
    }

    private void updateTableForFullProblem(final int rowFrom, final int rowTo,
                                           final int colFrom, final int colTo) {

        for (int i = rowFrom; i < rowTo; i++) {
            for (int j = colFrom; j < colTo; j++) {
                updateCell(i, j, prior[i][j], transition[i]);
            }
        }

    }

    public class Problem {
        private final byte[] hapSegment;
        private final int readStart;
        private final int readEnd;
        private final int hapStart;
        private final int hapEnd;
        private final int hashCode;
        private final boolean trailing;
        private final boolean leading;

        public Problem(final int start, final int end, final int hapStart,
                       final int hapEnd) {
            if (start < 0 || start > readBases.length) {
                throw new IllegalArgumentException("bad start index " + start);
            }
            if (end < start || end > readBases.length) {
                throw new IllegalArgumentException("bad end index " + end);
            }
            if (hapStart < 0 || hapStart > haplotypeBases.length) {
                throw new IllegalArgumentException("bad hap start index "
                        + hapStart);
            }
            if (hapEnd < hapStart || hapEnd > haplotypeBases.length) {
                throw new IllegalArgumentException("bad hap end index "
                        + hapEnd + " outside [" + hapStart + ","
                        + haplotypeBases);
            }

            hapSegment = Arrays.copyOfRange(haplotypeBases, hapStart, hapEnd);
            readStart = start;
            readEnd = end;
            this.hapStart = hapStart;
            this.hapEnd = hapEnd;
            trailing = readEnd == readBases.length;
            leading = readEnd == 0;

            hashCode = ((start * 31 + end) * 31 + Arrays.hashCode(hapSegment) * 31);
        }

        public int hashCode() {
            return hashCode;
        }

        public boolean equals(Object o) {
            if (o == this) {
                return true;
            } else if (o == null) {
                return false;
            } else if (o.getClass() != this.getClass()) {
                return false;
            } else {
                Problem p = (Problem) o;
                if (p.hashCode != this.hashCode) { // strange to fail and in
                    // that case no big tragedy
                    return false;
                } else {
                    if (p.readStart != this.readStart) {
                        return false;
                    } else if (p.readEnd != this.readEnd) {
                        return false;
                    } else if (Arrays.equals(hapSegment, p.hapSegment)) {
                        return true;
                    } else {
                        return false;
                    }
                }
            }
        }
    }


    public byte[] getReadBases() {
        return readBases;
    }

}
