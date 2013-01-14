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

package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.io.File;
import java.util.ArrayList;

/**
 * BasecallingReadModel represents the statistical models for
 * all bases in all cycles.  It allows for easy training via
 * the addTrainingPoint() method, and for the computation of
 * the 4x4 likelihood matrix or the 1x4 probability vector
 * (with contextual components marginalized out of the
 * likelihood matrix).
 *
 * @author Kiran Garimella
 */
public class BasecallingReadModel {
    private BasecallingBaseModel[] basemodels = null;
    private boolean correctForContext = true;

    /**
     * Constructs a BasecallingReadModel with space for a given read length.
     * 
     * @param readLength  the length of the reads to which this model will apply.
     */
    public BasecallingReadModel(int readLength) {
        initialize(readLength);
    }

    /**
     * Constructs a BasecallingReadModel and trains it using the specified training data.
     *
     * @param trainingData  a set of RawReads from which the model will be trained.
     */
    public BasecallingReadModel(ArrayList<RawRead> trainingData) {
        initialize(trainingData.get(0).getReadLength());

        train(trainingData);
    }

    /**
     * Initialize the model and set default parameters for each cycle appropriately.
     *
     * @param readLength  the length of the reads to which this model will apply.
     */
    public void initialize(int readLength) {
        basemodels = new BasecallingBaseModel[readLength];

        for (int cycle = 0; cycle < readLength; cycle++) {
            basemodels[cycle] = new BasecallingBaseModel(cycle != 0 && correctForContext);
        }
    }

    /**
     * Train the model using the specified training data.
     *
     * @param trainingData  a set of RawReads from which the model will be trained.
     */
    public void train(ArrayList<RawRead> trainingData) {
        for ( RawRead read : trainingData ) {
            addMeanPoints(read);
        }

        for ( RawRead read : trainingData ) {
            addCovariancePoints(read);
        }
    }

    /**
     * Add a training point for the mean intensity values per base and per cycle.
     *
     * @param cycle          the cycle number (0-based)
     * @param probMatrix     the probability matrix for the base
     * @param fourintensity  the four raw intensities for the base
     */
    public void addMeanPoint(int cycle, double[][] probMatrix, double[] fourintensity) {
        basemodels[cycle].addMeanPoint(probMatrix, fourintensity);
    }

    /**
     * Add a training point for the mean intensity values per base in all cycles.
     *
     * @param read  the raw read
     */
    public void addMeanPoints(RawRead read) {
        byte[] seqs = read.getSequence();
        byte[] quals = read.getQuals();
        short[][] ints = read.getIntensities();

        for (int cycle = 0; cycle < seqs.length; cycle++) {
            char basePrev = (char) ((cycle == 0) ? '.' : seqs[cycle - 1]);
            char baseCur = (char) seqs[cycle];
            double probCur = QualityUtils.qualToProb(quals[cycle]);

            double[][] probMatrix = getBaseProbabilityMatrix(cycle, basePrev, baseCur, probCur);

            double[] fourIntensity = new double[4];
            for (int channel = 0; channel < 4; channel++) {
                //fourIntensity[channel] = (double) ints[cycle][channel];
                fourIntensity[channel] = (double) ints[channel][cycle];
            }

            basemodels[cycle].addMeanPoint(probMatrix, fourIntensity);
        }
    }

    /**
     * Add a training point for the intensity covariance matrix per base and per cycle.
     *
     * @param cycle          the cycle number (0-based)
     * @param probMatrix     the probability matrix for the base
     * @param fourintensity  the four raw intensities for the base
     */
    public void addCovariancePoint(int cycle, double[][] probMatrix, double[] fourintensity) {
        basemodels[cycle].addCovariancePoint(probMatrix, fourintensity);
    }

    /**
     * Add a training point for the intensity covariance matrix per base in all cycles.
     *
     * @param read  the raw read
     */
    public void addCovariancePoints(RawRead read) {
        byte[] seqs = read.getSequence();
        byte[] quals = read.getQuals();
        short[][] ints = read.getIntensities();

        for (int cycle = 0; cycle < seqs.length; cycle++) {
            char basePrev = (char) ((cycle == 0) ? '.' : seqs[cycle - 1]);
            char baseCur = (char) seqs[cycle];
            double probCur = QualityUtils.qualToProb(quals[cycle]);

            double[][] probMatrix = getBaseProbabilityMatrix(cycle, basePrev, baseCur, probCur);

            double[] fourIntensity = new double[4];
            for (int channel = 0; channel < 4; channel++) {
                //fourIntensity[channel] = (double) ints[cycle][channel];
                fourIntensity[channel] = (double) ints[channel][cycle];
            }

            basemodels[cycle].addCovariancePoint(probMatrix, fourIntensity);
        }
    }

    /**
     * Compute the likelihoods that a given set of intensities yields each possible base.
     *
     * @param cycle          the cycle number (0-based)
     * @param fourintensity  the four raw intensities for the base
     * @return               the matrix of likelihoods
     */
    public double[][] computeLikelihoods(int cycle, double[] fourintensity) {
        return basemodels[cycle].computeLikelihoods(cycle, fourintensity);
    }

    /**
     * Compute the probabilities that a given set of intensities yields each possible base.
     *
     * @param cycle          the cycle number (0-based)
     * @param basePrev       the previous base
     * @param qualPrev       the previous base's quality score
     * @param fourintensity  the four raw intensities for the base
     * @return the probability distribution over the four base possibilities
     */
    public FourProb computeProbabilities(int cycle, char basePrev, byte qualPrev, double[] fourintensity) {
        double[][] likes = computeLikelihoods(cycle, fourintensity);

        double total = 0;

        for (int basePrevIndex = 0; basePrevIndex < likes.length; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                double prior = 1.0;
                if (correctForContext) {
                    double prob = QualityUtils.qualToProb(qualPrev);
                    if (basePrevIndex == BaseUtils.simpleBaseToBaseIndex(basePrev)) {
                        prior = prob;
                    } else {
                        prior = (1.0 - prob)/((double) (4*likes.length - 1));
                    }
                }
                likes[basePrevIndex][baseCurIndex] = prior*likes[basePrevIndex][baseCurIndex];
                total += likes[basePrevIndex][baseCurIndex];
            }
        }

        for (int basePrevIndex = 0; basePrevIndex < likes.length; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                likes[basePrevIndex][baseCurIndex] /= total;
            }
        }

        return new FourProb(likes);
    }

    /**
     * Call the bases in the given RawRead.
     *
     * @param read  the RawRead
     * @return  the basecalled read
     */
    public FourProbRead call(RawRead read) {
        FourProbRead fpr = new FourProbRead(read.getReadLength());
        
        for (int cycle = 0; cycle < read.getReadLength(); cycle++) {
            char basePrev = (char) ((cycle == 0) ? '.' : read.getSequence()[cycle - 1]);
            byte qualPrev = ((cycle == 0) ? 0 : read.getQuals()[cycle - 1]);

            double[] fourIntensity = new double[4];
            for (int channel = 0; channel < 4; channel++) {
                //fourIntensity[channel] = (double) read.getIntensities()[cycle][channel];
                fourIntensity[channel] = (double) read.getIntensities()[channel][cycle];
            }

            fpr.add(cycle, computeProbabilities(cycle, basePrev, qualPrev, fourIntensity));

        }

        return fpr;
    }

    /**
     * Return the probability matrix given the previous cycle's base, the current cycle's base, and the current base's probability.
     *
     * @param cycle     the cycle number (0-based)
     * @param basePrev  the previous base
     * @param baseCur   the current base
     * @param probCur   the probability of the current base
     * @return the probability matrix of the base
     */
    public double[][] getBaseProbabilityMatrix(int cycle, char basePrev, char baseCur, double probCur) {
        double[][] dist = new double[(correctForContext && cycle > 0) ? 4 : 1][4];

        int actualBasePrevIndex = (correctForContext && cycle > 0) ? BaseUtils.simpleBaseToBaseIndex(basePrev) : 0;
        int actualBaseCurIndex = BaseUtils.simpleBaseToBaseIndex(baseCur);

        if (actualBasePrevIndex == -1) { actualBasePrevIndex = BaseUtils.getRandomBaseIndex(); }
        if (actualBaseCurIndex == -1) { actualBaseCurIndex = BaseUtils.getRandomBaseIndex(); }

        double residualTheories = (double) (dist.length*dist[0].length - 1);

        for (int basePrevIndex = 0; basePrevIndex < dist.length; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < dist[basePrevIndex].length; baseCurIndex++) {
                dist[basePrevIndex][baseCurIndex] = (basePrevIndex == actualBasePrevIndex && baseCurIndex == actualBaseCurIndex) ? probCur : ((1.0 - probCur)/residualTheories);
            }
        }

        return dist;
    }

    /**
     * Write model parameters to disk.
     *
     * @param dir  the directory in which model parameters should be stored.
     */
    public void write(File dir) {
        for (int cycle = 0; cycle < basemodels.length; cycle++) {
            File outparam = new File(dir.getPath() + "/param." + cycle + ".r");
            basemodels[cycle].write(outparam);
        }
    }
}
