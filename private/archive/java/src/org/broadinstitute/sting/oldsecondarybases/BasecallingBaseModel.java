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

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * BasecallingBaseModel is a class that represents the statistical
 * model for all bases at a given cycle.  It allows for easy, one
 * pass training via the addTrainingPoint() method.  Once the model
 * is trained, computeLikelihoods will return the probability matrix
 * over previous cycle's base hypotheses and current cycle base
 * hypotheses (contextual prior is included in these likelihoods).
 *
 * @author Kiran Garimella
 */
public class BasecallingBaseModel {
    private double[][] counts;
    private DoubleMatrix1D[][] sums;
    private DoubleMatrix2D[][] unscaledCovarianceSums;

    private DoubleMatrix1D[][] means;
    private DoubleMatrix2D[][] inverseCovariances;
    private double[][] norms;

    private cern.jet.math.Functions F = cern.jet.math.Functions.functions;
    private Algebra alg;

    private boolean correctForContext = false;
    private int numTheories = 1;

    private boolean readyToCall = false;
    private boolean bustedCycle = false;

    /**
     * Constructor for BasecallingBaseModel.
     *
     * @param correctForContext  should we attempt to correct for contextual sequence effects?
     */
    public BasecallingBaseModel(boolean correctForContext) {
        this.correctForContext = correctForContext;
        this.numTheories = (correctForContext) ? 4 : 1;
        
        counts = new double[this.numTheories][4];

        sums = new DoubleMatrix1D[this.numTheories][4];
        unscaledCovarianceSums = new DoubleMatrix2D[this.numTheories][4];

        means = new DoubleMatrix1D[this.numTheories][4];
        inverseCovariances = new DoubleMatrix2D[this.numTheories][4];
        norms = new double[this.numTheories][4];

        for (int basePrevIndex = 0; basePrevIndex < this.numTheories; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                sums[basePrevIndex][baseCurIndex] = (DoubleFactory1D.dense).make(4);
                unscaledCovarianceSums[basePrevIndex][baseCurIndex] = (DoubleFactory2D.dense).make(4, 4);

                means[basePrevIndex][baseCurIndex] = (DoubleFactory1D.dense).make(4);
                inverseCovariances[basePrevIndex][baseCurIndex] = (DoubleFactory2D.dense).make(4, 4);
            }
        }

        alg = new Algebra();
    }


    /**
     * Add a single training point to the model to estimate the means.
     *
     * @param probMatrix     the matrix of probabilities for the base
     * @param fourIntensity  the four raw intensities for the base
     */
    public void addMeanPoint(double[][] probMatrix, double[] fourIntensity) {
        for (int basePrevIndex = 0; basePrevIndex < numTheories; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                double weight = probMatrix[basePrevIndex][baseCurIndex];                

                DoubleMatrix1D weightedChannelIntensities = (DoubleFactory1D.dense).make(fourIntensity);
                weightedChannelIntensities.assign(F.mult(weight));

                sums[basePrevIndex][baseCurIndex].assign(weightedChannelIntensities, F.plus);
                counts[basePrevIndex][baseCurIndex] += weight;
            }
        }

        readyToCall = false;
    }

    /**
     * Add a single training point to the model to estimate the covariances.
     *
     * @param probMatrix     the matrix of probabilities for the base
     * @param fourIntensity  the four raw intensities for the base
     */
    public void addCovariancePoint(double[][] probMatrix, double[] fourIntensity) {
        for (int basePrevIndex = 0; basePrevIndex < numTheories; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                double weight = probMatrix[basePrevIndex][baseCurIndex];

                DoubleMatrix1D mean = sums[basePrevIndex][baseCurIndex].copy();
                mean.assign(F.div(counts[basePrevIndex][baseCurIndex]));

                DoubleMatrix1D sub = (DoubleFactory1D.dense).make(fourIntensity);
                sub.assign(mean, F.minus);

                DoubleMatrix2D cov = (DoubleFactory2D.dense).make(4, 4);
                alg.multOuter(sub, sub, cov);

                cov.assign(F.mult(weight));
                unscaledCovarianceSums[basePrevIndex][baseCurIndex].assign(cov, F.plus);
            }
        }
        
        readyToCall = false;
    }

    /**
     * Precompute all the matrix inversions and determinants we'll need for computing the likelihood distributions.
     */
    public void prepareToCallBases() {
        for (int basePrevIndex = 0; basePrevIndex < numTheories; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                means[basePrevIndex][baseCurIndex] = sums[basePrevIndex][baseCurIndex].copy();
                means[basePrevIndex][baseCurIndex].assign(F.div(counts[basePrevIndex][baseCurIndex]));

                inverseCovariances[basePrevIndex][baseCurIndex] = unscaledCovarianceSums[basePrevIndex][baseCurIndex].copy();
                inverseCovariances[basePrevIndex][baseCurIndex].assign(F.div(counts[basePrevIndex][baseCurIndex]));

                if (MathUtils.compareDoubles(alg.det(inverseCovariances[basePrevIndex][baseCurIndex]), 0.0) == 0) {
                    bustedCycle = true;
                    readyToCall = true;

                    return;
                }
                
                DoubleMatrix2D invcov = alg.inverse(inverseCovariances[basePrevIndex][baseCurIndex]);

                inverseCovariances[basePrevIndex][baseCurIndex] = invcov;

                norms[basePrevIndex][baseCurIndex] = Math.pow(alg.det(invcov), 0.5)/Math.pow(2.0*Math.PI, 2.0);
            }
        }

        readyToCall = true;
    }

    /**
     * Compute the likelihood matrix for a base.
     *
     * @param cycle         the cycle we're calling right now
     * @param fourintensity the four intensities of the current cycle's base
     * @return              a 4x4 matrix of likelihoods, where the row is the previous cycle base hypothesis and
     *                      the column is the current cycle base hypothesis
     */
    public double[][] computeLikelihoods(int cycle, double[] fourintensity) {
        if (!readyToCall) {
            prepareToCallBases();
        }

        double[][] likedist = new double[numTheories][4];

        if (bustedCycle) {
            likedist[0][0] = 1.0;
        } else {
            for (int basePrevIndex = 0; basePrevIndex < numTheories; basePrevIndex++) {
                for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                    double norm = norms[basePrevIndex][baseCurIndex];

                    DoubleMatrix1D sub = (DoubleFactory1D.dense).make(fourintensity);
                    sub.assign(means[basePrevIndex][baseCurIndex], F.minus);

                    DoubleMatrix1D Ax = alg.mult(inverseCovariances[basePrevIndex][baseCurIndex], sub);
                    double exparg = -0.5*alg.mult(sub, Ax);

                    likedist[basePrevIndex][baseCurIndex] = norm*Math.exp(exparg);
                }
            }
        }

        return likedist;
    }

    /**
     * Write the model parameters to disk.
     *
     * @param outparam  the file in which the output parameters should be stored
     */
    public void write(File outparam) {
        try {
            PrintWriter writer = new PrintWriter(outparam);

            for (int basePrevIndex = 0; basePrevIndex < numTheories; basePrevIndex++) {
                for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                    writer.print("mean_" + BaseUtils.baseIndexToSimpleBase(baseCurIndex) + " = c(");
                    for (int channel = 0; channel < 4; channel++) {
                        writer.print(sums[basePrevIndex][baseCurIndex].getQuick(channel)/counts[basePrevIndex][baseCurIndex]);

                        if (channel < 3) {
                            writer.print(", ");
                        }
                    }
                    writer.println(");");

                    DoubleMatrix2D cov = unscaledCovarianceSums[basePrevIndex][baseCurIndex].copy();
                    cov.assign(F.div(counts[basePrevIndex][baseCurIndex]));

                    writer.print("cov_" + BaseUtils.baseIndexToSimpleBase(baseCurIndex) + " = matrix(c(");
                    for (int channel1 = 0; channel1 < 4; channel1++) {
                        for (int channel2 = 0; channel2 < 4; channel2++) {
                            writer.print(cov.get(channel2, channel1) + (channel1 == 3 && channel2 == 3 ? "" : ","));
                        }
                    }
                    writer.println("), nr=4, nc=4);\n");
                }
            }

            writer.close();
        } catch (IOException e) {
        }
    }
}
