/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.siblingibd;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.MathUtils;

import java.util.*;

/**
 * Created by cwhelan on 5/19/14.
 *
 * Models IBD state across the genomes of two siblings using an HMM.
 *
 */
public class QuartetIBDStateHMM {
    final protected static Logger logger = Logger.getLogger(QuartetIBDStateHMM.class);

    private final List<IBDObservation> observations = new ArrayList<>();
    private final List<VariantContext> sites = new ArrayList<>();

    // these constants tuned to be arbitrarily low based on experiments conducted on
    // a human whole genome data set; we found that the system performs best when
    // it is discouraged from switching from one state to another.
    private static final double DEFAULT_STATE_CHANGE_PROB = 1.0 / 10000000000000000.0;
    private static final double DEFAULT_ERROR_PROB = 1.0 / 1000000000000000000.0;

    // emission probabilities for the HMM; we found that uniform emission probabilities
    // across valid or invalid observations performed well on a human test data set,
    // with a slight tweak to give extra probability to IBD_ZERO_OR_TWO observations in
    // IBD1 states - a common error mode due to overcalling of het genotypes
    public static final double NON_ERROR_EMISSION_PROBABILITY = .97;
    public static final double ERROR_EMISSION_PROBABILITY = .03;
    public static final double IBD1_HET_ERROR_EMISSION_PROBABILITY = .02;
    public static final double IDB1_OTHER_ERROR_EMISSION_PROBABILITY = .01;

    private String chrom = null;

    private EnumMap<IBDState, Double> initialProbs;
    private EnumMap<IBDState, Double[]> transitionProbs;

    public QuartetIBDStateHMM() {
        initProbs(DEFAULT_STATE_CHANGE_PROB, DEFAULT_ERROR_PROB);
    }

    public QuartetIBDStateHMM(final double recombinationProb, final double errorProb) {
        initProbs(recombinationProb, errorProb);
    }

    /**
     * Initialize initial and transition probabilities
     *
     * HMM probabilities are modeled simplistically, with uniform initial probabilities and
     * uniform probabilities for switching states
     * @param recombinationProb
     * @param errorProb
     */
    private void initProbs(final double recombinationProb, final double errorProb) {

        initialProbs = new EnumMap<>(IBDState.class);
        initialProbs.put(IBDState.ZERO, Math.log10(.33));
        initialProbs.put(IBDState.ONE, Math.log10(.33));
        initialProbs.put(IBDState.TWO, Math.log10(.33));

        transitionProbs = new EnumMap<>(IBDState.class);
        transitionProbs.put(IBDState.ZERO, new Double[]{Math.log10((1 - recombinationProb) - (errorProb / 2)), Math.log10(recombinationProb - (errorProb / 2)), Math.log10(errorProb)});
        transitionProbs.put(IBDState.ONE, new Double[]{Math.log10((recombinationProb / 2)), Math.log10(1 - recombinationProb), Math.log10(recombinationProb / 2)});
        transitionProbs.put(IBDState.TWO, new Double[]{Math.log10(errorProb), Math.log10(recombinationProb - (errorProb / 2)), Math.log10((1 - recombinationProb) - (errorProb / 2))});
    }

    /**
     * Add an observation to the HMM. If we are ready to run the model (ie if we've changed chromosome, execute Viterbi
     * and return an iterator over the predicted IBD state for each locus for which we've seen an observation. Otherwise,
     * return an empty iterator.
     *
     * @param vc
     * @param observation
     * @return
     */
    public Iterator<IBDLocus> addObservation(final VariantContext vc, final IBDObservation observation) {
        if (this.chrom != null && !vc.getChr().equals(this.chrom)) {
            final List<IBDLocus> results = runModelAndCoallateResults();
            this.chrom = vc.getChr();
            observations.clear();
            sites.clear();
            observations.add(observation);
            sites.add(vc);

            return results.iterator();
        } else {
            observations.add(observation);
            sites.add(vc);
            this.chrom = vc.getChr();
            return Collections.<IBDLocus>emptyList().iterator();
        }
    }

    /**
     * Execute Viterbi for the final chromosome and return the results
     * @return
     */
    public Iterator<IBDLocus> runLastChrom() {
        return runModelAndCoallateResults().iterator();
    }

    private List<IBDLocus> runModelAndCoallateResults() {
        final ModelResults modelResults = runModel();
        final IBDState[] states = modelResults.viterbiStates;
        final List<IBDLocus> results = new ArrayList<>(observations.size());
        for (int i = 0; i < states.length; i++) {
            results.add(new IBDLocus(sites.get(i), states[i], modelResults.statePosteriors[i]));
        }
        return results;
    }

    protected ModelResults runModel() {
        final int numObservations = observations.size();
        if (numObservations == 0) {
            final ModelResults nullResult = new ModelResults();
            nullResult.viterbiStates = new IBDState[0];
            nullResult.statePosteriors = new double[0][0];
            return nullResult;
        }

        final double[][] viterbi = new double[numObservations][IBDState.values().length];

        // forward, backward, and posterior probabilities for forward-backward algorithm
        final double[][] alphas = new double[numObservations][IBDState.values().length];
        final double[][] betas = new double[numObservations][IBDState.values().length];
        final double[][] gammas = new double[numObservations][IBDState.values().length];

        final short[][] path = new short[numObservations][IBDState.values().length];

        // forward and viterbi initialization
        for (final IBDState i : IBDState.values()) {
            viterbi[0][i.ordinal()] = initialProbs.get(i) + getEmission(i, observations.get(0));
            alphas[0][i.ordinal()] = initialProbs.get(i) + getEmission(i, observations.get(0));
        }

        // forward and viterbi
        for (int t = 1; t < numObservations; t++) {
            for (final IBDState j : IBDState.values()) {

                double maxProb = Double.NEGATIVE_INFINITY;
                int argMaxProb = -1;

                double sumProb = Double.NEGATIVE_INFINITY;

                for (final IBDState i : IBDState.values()) {
                    final double prob = viterbi[t - 1][i.ordinal()] + transitionProbs.get(i)[j.ordinal()] + getEmission(j, observations.get(t));
                    if (prob > maxProb) {
                        maxProb = prob;
                        argMaxProb = i.ordinal();
                    }

                    sumProb = MathUtils.approximateLog10SumLog10(sumProb, alphas[t - 1][i.ordinal()] + transitionProbs.get(i)[j.ordinal()] + getEmission(j, observations.get(t)));
                }
                viterbi[t][j.ordinal()] = maxProb;
                path[t][j.ordinal()] = (short) argMaxProb;

                alphas[t][j.ordinal()] = sumProb;
            }
        }

        // backward initialization
        for (final IBDState s : IBDState.values()) {
            betas[numObservations-1][s.ordinal()] = 0;
        }

        // backward
        for (int t = numObservations - 2; t >= 0; t--) {
            for (final IBDState i : IBDState.values()) {
                double sumProb = Double.NEGATIVE_INFINITY;
                for (final IBDState j : IBDState.values()) {
                    sumProb = MathUtils.approximateLog10SumLog10(sumProb, betas[t + 1][j.ordinal()] + transitionProbs.get(i)[j.ordinal()] + getEmission(j, observations.get(t + 1)));
                }
                betas[t][i.ordinal()] = sumProb;
            }
        }

        // backtrack through viterbi to find most likely states
        final IBDState[] states = new IBDState[numObservations];
        double maxProb = Double.NEGATIVE_INFINITY;
        int argMaxProb = -1;
        for (final IBDState s : IBDState.values()) {
            if (viterbi[numObservations-1][s.ordinal()] > maxProb) {
                maxProb = viterbi[numObservations-1][s.ordinal()];
                argMaxProb = s.ordinal();
            }
        }
        states[numObservations - 1] = IBDState.values()[argMaxProb];
        for (int i = numObservations - 1; i >= 1; i--) {
            states[i-1] = IBDState.values()[path[i][states[i].ordinal()]];
        }

        // gammas give posteriors of being in each state
        for (int t = 0; t < numObservations; t++) {
            double totalProb = Double.NEGATIVE_INFINITY;
            for (final IBDState i : IBDState.values()) {
                final double gamma = alphas[t][i.ordinal()] + betas[t][i.ordinal()];
                totalProb = MathUtils.approximateLog10SumLog10(totalProb, gamma);
                gammas[t][i.ordinal()] = gamma;
            }

            // normalize
            for (final IBDState j : IBDState.values()) {
                gammas[t][j.ordinal()] = gammas[t][j.ordinal()] - totalProb;
            }
        }

        final ModelResults modelResults = new ModelResults();
        modelResults.viterbiStates = states;
        modelResults.statePosteriors = gammas;
        return modelResults;
    }


    /**
     * Emission probabilities are simplistically modeled also, with (almost) uniform probabilities
     * for each allowable IBD observation in a given state and for error observations. A small amount
     * of extra probability is given to ZERO_OR_TWO observations (which arise when the quartet is all Het)
     * in IBD1 regions, which seem to be the most frequent error mode, likely arising from segmental duplications
     * or high-copy number regions
     *
     * @param s
     * @param o
     * @return
     */
    private Double getEmission(final IBDState s, final IBDObservation o) {
        switch (s) {
            case ZERO:
                if (o == IBDObservation.ZERO || o == IBDObservation.ZERO_OR_ONE || o == IBDObservation.ZERO_OR_TWO) {
                    return Math.log10(NON_ERROR_EMISSION_PROBABILITY / 3);
                } else {
                    return Math.log10(ERROR_EMISSION_PROBABILITY / 3);
                }
            case ONE:
                if (o == IBDObservation.ONE || o == IBDObservation.ZERO_OR_ONE || o == IBDObservation.ONE_OR_TWO) {
                    return Math.log10(NON_ERROR_EMISSION_PROBABILITY / 3);
                } else if (o == IBDObservation.ZERO_OR_TWO) {
                    return Math.log10(IBD1_HET_ERROR_EMISSION_PROBABILITY);
                } else {
                    return Math.log10(IDB1_OTHER_ERROR_EMISSION_PROBABILITY / 2);
                }
            case TWO:
                if (o == IBDObservation.TWO || o == IBDObservation.ZERO_OR_TWO || o == IBDObservation.ONE_OR_TWO) {
                    return Math.log10(NON_ERROR_EMISSION_PROBABILITY / 3);
                } else {
                    return Math.log10(ERROR_EMISSION_PROBABILITY / 3);
                }
        }
        return Double.NEGATIVE_INFINITY;
    }

    protected class ModelResults {
        IBDState[] viterbiStates;
        double[][] statePosteriors;
    }
}
