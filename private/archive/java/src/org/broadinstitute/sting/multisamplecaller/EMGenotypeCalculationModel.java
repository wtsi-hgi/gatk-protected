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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;

public abstract class EMGenotypeCalculationModel extends GenotypeCalculationModel {

    // We need to set a limit on the EM iterations in case something flukey goes on
    protected static final int MAX_EM_ITERATIONS = 8;

    // We consider the EM stable when the MAF doesn't change more than 1/10,000
    protected static final double EM_STABILITY_METRIC = 1e-4;

    protected EMGenotypeCalculationModel() {}

    public VariantCallContext callLocus(RefMetaDataTracker tracker, char ref, GenomeLoc loc, Map<String, StratifiedAlignmentContext> contexts, DiploidGenotypePriors priors) {

        // run the EM calculation
        EMOutput overall = runEM(ref, contexts, priors, StratifiedAlignmentContext.StratifiedContextType.COMPLETE);

        double PofD = Math.pow(10, overall.getPofD());
        double PofNull = Math.pow(10, overall.getPofNull());
        double sum = PofD + PofNull;

        // calculate the phred-scaled confidence score
        double phredScaledConfidence;
        boolean bestIsRef = false;
        if ( PofD > PofNull ) {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(PofNull / sum);
        } else {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(PofD / sum);
            bestIsRef = true;
        }

        // are we above the lod threshold for emitting calls (and not in all-bases mode)?
        if ( !ALL_BASE_MODE && ((!GENOTYPE_MODE && bestIsRef) || phredScaledConfidence < CONFIDENCE_THRESHOLD) )
            return new VariantCallContext(phredScaledConfidence >= CONFIDENCE_THRESHOLD);

        // generate the calls
        List<Genotype> calls = genotypeCallsFromGenotypeLikelihoods(overall, ref, contexts);

        VariationCall locusdata = GenotypeWriterFactory.createSupportedCall(OUTPUT_FORMAT, ref, loc, bestIsRef ? Variation.VARIANT_TYPE.REFERENCE : Variation.VARIANT_TYPE.SNP);
        if ( locusdata != null ) {
            if ( locusdata instanceof ConfidenceBacked ) {
                ((ConfidenceBacked)locusdata).setConfidence(phredScaledConfidence);
            }
            if ( locusdata instanceof IDBacked ) {
                rodDbSNP dbsnp = getDbSNP(tracker);
                if ( dbsnp != null )
                    ((IDBacked)locusdata).setID(dbsnp.getRS_ID());
            }
            if ( locusdata instanceof SLODBacked && REPORT_SLOD ) {

                // calculate strand score

                double lod = overall.getPofD() - overall.getPofNull();
                logger.debug(String.format("LOD=%f", lod));

                EMOutput forward = runEM(ref, contexts, priors, StratifiedAlignmentContext.StratifiedContextType.FORWARD);
                EMOutput reverse = runEM(ref, contexts, priors, StratifiedAlignmentContext.StratifiedContextType.REVERSE);
                double forwardLod = (forward.getPofD() + reverse.getPofNull()) - overall.getPofNull();
                double reverseLod = (reverse.getPofD() + forward.getPofNull()) - overall.getPofNull();
                logger.debug("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);
                double strandScore = Math.max(forwardLod - lod, reverseLod - lod);
                logger.debug(String.format("SLOD=%f", strandScore));
                // rescale by a factor of 10
                strandScore *= 10.0;

                ((SLODBacked)locusdata).setSLOD(strandScore);
            }
            locusdata.setNonRefAlleleFrequency(overall.getMAF());

            // finally, associate the Variation with the Genotypes
            locusdata.setGenotypeCalls(calls);
            for ( Genotype call : calls )
                ((GenotypeCall)call).setVariation(locusdata);
        }
        return new VariantCallContext(locusdata, calls, phredScaledConfidence >= CONFIDENCE_THRESHOLD);
    }

    protected List<Genotype> genotypeCallsFromGenotypeLikelihoods(EMOutput results, char ref, Map<String, StratifiedAlignmentContext> contexts) {
        HashMap<String, GenotypeLikelihoods> GLs = results.getGenotypeLikelihoods();

        ArrayList<Genotype> calls = new ArrayList<Genotype>();
        int variantCalls = 0;

        for ( String sample : GLs.keySet() ) {

            // create the call
            AlignmentContext context = contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE);
            GenotypeCall call = GenotypeWriterFactory.createSupportedGenotypeCall(OUTPUT_FORMAT, ref, context.getLocation());

            // set the genotype and confidence
            double[] posteriors = GLs.get(sample).getPosteriors();
            Integer sorted[] = Utils.SortPermutation(posteriors);
            DiploidGenotype bestGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 1]];
            DiploidGenotype nextGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 2]];
            call.setNegLog10PError(posteriors[bestGenotype.ordinal()] - posteriors[nextGenotype.ordinal()]);
            call.setGenotype(bestGenotype);

            if ( call instanceof ReadBacked ) {
                ReadBackedPileup pileup = contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
                ((ReadBacked)call).setPileup(pileup);
            }
            if ( call instanceof SampleBacked ) {
                ((SampleBacked)call).setSampleName(sample);
            }
            if ( call instanceof LikelihoodsBacked ) {
                ((LikelihoodsBacked)call).setLikelihoods(GLs.get(sample).getLikelihoods());
            }
            if ( call instanceof PosteriorsBacked ) {
                ((PosteriorsBacked)call).setPosteriors(posteriors);
            }

            calls.add(call);

            // increment the variant count if it's non-ref
            if ( call.isVariant(ref) )
                variantCalls++;
        }

        // if everyone is ref, don't emit any calls
        if ( variantCalls == 0 )
            calls.clear();

        return calls;
    }

    public EMOutput runEM(char ref, Map<String, StratifiedAlignmentContext> contexts, DiploidGenotypePriors priors, StratifiedAlignmentContext.StratifiedContextType contextType) {

        // initialize the allele frequencies
        initializeAlleleFrequencies(contexts.size(), ref);

        // initialize the genotype likelihoods
        initializeGenotypeLikelihoods(ref, contexts, priors, contextType);

        // The EM loop:
        //   we want to continue until the calculation is stable, but we need some max on the number of iterations
        int iterations = 0;
        boolean EM_IS_STABLE;

        do {
            calculateAlleleFrequencyPosteriors();

            applyAlleleFrequencyToGenotypeLikelihoods();

            EM_IS_STABLE = isStable();

        } while ( ++iterations < MAX_EM_ITERATIONS && !EM_IS_STABLE );

        logger.debug("EM loop took " + iterations + " iterations");

        return computePofF(ref);
    }

    protected abstract void initializeAlleleFrequencies(int numSamplesInContext, char ref);
    protected abstract void initializeGenotypeLikelihoods(char ref, Map<String, StratifiedAlignmentContext> contexts, DiploidGenotypePriors priors, StratifiedAlignmentContext.StratifiedContextType contextType);
    protected abstract void calculateAlleleFrequencyPosteriors();
    protected abstract void applyAlleleFrequencyToGenotypeLikelihoods();
    protected abstract boolean isStable();
    protected abstract EMOutput computePofF(char ref);


    /**
     * A class to keep track of the EM output
     */
    protected class EMOutput {
        private double pD, pNull, pF, MAF;
        private HashMap<String, GenotypeLikelihoods> GLs;

        EMOutput(double pD, double pNull, double pF, double MAF, HashMap<String, GenotypeLikelihoods> GLs) {
            this.pD = pD;
            this.pNull = pNull;
            this.pF = pF;
            this.MAF = MAF;
            this.GLs = GLs;
        }

        public double getPofD() { return pD; }
        public double getPofNull() { return pNull; }
        public double getPofF() { return pF; }
        public double getMAF() { return MAF; }
        public HashMap<String, GenotypeLikelihoods> getGenotypeLikelihoods() { return GLs; }
    }
}