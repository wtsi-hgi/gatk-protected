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
import org.broadinstitute.variant.utils.BaseUtils;

import java.util.*;

public class PointEstimateGenotypeCalculationModel extends EMGenotypeCalculationModel {

    protected PointEstimateGenotypeCalculationModel() {}

    // the allele frequencies
    private double[] alleleFrequencies = new double[4];
    private double[] oldAlleleFrequencies;

    // the GenotypeLikelihoods map
    private HashMap<String, GenotypeLikelihoods> GLs = new HashMap<String, GenotypeLikelihoods>();

    // Allele frequency initialization values from the original MSG code (so we can be consistent)
    private static final double NON_REF = 0.0005002502;  // heterozygosity / (2 * sqrt(1-heterozygosity)
    private static final double REF = 0.9994999;         //sqrt(1-heterozygosity)


    // overload this method so we can special-case the single sample
    public VariantCallContext callLocus(RefMetaDataTracker tracker, char ref, GenomeLoc loc, Map<String, StratifiedAlignmentContext> contexts, DiploidGenotypePriors priors) {

        // we don't actually want to run EM for single samples
        if ( samples.size() == 1 ) {

            // get the context for the sample
            String sample = samples.iterator().next();
            StratifiedAlignmentContext sampleContext = contexts.get(sample);

            // if there were no good bases, the context wouldn't exist
            if ( sampleContext == null )
                return null;

            // get the genotype likelihoods
            Pair<ReadBackedPileup, GenotypeLikelihoods> discoveryGL = getSingleSampleLikelihoods(sampleContext, priors, StratifiedAlignmentContext.StratifiedContextType.COMPLETE);

            // find the index of the best genotype
            double[] normPosteriors = discoveryGL.second.getNormalizedPosteriors();
            Integer sortedNormPosteriors[] = Utils.SortPermutation(normPosteriors);
            int bestIndex = sortedNormPosteriors[sortedNormPosteriors.length - 1];

            // flag to determine if ref is the best call (not necessary in genotype mode)
            boolean bestIsRef = false;

            // calculate the phred-scaled confidence score
            double phredScaledConfidence;
            if ( GENOTYPE_MODE ) {
                phredScaledConfidence = QualityUtils.phredScaleErrorRate(1.0 - normPosteriors[bestIndex]);
            } else {
                int refIndex = DiploidGenotype.createHomGenotype(ref).ordinal();
                bestIsRef = (refIndex == bestIndex);
                double pError = (bestIsRef ? 1.0 - normPosteriors[refIndex] : normPosteriors[refIndex]);
                phredScaledConfidence = QualityUtils.phredScaleErrorRate(pError);
            }

            // are we above the lod threshold for emitting calls (and not in all-bases mode)?
            if ( !ALL_BASE_MODE && ((!GENOTYPE_MODE && bestIsRef) || phredScaledConfidence < CONFIDENCE_THRESHOLD) )
                return new VariantCallContext(phredScaledConfidence >= CONFIDENCE_THRESHOLD);

            // we can now create the genotype call object
            GenotypeCall call = GenotypeWriterFactory.createSupportedGenotypeCall(OUTPUT_FORMAT, ref, loc);

            // set the genotype and confidence
            double[] posteriors = discoveryGL.second.getPosteriors();
            Integer sorted[] = Utils.SortPermutation(posteriors);
            DiploidGenotype bestGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 1]];
            DiploidGenotype nextGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 2]];
            call.setNegLog10PError(posteriors[bestGenotype.ordinal()] - posteriors[nextGenotype.ordinal()]);
            call.setGenotype(bestGenotype);

            if ( call instanceof ReadBacked ) {
                ((ReadBacked)call).setPileup(discoveryGL.first);
            }
            if ( call instanceof SampleBacked ) {
                ((SampleBacked)call).setSampleName(sample);
            }
            if ( call instanceof LikelihoodsBacked ) {
                ((LikelihoodsBacked)call).setLikelihoods(discoveryGL.second.getLikelihoods());
            }
            if ( call instanceof PosteriorsBacked ) {
                ((PosteriorsBacked)call).setPosteriors(posteriors);
            }

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
                locusdata.setGenotypeCalls(Arrays.asList((Genotype)call));
            }
            
            call.setVariation(locusdata);

            return new VariantCallContext(locusdata, Arrays.asList((Genotype)call), phredScaledConfidence >= CONFIDENCE_THRESHOLD);
        }

        return super.callLocus(tracker, ref, loc, contexts, priors);
    }

    private Pair<ReadBackedPileup, GenotypeLikelihoods> getSingleSampleLikelihoods(StratifiedAlignmentContext sampleContext, DiploidGenotypePriors priors, StratifiedAlignmentContext.StratifiedContextType contextType) {
        // create the pileup
        AlignmentContext myContext = sampleContext.getContext(contextType);
        ReadBackedPileup pileup = myContext.getBasePileup();

        // create the GenotypeLikelihoods object
        GenotypeLikelihoods GL = new GenotypeLikelihoods(baseModel, priors, defaultPlatform);
        GL.add(pileup, true);
        return new Pair<ReadBackedPileup, GenotypeLikelihoods>(pileup, GL);
    }

    protected void initializeAlleleFrequencies(int numSamplesInContext, char ref) {
        for (int i = 0; i < 4; i++)
            alleleFrequencies[i] = NON_REF;
        alleleFrequencies[BaseUtils.simpleBaseToBaseIndex(ref)] = REF;

        for (int i = 0; i < 4; i++)
            logger.debug("Initial allele frequency for " + BaseUtils.baseIndexToSimpleBase(i) + ": " + alleleFrequencies[i]);
    }

    protected void initializeGenotypeLikelihoods(char ref, Map<String, StratifiedAlignmentContext> contexts, DiploidGenotypePriors priors, StratifiedAlignmentContext.StratifiedContextType contextType) {
        GLs.clear();

        DiploidGenotypePriors AFPriors = calculateAlleleFreqBasedPriors(alleleFrequencies);

        for ( String sample : contexts.keySet() ) {
            StratifiedAlignmentContext context = contexts.get(sample);
            ReadBackedPileup pileup = context.getContext(contextType).getBasePileup();

            // create the GenotypeLikelihoods object
            GenotypeLikelihoods GL = new GenotypeLikelihoods(baseModel, AFPriors, defaultPlatform);
            GL.add(pileup, true);

            GLs.put(sample, GL);
        }
    }

    private static DiploidGenotypePriors calculateAlleleFreqBasedPriors(double[] alleleFreqs) {
        // convert to log-space
        double[] log10Freqs = new double[4];
        for (int i = 0; i < 4; i++)
            log10Freqs[i] = Math.log10(alleleFreqs[i]);

        double[] alleleFreqPriors = new double[10];

        // this is the Hardy-Weinberg based allele frequency (p^2, q^2, 2pq)
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            alleleFreqPriors[g.ordinal()] = log10Freqs[BaseUtils.simpleBaseToBaseIndex(g.base1)] + log10Freqs[BaseUtils.simpleBaseToBaseIndex(g.base2)];
            // add a factor of 2 for the 2pq case
            if ( g.isHet() )
                alleleFreqPriors[g.ordinal()] += Math.log10(2);
        }

        return new DiploidGenotypePriors(alleleFreqPriors);
    }

    protected void calculateAlleleFrequencyPosteriors() {
        // initialization
        oldAlleleFrequencies = alleleFrequencies.clone();
        for (int i = 0; i < 4; i++)
            alleleFrequencies[i] = 0.0;

        for ( GenotypeLikelihoods GL : GLs.values() ) {
            double[] normalizedPosteriors = GL.getNormalizedPosteriors();

            // calculate the posterior weighted frequencies for this sample
            double[] personalAllelePosteriors = new double[4];
            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                double posterior = normalizedPosteriors[g.ordinal()] / 2.0;   // each base gets half the probability
                personalAllelePosteriors[BaseUtils.simpleBaseToBaseIndex(g.base1)] += posterior;
                personalAllelePosteriors[BaseUtils.simpleBaseToBaseIndex(g.base2)] += posterior;
            }

            for (int i = 0; i < 4; i++)
                alleleFrequencies[i] += personalAllelePosteriors[i];
        }

        // normalize
        double sum = 0.0;
        for (int i = 0; i < 4; i++)
            sum += alleleFrequencies[i];
        for (int i = 0; i < 4; i++)
            alleleFrequencies[i] /= sum;

        for (int i = 0; i < 4; i++)
            logger.debug("New allele frequency for " + BaseUtils.baseIndexToSimpleBase(i) + ": " + alleleFrequencies[i]);
    }

    protected void applyAlleleFrequencyToGenotypeLikelihoods() {
        DiploidGenotypePriors AFPriors = calculateAlleleFreqBasedPriors(alleleFrequencies);
        for ( GenotypeLikelihoods GL : GLs.values() )
            GL.setPriors(AFPriors);
    }

    protected boolean isStable() {
        // We consider the EM stable when the MAF doesn't change more than EM_STABILITY_METRIC
        double AF_delta = 0.0;
        for (int i = 0; i < 4; i++)
            AF_delta += Math.abs(oldAlleleFrequencies[i] - alleleFrequencies[i]);

        return (AF_delta < EM_STABILITY_METRIC);
    }

    protected EMOutput computePofF(char ref) {
        // some debugging output
        for ( String sample : GLs.keySet() )
            logger.debug("GenotypeLikelihoods for sample " + sample + ": " + GLs.get(sample).toString());

        // compute pD and pNull without allele frequencies
        double pD = compute_pD(GLs);
        double pNull = compute_pNull(ref, GLs);
        logger.debug("Original pD=" + pD + ", pNull=" + pNull);

        // compute p0
        double pVar = 0.0;
        for (int i = 1; i < GLs.size(); i++)
            pVar += heterozygosity/(double)i;
        double p0 = Math.log10(1.0 - pVar);

        // compute actual priors: theta / MAF
        double MAF;
        Integer[] sortedIndexes = Utils.SortPermutation(alleleFrequencies);
        if ( sortedIndexes[3] != BaseUtils.simpleBaseToBaseIndex(ref) )
            MAF = alleleFrequencies[sortedIndexes[3]];
        else
            MAF = alleleFrequencies[sortedIndexes[2]];

        //  compute pF
        double pF;
        double expectedChromosomes = 2.0 * (double)GLs.size() * MAF;
        if ( expectedChromosomes < 1.0 )
            pF = p0;
        else
            pF = Math.log10(heterozygosity / expectedChromosomes);
        logger.debug("p0=" + p0 + ", pF=" + pF);

        pD += pF;
        pNull += p0;
        logger.debug("Final pD=" + pD + ", pNull=" + pNull);

        return new EMOutput(pD, pNull, pF, MAF, GLs);
    }

    private static double compute_pD(HashMap<String, GenotypeLikelihoods> GLs) {
        double pD = 0.0;
        for ( GenotypeLikelihoods GL : GLs.values() ) {
            double sum = 0.0;
            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                sum += Math.pow(10, GL.getPosterior(g));
            }
            pD += Math.log10(sum);
        }
        return pD;
    }

    private static double compute_pNull(char ref, HashMap<String, GenotypeLikelihoods> GLs) {
        // compute null likelihoods
        double[] alleleLikelihoods = new double[4];
        for (int i = 0; i < 4; i++)
            alleleLikelihoods[i] = 1e-6/3.0;
        alleleLikelihoods[BaseUtils.simpleBaseToBaseIndex(ref)] = 1.0-1e-6;
        DiploidGenotypePriors AFPriors = calculateAlleleFreqBasedPriors(alleleLikelihoods);

        HashMap<String, GenotypeLikelihoods> GL_null = new HashMap<String, GenotypeLikelihoods>();
        try {
            for ( String sample : GLs.keySet() ) {
                GenotypeLikelihoods GL = (GenotypeLikelihoods)GLs.get(sample).clone();
                GL.setPriors(AFPriors);
                GL_null.put(sample, GL);
            }
        } catch (CloneNotSupportedException e) {
            throw new StingException("Clone() not supported for given GenotypeLikelihoods subclass?");
        }

        return compute_pD(GL_null);
    }
}