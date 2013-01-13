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

package org.broadinstitute.sting.gatk.walkers;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.MutableVariantContext;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.StandardVCFWriter;
import org.broad.tribble.vcf.VCFConstants;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.gatk.walkers.varianteval.MendelianViolationEvaluator;

import java.util.*;
import java.io.File;

/**
 * Implements an (STILL BEING TESTED) algorithm for calling SNPs in trios
 */

@By(DataSource.REFERENCE)
//@Requires(value={DataSource.REFERENCE, DataSource.REFERENCE_BASES, DataSource.READS},referenceMetaData={@RMD(name="sites",type= VariationRod.class)})
@Allows({DataSource.READS, DataSource.REFERENCE})
//, @RMD(name="parent1",type= VariationRod.class), @RMD(name="parent2",type= VariationRod.class)})
public class TrioGenotyperWalker extends RefWalker<VariantContext, Integer>{
    @Argument(shortName="mom", doc="", required=true)
    protected String mom;

    @Argument(shortName="dad", doc="", required=true)
    protected String dad;

    @Argument(shortName="kid", doc="", required=true)
    protected String kid;

    @Argument(shortName="log10PriorOfDeNovoOfTrueVariant", doc="", required=false)
    double LOG10_MENDEL_VIOLATION_PRIOR = -5;   // 30 in 3B bases

    @Argument(shortName = "varout", doc = "File to which variants should be written", required = true)
    public String vcfOutputFile = null;

    @ArgumentCollection
    private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    UnifiedGenotyperEngine UGEngine = null;
    private List<String> FAMILY_MEMBERS;
    private VCFWriter writer = null;

    public void initialize() {
        UGEngine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, null);
        // initialize the header
        FAMILY_MEMBERS = Arrays.asList(mom, dad, kid);

        // initialize the writer
        writer = new StandardVCFWriter(new File(vcfOutputFile));
    }

    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        VariantContext vc = tracker.getVariantContext(ref, "variants", EnumSet.of(VariantContext.Type.SNP), context.getLocation(), true);
        
        if ( vc != null && vc.isPolymorphic() ) {
            if ( ! vc.hasGenotypes(FAMILY_MEMBERS) )
                throw new StingException("variants file does not contain genotypes for everyone in family: " + FAMILY_MEMBERS);

            VariantCallContext call = UGEngine.runGenotyper(tracker, ref, context);

            // is call ever be null?
            vc = annotateTrioPrior(vc, call.vc);

            return vc;
        } else {
            return null;
        }
    }

    private VariantContext annotateTrioPrior(VariantContext vcIn, VariantContext call) {
        Genotype momG = call.getGenotype(mom);
        Genotype dadG = call.getGenotype(dad);
        Genotype kidG = call.getGenotype(kid);

        double log10POfGenotype = Double.MIN_VALUE;
        if ( MendelianViolationEvaluator.isViolation(call, momG, dadG, kidG) ) {
            Allele R = call.getReference();
            Allele A = call.getAlternateAllele(0);

            List<List<Allele>> possibleGenotypes = Arrays.asList(Arrays.asList(R,R), Arrays.asList(R,A), Arrays.asList(A,A));

            double[] L = new double[3 * 3 * 3];
            int i = 0, bestIndex = 0;
            double log10LOfBestGenotypes = Integer.MIN_VALUE;
            for ( List<Allele> momPG : possibleGenotypes ) {
                for ( List<Allele> dadPG : possibleGenotypes ) {
                    for ( List<Allele> kidPG : possibleGenotypes ) {
                        double log10LOfG = genotypeL(momPG, momG) + genotypeL(dadPG, dadG) + genotypeL(kidPG, kidG);
                        boolean isViolation = MendelianViolationEvaluator.isViolation(call, momPG, dadPG, kidPG);
                        double log10prior = isViolation ? LOG10_MENDEL_VIOLATION_PRIOR : 0;
                        L[i] = log10LOfG + log10prior;

                        if ( log10LOfG > log10LOfBestGenotypes ) {
                            bestIndex = i;
                            log10LOfBestGenotypes = log10LOfG;
                        }
                        logger.debug(String.format("%10s %10s => %10s : %b\t%.2f\t\t%.2f\t\t%3.2f", momPG, dadPG, kidPG, isViolation, log10LOfG, log10prior, L[i]));
                        i++;
                    }
                }
            }

            double[] posteriors = MathUtils.normalizeFromLog10(L, true);
            log10POfGenotype = posteriors[bestIndex];
        }
        //log10POfViolation = Math.min(log10POfViolation, 0);

        double Q = QualityUtils.phredScaleCorrectRate(Math.pow(10, log10POfGenotype));
        logger.debug(String.format("log10 P of best genotype log10 post = %.2f, Q = %.2f", log10POfGenotype, Q));
        MutableVariantContext mvc = new MutableVariantContext(vcIn);
        mvc.putAttribute("MVQ", Q);
        return new VariantContext(mvc);
    }

    /**
     * Isolate the rest of the walker from the code to get genotype likelihood values for allele A/B in genotypeCall
     * @param alleles
     * @param genotypeCall
     * @return
     */
    private double genotypeL( List<Allele> alleles, Genotype genotypeCall ) {
        String postTriplet = (String)genotypeCall.getAttribute(VCFConstants.GENOTYPE_LIKELIHOODS_KEY);
        if ( postTriplet == null )
            throw new StingException("BUG: TrioGenotyperWalker expected genotype likelihood triplets " + VCFConstants.GENOTYPE_LIKELIHOODS_KEY);

        // calculate the offset -- AA => 0, AB => 1, BB => 2
        int i = 0;
        for ( Allele a : alleles )
            i += a.isNonReference() ? 1 : 0;

        // convert the corresponding GL field to a double
        String log10LStrings[] = postTriplet.split(",");
        return Double.valueOf(log10LStrings[i]);
    }

    public Integer reduceInit() { return 0; }
    public Integer reduce(VariantContext vc, Integer a) {
        if ( vc != null ) {
            if ( a == 0 )
                writer.writeHeader(VariantContextAdaptors.createVCFHeader(null, vc));

            writer.add(vc, (byte)'.');
            a++;
        }

        return a;
    }

    public void onTraversalDone(Integer result) {} // Don't print the reduce result
}