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

import org.apache.commons.lang.NotImplementedException;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.variant.utils.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.variantcontext.*;

import java.util.*;

import org.broadinstitute.variant.vcf.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 8/30/11
 * Time: 10:08 AM
 * To change this template use File | Settings | File Templates.
 */
public class UGBoundAF extends RodWalker<VariantContext,Integer> {

    @Output(shortName="vcf",fullName="VCF",doc="file to write to",required=true)
    VCFWriter writer;

    @Input(shortName="V",fullName="Variants",doc="variant tracks to use in calculation",required=true)
    List<RodBinding<VariantContext>> variants;

    private static double EPS_LOWER_LIMIT = Math.pow(10,-6.0);

    private HashMap<Integer,Pair<Double,Double>> epsilonPosteriorCache = new HashMap<Integer,Pair<Double,Double>>(8192);
    private HashMap<Integer,Double> logAC0Cache = new HashMap<Integer, Double>(8192);
    private int QUANTIZATION_FACTOR = 1000;


    public void initialize() {
        Set<VCFHeaderLine> allHeaderLines = new HashSet<VCFHeaderLine>(1024);
        for ( RodBinding<VariantContext> v : variants ) {
            String trackName = v.getName();
            Map<String, VCFHeader> vcfHeaders = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(trackName));
            Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>(vcfHeaders.get(trackName).getMetaData());
        }
        allHeaderLines.add(new VCFInfoHeaderLine("AFB",2,VCFHeaderLineType.Float,"The 95% bounds on the allele "+
                "frequency. First value=95% probability AF>x. Second value=95% probability AF<x."));
        writer.writeHeader(new VCFHeader(allHeaderLines));
    }
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext unused ) {
        List<VariantContext> allVariants = tracker.getValues(variants);
        if ( allVariants.size() == 0 ) {
            return null;
        }

        List<Allele> alternateAlleles = getAllAlternateAlleles(allVariants);
        VariantContextBuilder builder = new VariantContextBuilder(allVariants.get(0).subContextFromSamples(new TreeSet<String>()));
        if ( alternateAlleles.size() > 1 ) {
            logger.warn("Multiple Segregating Variants at position "+ref.getLocus().toString());
            alternateAlleles.add(allVariants.get(0).getReference());
            builder.alleles(alternateAlleles);
            builder.filters(String.format("MULTIPLE_SEGREGATING[%s]", Utils.join(",",alternateAlleles)));
        } else {
            // get all the genotype likelihoods
            GenotypesContext context = GenotypesContext.create();
            int numNoCall = 0;
            for ( VariantContext v : allVariants ) {
                numNoCall += v.getNoCallCount();
                context.addAll(v.getGenotypes());
            }
            builder.attribute("AFB",boundAlleleFrequency(getACPosteriors(context)));
        }

        return builder.make();
    }

    private List<Allele> getAllAlternateAlleles(List<VariantContext> variants) {
        List<Allele> alleles = new ArrayList<Allele>(3); // some overhead
        for ( VariantContext v : variants ) {
            alleles.addAll(v.getAlternateAlleles());
        }
        return alleles;
    }

    @Override
    public Integer reduce(VariantContext value, Integer sum) {
        if ( value == null )
            return sum;
        writer.add(value);
        return ++sum;
    }

    private int N_ITERATIONS = 1;
    private double[] getACPosteriors(GenotypesContext gc) {
        // note this uses uniform priors (!)

        double[][] zeroPriors = new double[1][1+2*gc.size()];
        AlleleFrequencyCalculationResult result = new AlleleFrequencyCalculationResult(2,2*gc.size());
        // todo -- allow multiple alleles here
        for ( int i = 0; i < N_ITERATIONS; i ++ ) {
            ExactAFCalculationModel.linearExactMultiAllelic(gc, 2, zeroPriors, result, false);
        }

        return result.log10AlleleFrequencyPosteriors[0];
    }

    private String boundAlleleFrequency(double[] ACLikelihoods) {
        // note that no-calls are unnecessary: the ML likelihoods take nocalls into account as 0,0,0 GLs
        // thus, for sites with K 100,40,0 likelihoods and M no-calls, the likelihoods will be
        // agnostic between 2*K alleles through 2*(K+M) alleles - exactly what we want to marginalize over

        // want to pick a lower limit x and upper limit y such that
        // int_{f = x to y} sum_{c = 0 to 2*AN} P(AF=f | c, AN) df = 0.95
        // int_{f=x to y} calculateAFPosterior(f) df = 0.95
        // and that (y-x) is minimized

        // this is done by quantizing [0,1] into small bins and, since the distribution is
        // unimodal, greedily adding them until the probability is >= 0.95

        throw new ReviewedStingException("This walker is unsupported, and is not fully implemented", new NotImplementedException("bound allele frequency not implemented"));
    }

    private double calculateAFPosterior(double[] likelihoods, double af) {
        double[] probLiks = new double[likelihoods.length];
        for ( int c = 0; c < likelihoods.length; c++) {
            probLiks[c] = calculateAFPosterior(c,likelihoods.length,af);
        }

        return MathUtils.log10sumLog10(probLiks);
    }

    private double calculateAFPosterior(int ac, int n, double af) {
        // evaluate the allele frequency posterior distribution at AF given AC observations of N chromosomes
        switch ( ac ) {
            case 0:
                return logAC0Coef(n) + n*Math.log10(1 - af) - Math.log10(af);
            case 1:
                return Math.log10(n) + (n-1)*Math.log10(1-af) - n*Math.log10(1-EPS_LOWER_LIMIT);
            case 2:
                return Math.log10(n) + Math.log10(n-1) + Math.log10(af) + (n-2)*Math.log10(1-af) - Math.log10(1-(n-1)*EPS_LOWER_LIMIT) - (n-1)*Math.log10(EPS_LOWER_LIMIT);
            default:
                return  (ac-1)*Math.log10(af)+ac*Math.log10( (double) n-ac)-(n-ac)*af*Math.log10(Math.E) - MathUtils.log10Gamma(ac);
        }
    }

    private double logAC0Coef(int an) {
        if ( ! logAC0Cache.containsKey(an) ) {
            double coef = -Math.log10(EPS_LOWER_LIMIT);
            for ( int k = 1; k <= an; k++ ) {
                // note this should typically just be
                // term = ( 1 - Math.pow(EPS_LOWER_LIMIT,k) ) * MathUtils.binomialCoefficient(an,k) / k
                // but the 1-E term will just be 1, so we do the following to mitigate this problem
                double binom = MathUtils.binomialCoefficient(an,k);
                double eps_correction = EPS_LOWER_LIMIT*Math.pow(binom,1/k);
                double term = binom/k - Math.pow(eps_correction,k);
                if ( k % 2 == 0 ) {
                    coef += term;
                }  else {
                    coef -= term;
                }
            }

            logAC0Cache.put(an,coef);
        }

        return logAC0Cache.get(an);
    }

    private double adaptiveSimpson(double[] likelihoods, double start, double stop, double err, int cap) {
        double mid = (start + stop)/2;
        double size = stop-start;
        double fa = calculateAFPosterior(likelihoods,start);
        double fb = calculateAFPosterior(likelihoods,mid);
        double fc = calculateAFPosterior(likelihoods,stop);
        double s = (size/6)*(fa + 4*fc + fb);
        double h = simpAux(likelihoods,start,stop,err,s,fa,fb,fc,cap);
        return h;
    }

    private double simpAux(double[] likelihoods, double a,double b,double eps,double s,double fa,double fb,double fc,double cap){
        if ( s == 0 )
                return -300.0;
        double c = ( a + b )/2;
        double h = b-a;
        double d = (a + c)/2;
        double e = (c + b)/2;
        double fd = calculateAFPosterior(likelihoods, d);
        double fe = calculateAFPosterior(likelihoods, e);
        double s_l = (h/12)*(fa + 4*fd + fc);
        double s_r = (h/12)*(fc + 4*fe + fb);
        double s_2 = s_l + s_r;
        if ( cap <= 0 || Math.abs(s_2 - s) <= 15*eps ){
            return Math.log10(s_2 + (s_2 - s)/15.0);
        }

        return MathUtils.approximateLog10SumLog10(simpAux(likelihoods,a,c,eps/2,s_l,fa,fc,fd,cap-1),simpAux(likelihoods, c, b, eps / 2, s_r, fc, fb, fe, cap - 1));
    }
}
