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

package org.broadinstitute.sting.multisamplecaller;

import org.broadinstitute.sting.utils.MathUtils;

import static java.lang.Math.log10;
import static java.lang.Math.pow;
import java.lang.Cloneable;

public class ClassicGenotypeLikelihoods implements Cloneable {
    // precalculate these for performance (pow/log10 is expensive!)
    private static final double[] oneMinusData = new double[Byte.MAX_VALUE];

    static {
        for (int qual = 0; qual < Byte.MAX_VALUE; qual++) {
            oneMinusData[qual] = log10(1.0 - pow(10, (qual / -10.0)));
            //oneMinusData[qual] = log10(1.0 - QualityUtils.qualToProb(qual));
        }
    }

    private static double getOneMinusQual(final byte qual) {
        return oneMinusData[qual];
    }

    private static final double[] oneHalfMinusData = new double[Byte.MAX_VALUE];

    static {
        for (int qual = 0; qual < Byte.MAX_VALUE; qual++) {
            oneHalfMinusData[qual] = log10(0.5 - pow(10, (qual / -10.0)) / 2.0);
            //oneHalfMinusData[qual] = log10(0.5 - QualityUtils.qualToProb(qual) / 2.0);
        }
    }

    private static double getOneHalfMinusQual(final byte qual) {
        return oneHalfMinusData[qual];
    }

    public double[] likelihoods;
    public String[] genotypes;
	public int coverage;

    // The genotype priors;
    private double priorHomRef;
    private double priorHet;
    private double priorHomVar;
    public String[] sorted_genotypes;
    public double[] sorted_likelihoods;
    double ref_likelihood = Double.NaN;
	private IndelLikelihood indel_likelihood;

    // Store the 2nd-best base priors for on-genotype primary bases
    //private HashMap<String, Double> onNextBestBasePriors = new HashMap<String, Double>();

    // Store the 2nd-best base priors for off-genotype primary bases
    //private HashMap<String, Double> offNextBestBasePriors = new HashMap<String, Double>();

    private static double[] p2ndon = {0.000, 0.302, 0.366, 0.142, 0.000, 0.548, 0.370, 0.000, 0.319, 0.000};
    private static double[] p2ndoff = {0.480, 0.769, 0.744, 0.538, 0.575, 0.727, 0.768, 0.589, 0.762, 0.505};

    public ClassicGenotypeLikelihoods() {
        initialize(1.0 - 1e-3, 1e-3, 1e-5, p2ndon, p2ndoff);
    }

    public ClassicGenotypeLikelihoods(boolean foo) {
    }

    public ClassicGenotypeLikelihoods(double priorHomRef, double priorHet, double priorHomVar) {
        initialize(priorHomRef, priorHet, priorHomVar, p2ndon, p2ndoff);
    }

    public ClassicGenotypeLikelihoods(double priorHomRef, double priorHet, double priorHomVar, double[] p2ndon, double[] p2ndoff) {
        initialize(priorHomRef, priorHet, priorHomVar, p2ndon, p2ndoff);
    }

    public ClassicGenotypeLikelihoods clone() {
        ClassicGenotypeLikelihoods c = new ClassicGenotypeLikelihoods(false);
        c.likelihoods = this.likelihoods.clone();
        c.genotypes = this.genotypes.clone();
	    c.coverage = this.coverage;

        // The genotype priors;
	    c.priorHomRef = this.priorHomRef;
	    c.priorHet = this.priorHet;
	    c.priorHomVar = this.priorHomVar;
        //public String[] sorted_genotypes;
        //public double[] sorted_likelihoods;
        //double ref_likelihood = Double.NaN;
	    //private IndelLikelihood indel_likelihood;
	    return c;
    }

    private void initialize(double priorHomRef, double priorHet, double priorHomVar, double[] p2ndon, double[] p2ndoff) {
        this.priorHomRef = priorHomRef;
        this.priorHet = priorHet;
        this.priorHomVar = priorHomVar;

        likelihoods = new double[10];
        genotypes = new String[10];
		coverage = 0;

		for (int i = 0; i < likelihoods.length; i++) { likelihoods[i] = Math.log10(0.1); }

        genotypes[0] = "AA";
        genotypes[1] = "AC";
        genotypes[2] = "AG";
        genotypes[3] = "AT";
        genotypes[4] = "CC";
        genotypes[5] = "CG";
        genotypes[6] = "CT";
        genotypes[7] = "GG";
        genotypes[8] = "GT";
        genotypes[9] = "TT";

        //for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
        //    onNextBestBasePriors.put(genotypes[genotypeIndex], p2ndon[genotypeIndex]);
        //    offNextBestBasePriors.put(genotypes[genotypeIndex], p2ndoff[genotypeIndex]);
        //}
    }

    public double getHomRefPrior() {
        return priorHomRef;
    }

    public void setHomRefPrior(double priorHomRef) {
        this.priorHomRef = priorHomRef;
    }

    public double getHetPrior() {
        return priorHet;
    }

    public void setHetPrior(double priorHet) {
        this.priorHet = priorHet;
    }

    public double getHomVarPrior() {
        return priorHomVar;
    }

    public void setHomVarPrior(double priorHomVar) {
        this.priorHomVar = priorHomVar;
    }

//    public double[] getOnGenotypeSecondaryPriors() {
//        double[] p2ndon = new double[10];
//
//        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
//            p2ndon[genotypeIndex] = onNextBestBasePriors.get(genotypes[genotypeIndex]);
//        }
//
//        return p2ndon;
//    }
//
//    public void setOnGenotypeSecondaryPriors(double[] p2ndon) {
//        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
//            onNextBestBasePriors.put(genotypes[genotypeIndex], p2ndon[genotypeIndex]);
//        }
//    }
//
//    public double[] getOffGenotypeSecondaryPriors() {
//        double[] p2ndoff = new double[10];
//
//        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
//            p2ndoff[genotypeIndex] = offNextBestBasePriors.get(genotypes[genotypeIndex]);
//        }
//
//        return p2ndoff;
//    }
//
//    public void setOffGenotypeSecondaryPriors(double[] p2ndoff) {
//        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
//            offNextBestBasePriors.put(genotypes[genotypeIndex], p2ndoff[genotypeIndex]);
//        }
//    }

    public void add(char ref, char read, byte qual) 
	{ 
		if (qual <= 0) { qual = 1; }

		if (coverage == 0)
		{
			for (int i = 0; i < likelihoods.length; i++)
			{
				likelihoods[i] = 0;
			}
		}
		double sum = 0.0;
        for (int i = 0; i < genotypes.length; i++) 
		{
			double likelihood = calculateAlleleLikelihood(ref, read, genotypes[i], qual); 
            likelihoods[i] += likelihood;
        }
		coverage += 1;
    }

    public void add(char ref, char read, byte qual, ConfusionMatrix matrix, String platform) 
	{ 
		if (qual <= 0) { qual = 1; }
		if (platform == null) { platform = "ILLUMINA"; }
		if (read == 'N') { return; }

		if (coverage == 0)
		{
			for (int i = 0; i < likelihoods.length; i++)
			{
				likelihoods[i] = 0;
			}
		}
		double sum = 0.0;
        for (int i = 0; i < genotypes.length; i++) 
		{
        	char h1 = genotypes[i].charAt(0);
        	char h2 = genotypes[i].charAt(1);

			double p1 = matrix.lookup(platform, read, h1);
			double p2 = matrix.lookup(platform, read, h2);

			double likelihood = calculateAlleleLikelihood(ref, read, genotypes[i], qual, p1, p2); 

			//System.out.printf("DBG: %c %c %s %d %f %f %f\n", ref, read, genotypes[i], qual, p1, p2, likelihood);

            likelihoods[i] += likelihood;
        }
		coverage += 1;
    }

    private double calculateAlleleLikelihood(char ref, char read, String genotype, byte qual) {
        if (qual == 0) {
            qual = 1;
        } // zero quals are wrong.

        char h1 = genotype.charAt(0);
        char h2 = genotype.charAt(1);

        double p_base;

        if ((h1 == h2) && (h1 == read)) {
            // hom
            p_base = getOneMinusQual(qual);
        } else if ((h1 != h2) && ((h1 == read) || (h2 == read))) {
            // het
            p_base = getOneHalfMinusQual(qual);
        } else {
            // error
            p_base = qual / -10.0;
        }

        return p_base;
    }

	public void TEST()
	{
		double p_A2A = 1.00;
		double p_T2T = 1.00;

		double p_A2T = 0.75;
		double p_T2A = 0.25;

		char ref = 'A';

		System.out.printf("\tA\tT\n");
		System.out.printf("A\t%.02f\t%.02f\n", p_A2A, p_A2T);
		System.out.printf("T\t%.02f\t%.02f\n", p_T2A, p_T2T);
		System.out.printf("\n");

		System.out.printf("P(A,Q20|AA) = %f\n", calculateAlleleLikelihood(ref, 'A', "AA", (byte)20, p_A2A, p_A2A));
		System.out.printf("P(A,Q20|AT) = %f\n", calculateAlleleLikelihood(ref, 'A', "AT", (byte)20, p_A2A, p_A2T));
		System.out.printf("P(A,Q20|TT) = %f\n", calculateAlleleLikelihood(ref, 'A', "TT", (byte)20, p_A2T, p_A2T));

		System.out.printf("P(T,Q20|AA) = %f\n", calculateAlleleLikelihood(ref, 'T', "AA", (byte)20, p_T2A, p_T2A));
		System.out.printf("P(T,Q20|AT) = %f\n", calculateAlleleLikelihood(ref, 'T', "AT", (byte)20, p_T2A, p_T2T));
		System.out.printf("P(T,Q20|TT) = %f\n", calculateAlleleLikelihood(ref, 'T', "TT", (byte)20, p_T2T, p_T2T));

		//System.exit(0);
	}

    private double calculateAlleleLikelihood(char ref, char read, String genotype, byte qual, double p1, double p2) {
        if (qual == 0) {
            qual = 1;
        } // zero quals are wrong.

        char h1 = genotype.charAt(0);
        char h2 = genotype.charAt(1);

		double perr = Math.pow(10.0,qual/-10.0);

        double p_base = 0;

		if (read == h1)
		{
			p_base += (1.0 - perr);
		}
		else
		{
			p_base += perr * p1;
		}

		if (read == h2)
		{
			p_base += (1.0 - perr);
		}
		else
		{
			p_base += perr * p2;
		}

		p_base = Math.log10(p_base/2.0);

        return p_base;
    }

    public void sort() {
        Integer[] permutation = MathUtils.sortPermutation(likelihoods);

        Integer[] reverse_permutation = new Integer[permutation.length];
        for (int i = 0; i < reverse_permutation.length; i++) {
            reverse_permutation[i] = permutation[(permutation.length - 1) - i];
        }

        sorted_genotypes = MathUtils.permuteArray(genotypes, reverse_permutation);
        sorted_likelihoods = MathUtils.permuteArray(likelihoods, reverse_permutation);
    }

    public String toString(char ref) {
        this.sort();
		double sum = 0;
        String s = String.format("%s %f %f ", this.BestGenotype(), this.LodVsNextBest(), this.LodVsRef(ref));
        for (int i = 0; i < sorted_genotypes.length; i++) {
            if (i != 0) {
                s = s + " ";
            }
            s = s + sorted_genotypes[i] + ":" + String.format("%.2f", sorted_likelihoods[i]);
			sum += Math.pow(10,sorted_likelihoods[i]);
        }
		s = s + String.format(" %f", sum);
        return s;
    }

    public void ApplyPrior(char ref, double[] allele_likelihoods)
	{
		int k = 0;
		for (int i = 0; i < 4; i++)
		{ 
			for (int j = i; j < 4; j++)
			{
				if (i == j) 
				{
					this.likelihoods[k] += Math.log10(allele_likelihoods[i]) + Math.log10(allele_likelihoods[j]);
				}
				else
				{
					this.likelihoods[k] += Math.log10(allele_likelihoods[i]) + Math.log10(allele_likelihoods[j]) + Math.log10(2);
				}
				k++;
			}
		}
		this.sort();
	}

    public void ApplyPrior(char ref, char alt, double p_alt) {
        for (int i = 0; i < genotypes.length; i++) {
            if ((p_alt == -1) || (p_alt <= 1e-6)) {
                if ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == ref)) {
                    // hom-ref
                    likelihoods[i] += Math.log10(priorHomRef);
                } else if ((genotypes[i].charAt(0) != ref) && (genotypes[i].charAt(1) != ref)) {
                    // hom-nonref
                    likelihoods[i] += Math.log10(priorHomVar);
                } else {
                    // het
                    likelihoods[i] += Math.log10(priorHet);
                }
                if (Double.isInfinite(likelihoods[i])) {
                    likelihoods[i] = -1000;
                }
            } else {
                if ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == ref)) {
                    // hom-ref
                    likelihoods[i] += 2.0 * Math.log10(1.0 - p_alt);
                } else if ((genotypes[i].charAt(0) == alt) && (genotypes[i].charAt(1) == alt)) {
                    // hom-nonref
                    likelihoods[i] += 2.0 * Math.log10(p_alt);
                } else if (((genotypes[i].charAt(0) == alt) && (genotypes[i].charAt(1) == ref)) ||
                        ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == alt))) {
                    // het
                    likelihoods[i] += Math.log10((1.0 - p_alt) * p_alt * 2.0);
                } else {
                    // something else (noise!)
                    likelihoods[i] += Math.log10(1e-5);
                }

                if (Double.isInfinite(likelihoods[i])) {
                    likelihoods[i] = -1000;
                }
            }
        }
        this.sort();
    }

	public void ApplyWeight(double weight)
	{
	    double log10weight = Math.log10(weight);
		for (int i = 0; i < genotypes.length; i++) { likelihoods[i] += log10weight; }
		this.sort();
	}

//    public void applySecondBaseDistributionPrior(String primaryBases, String secondaryBases) {
//        for (int genotypeIndex = 0; genotypeIndex < genotypes.length; genotypeIndex++) {
//            char firstAllele = genotypes[genotypeIndex].charAt(0);
//            char secondAllele = genotypes[genotypeIndex].charAt(1);
//
//            int offIsGenotypic = 0;
//            int offTotal = 0;
//
//            int onIsGenotypic = 0;
//            int onTotal = 0;
//
//            for (int pileupIndex = 0; pileupIndex < primaryBases.length(); pileupIndex++) {
//                char primaryBase = primaryBases.charAt(pileupIndex);
//
//                if (secondaryBases != null) {
//                    char secondaryBase = secondaryBases.charAt(pileupIndex);
//
//                    if (primaryBase != firstAllele && primaryBase != secondAllele) {
//                        if (secondaryBase == firstAllele || secondaryBase == secondAllele) {
//                            offIsGenotypic++;
//                        }
//                        offTotal++;
//                    } else {
//                        if (secondaryBase == firstAllele || secondaryBase == secondAllele) {
//                            onIsGenotypic++;
//                        }
//                        onTotal++;
//                    }
//                }
//            }
//
//            double offPrior = MathUtils.binomialProbability(offIsGenotypic, offTotal, offNextBestBasePriors.get(genotypes[genotypeIndex]));
//            double onPrior = MathUtils.binomialProbability(onIsGenotypic, onTotal, onNextBestBasePriors.get(genotypes[genotypeIndex]));
//
//            likelihoods[genotypeIndex] += Math.log10(offPrior) + Math.log10(onPrior);
//        }
//        this.sort();
//    }

    public double LodVsNextBest() {
        this.sort();
        return sorted_likelihoods[0] - sorted_likelihoods[1];
    }

    public double LodVsRef(char ref) {
        if ((this.BestGenotype().charAt(0) == ref) && (this.BestGenotype().charAt(1) == ref)) {
            ref_likelihood = sorted_likelihoods[0];
            return (-1.0 * this.LodVsNextBest());
        } else {
            for (int i = 0; i < genotypes.length; i++) {
                if ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == ref)) {
                    ref_likelihood = likelihoods[i];
                }
            }
        }
        return sorted_likelihoods[0] - ref_likelihood;
    }

    public String BestGenotype() {
        this.sort();
        return this.sorted_genotypes[0];
    }

    public double BestPosterior() {
        this.sort();
        return this.sorted_likelihoods[0];
    }

	public double RefPosterior(char ref)
	{
		this.LodVsRef(ref);
		return this.ref_likelihood;
	}

	public void addIndelLikelihood(IndelLikelihood indel_likelihood) { this.indel_likelihood = indel_likelihood; }
	public IndelLikelihood getIndelLikelihood() { return this.indel_likelihood; }

}
