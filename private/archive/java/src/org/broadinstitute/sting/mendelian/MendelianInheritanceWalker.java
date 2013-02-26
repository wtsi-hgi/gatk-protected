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


import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.Genotype;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.GenotypingCallStats;
import org.broadinstitute.sting.utils.TrioConcordanceRecord;
import org.broadinstitute.sting.utils.GenotypeUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.List;

//@Requires(value=DataSource.REFERENCE,referenceMetaData={@RMD(name="mother",type=rodSAMPileup.class),
//                                                        @RMD(name="father",type=rodSAMPileup.class),
//                                                        @RMD(name="daughter",type=rodSAMPileup.class)})
public class MendelianInheritanceWalker  extends RefWalker<TrioConcordanceRecord, TrioConcordanceRecord> {

	@Argument(fullName="point_consensus_cutoff", shortName="XPC", doc="confidence cutoff for consensus in point genotype", required=true ) public Double POINT_CONS_CUTOFF;
	@Argument(fullName="point_variant_cutoff", shortName="XPV", doc="confidence cutoff for variant (snp) in point genotype", required=true ) public Double POINT_VAR_CUTOFF;
	@Argument(fullName="indel_consensus_cutoff", shortName="XIC", doc="confidence cutoff for consensus in indel genotype", required=true ) public Double INDEL_CONS_CUTOFF;
	@Argument(fullName="indel_variant_cutoff", shortName="XIV", doc="confidence cutoff for variant (snp) in indel genotype", required=true ) public Double INDEL_VAR_CUTOFF;
	@Argument(fullName="log_concordant", shortName="LC",doc="If set, all trio-concordant sites will be logged at level INFO") public boolean LOG_CONCORDANT; 
	@Argument(fullName="log_discordant", shortName="LD",doc="If set, all trio-discordant sites will be logged at level INFO") public boolean LOG_DISCORDANT; 
	@Argument(fullName="default_reference_calls",shortName="DRC",
			doc="If set,  any position where the specified genotype is NOT explicitly specified, while the other is provided, is considered to be an implicit confident 'reference' (no-indel or no-snp) call") 
			public boolean defCalls;
	@Argument(fullName="variant_type",
					      shortName="VT",
					      doc="Assess concordance for the variants of the specified type, INDEL or POINT. If genotype track(s) provide both types, the requested one will be selected",
					      required=true) 
					      public String VTYPE_STR; 
	
	private static Logger logger = Logger.getLogger(MendelianInheritanceWalker.class);	
	private final static String star = new String("*");
	private GenotypeUtils.VariantType VARIANT_TYPE;
	
	@Override
	public TrioConcordanceRecord map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
				
//		String outLine = new String(context.getLocation() + " REF: "+ref + " RODS:" + rodData.getAllRods().size());
		
		ReferenceOrderedDatum rodMom = rodData.lookup("mother", null);
		ReferenceOrderedDatum rodDad = rodData.lookup("father", null);
		ReferenceOrderedDatum rodKid = rodData.lookup("daughter", null);
		
		Genotype mom = GenotypeUtils.extractGenotype(rodMom,VARIANT_TYPE,defCalls);
		Genotype dad = GenotypeUtils.extractGenotype(rodDad,VARIANT_TYPE,defCalls);
		Genotype kid = GenotypeUtils.extractGenotype(rodKid,VARIANT_TYPE,defCalls);

		return assessGenotypesInTrio(mom, dad, kid);
	}
	

/**	
 * 	@Override
 * @see org.broadinstitute.sting.gatk.walkers.Walker#initialize()
 */
	public void initialize() {
		super.initialize();
		VARIANT_TYPE = GenotypeUtils.VariantType.valueOf(VTYPE_STR.toUpperCase());
	};


	/** Takes a single genotype object and returns properly filled new assessment object (covered/assessed/ref/variant set to 0/1 
	 * according to what the genotype says)  
	 * @param g
	 * @return
	 */
	protected GenotypingCallStats assessGenotype(Genotype g, GenotypingCallStats stats) {
		
		if ( g != null ) stats.covered = 1;
		
		if ( hasCall(g)) {
			stats.assessed = 1;
			if ( g.isReference() ) stats.ref = 1;
			else {
				stats.variant = 1;
				if ( ! g.isBiallelic() ) stats.non_biallelic_variant = 1;
			}
		}
		return stats;
	}
	
	public TrioConcordanceRecord assessGenotypesInTrio(Genotype gMom, Genotype gDad, Genotype gKid) {		

		TrioConcordanceRecord t = new TrioConcordanceRecord();
		
//		String outLine = new String(context.getLocation() + " REF: "+ref + " RODS:" + rodData.getAllRods().size());

		// first get separate stats on each individual
		assessGenotype(gMom,t.mom);
		assessGenotype(gDad,t.dad);
		assessGenotype(gKid,t.kid);
		
		//		if ( hasCall(mom) && mom.isIndel() ) System.out.println("GOT INDEL: "+mom.toString());
		
		if ( t.mom.covered == 0 || t.dad.covered == 0 || t.kid.covered== 0 ) return t; // current position is not covered in at least one individual, there's nothing else to do
		t.trio.covered = 1; // ok, base covered in a trio (e.g. in all individuals)
		
		if ( t.mom.assessed != 0 &&  t.dad.assessed != 0 && t.kid.assessed != 0 ) {
			t.trio.assessed = 1; // assessed in trio = assessed in each individual
		} else return t; // at least one individual is not assessed, nothing left to do
			
		// we are here only if everyone is assessed
						
		// NOTE: "consistent_ref" in individuals is counted only when all 3 are assessed AND all 3 are ref (i.e. ref call in a trio)
		if( gMom.isReference() && gDad.isReference() && gKid.isReference() ) { // everyone is a ref
			t.trio.ref = t.trio.consistent_ref = 1;
			t.mom.consistent_ref = 1;
			t.dad.consistent_ref = 1;
			t.kid.consistent_ref = 1;
			return t; // done
		}

		// by now we know that there's a variant in at least one of the individuals

		t.trio.variant = 1;

		if ( t.mom.non_biallelic_variant == 1 || t.dad.non_biallelic_variant == 1 || t.kid.non_biallelic_variant == 1 ) {
			t.trio.non_biallelic_variant = 1;
			return t;
		}

		String kid_allele_1 = gKid.getFWDAlleles().get(0);
		String kid_allele_2 = gKid.getFWDAlleles().get(1);
		List<String> mom_alleles = gMom.getFWDAlleles();
		List<String> dad_alleles = gDad.getFWDAlleles();
		
		// warning: no special processing for X/Y chromosomes yet; not an issue for daughter

			
			
		if ( mom_alleles.contains(kid_allele_1) && dad_alleles.contains(kid_allele_2) ||
					mom_alleles.contains(kid_allele_2) && dad_alleles.contains(kid_allele_1) ) {
				t.trio.consistent_variant = 1;
				if ( ! gMom.isReference() ) { 
					t.mom.consistent_variant = 1 ;
					if ( ! gKid.isReference()  ) t.mom_passed_variant = 1;
				}
				if ( ! gDad.isReference() ) {
					t.dad.consistent_variant = 1 ;
					if ( ! gKid.isReference()  ) t.dad_passed_variant = 1;
				}
				if ( ! gKid.isReference() ) t.kid.consistent_variant = 1 ; 
				
				if ( LOG_CONCORDANT ) {
					logger.info("consistent variant at "+gMom.getLocation() + 
							"("+gMom.getFWDRefBases()+")  mom: " + genotypeString(gMom) + " dad: " +
							genotypeString(gDad) + " kid: " + genotypeString(gKid)
							);
				}
		} else {	
				// we are inconsistent. let's see what happened:
				
			   t.trio.inconsistent_variant = 1;

			   if ( ! gMom.isReference() ) t.mom.inconsistent_variant = 1 ;
			   if ( ! gDad.isReference() ) 	t.dad.inconsistent_variant = 1 ;
			   if ( ! gKid.isReference()  ) t.kid.inconsistent_variant = 1;

			   if ( gKid.isReference() && ( ! gMom.isReference() && gMom.isHom() || ! gDad.isReference() && gDad.isHom() ) ) t.missing_variant_in_kid = 1;
			   if ( ! gKid.isReference() && gMom.isReference() && gDad.isReference()  ) t.missing_variant_in_parents = 1;
			   if ( ! gKid.isReference() && ( ! gMom.isReference() || ! gDad.isReference() )  ) t.nonmatching_variant_in_kid = 1;
				
				if ( LOG_DISCORDANT ) {
					logger.info("inconsistent variant at "+gMom.getLocation() + 
							"("+gMom.getFWDRefBases()+")  mom: " + genotypeString(gMom) + " dad: " +
							genotypeString(gDad) + " kid: " + genotypeString(gKid)
							);
				}

		}
				
		return t;
		
	}
	
	
	protected String alleleString(Genotype g, int n) {
		if ( g.getFWDAlleles().get(n).length() == 0 ) return star;
		return g.getFWDAlleles().get(n);
	}
	
	protected String genotypeString(Genotype g) {
		return alleleString(g, 0) +"/"+alleleString(g, 1);
	}

	@Override
	public TrioConcordanceRecord reduce(TrioConcordanceRecord value, TrioConcordanceRecord sum) {
		return sum.add(value);
	}

	@Override
	public TrioConcordanceRecord reduceInit() {
		
		return new TrioConcordanceRecord();
	}

	boolean hasCall(Genotype g) {
		if ( g == null ) return false; // there's no call if there's no rod data available, duh!

		if ( g.isReference() ) {
			if ( g.isPointGenotype() ) return g.getConsensusConfidence() >= POINT_CONS_CUTOFF ;
			else return g.getConsensusConfidence() >= INDEL_CONS_CUTOFF ;
		}
		else { // it's a variant
			if ( g.isPointGenotype() ) return g.getVariantConfidence() >= POINT_VAR_CUTOFF ;
			else return g.getVariantConfidence() >= INDEL_VAR_CUTOFF ;
		}
		
	}
	
	
	/*
	protected String shortLine(Genotype av) {
		
		if ( av == null )  return new String( "<NULL>");
		
		if ( av.isReference() ) return new String ("*");
		
		List<String> alleles = av.getFWDAlleles();
		
		if ( av.isSNP() ) {
            if  ( alleles.get(0).charAt(0) == av.getRef() ) return alleles.get(1);
            else return alleles.get(0);
		}
		if ( av.isInsertion() ) {
            if ( alleles.get(0).length() == 0 ) return new String('+'+alleles.get(1));
            else return new String('+'+alleles.get(0));
		}
		if ( av.isDeletion() ) {
            if ( alleles.get(0).length() == 0 ) return new String ('-'+alleles.get(0));
            else return new String('-'+alleles.get(1));
		}
		if ( av.isIndel() ) {
			return new String('?'+alleles.get(0));
		}
		return new String("<missing data!!!>");
	}
	*/
}
