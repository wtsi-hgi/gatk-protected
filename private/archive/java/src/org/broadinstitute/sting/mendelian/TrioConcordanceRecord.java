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

package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.utils.Utils;

/**
 * This class is a trivial wrapper for keeping together and passing around counts of different possible outcomes of 
 * comparisons in a trio.  Mendelian walker uses this class to classify/accumulate events such as consistent snp, inconsistent snp 
 * (e.g. only in a kid, but not in parents), loci with no calls etc. 
 * @author asivache
 *
 */
public class TrioConcordanceRecord {
	
	public GenotypingCallStats mom;
	public GenotypingCallStats dad;
	public GenotypingCallStats kid;
	
	public GenotypingCallStats trio;

//	public long mom_assessed_ref; // number of ref calls in mother on positions assessed *in all 3 individuals*
//	public long dad_assessed_ref; // ditto
//	public long kid_assessed_ref; 
//	public int mom_assessed_variant; // number of variant calls in mother on  positions assessed *in all 3 individuals*
//	public int dad_assessed_variant; // ditto
//	public int kid_assessed_variant; 
	public int missing_variant_in_kid;
	public int nonmatching_variant_in_kid;
	public int missing_variant_in_parents;
	public int mom_passed_variant;
	public int dad_passed_variant;

	//	public long consistent_ref = 0; // number of assessed loci, where all 3 people have homogeneous reference allele
//	public int consistent_variant = 0; // number of assessed loci where a variant is observed in at least one individual and genotyping calls are consistent between the trio members
//	public int inconsistent_variant = 0; // number of assessed loci where a variant is observed in at least one individual and genotyping calls are inconsistent
//	public int missing_variant_in_parents = 0; // number of inconsistent variants (see above), where parent(s) have a variant but the kid does not while she should
//	public int missing_variant_in_kid = 0; // number of inconsistent variants (see above), where kid has a snp but the parents do not while they should
//	public int consistent_variant_passed = 0; // variants that are consistent and *passed* (i.e. present in kid and one of the parents)
//	public int non_biallelic_variant = 0; // number of variant calls that are not biallelic
//	public long unclassified_events = 0;
	
	public TrioConcordanceRecord() {
		mom = new GenotypingCallStats();
		dad = new GenotypingCallStats();
		kid = new GenotypingCallStats();
		trio = new GenotypingCallStats();
	}
	
	public TrioConcordanceRecord add(TrioConcordanceRecord other) {
		
		this.mom.add(other.mom);
		this.dad.add(other.dad);
		this.kid.add(other.kid);

		this.trio.add(other.trio);

//		this.mom_assessed_ref += other.mom_assessed_ref;
//		this.dad_assessed_ref += other.dad_assessed_ref;
//		this.kid_assessed_ref += other.kid_assessed_ref;
//		this.mom_assessed_variant += other.mom_assessed_variant;
//		this.dad_assessed_variant += other.dad_assessed_ref;
//		this.kid_assessed_variant += other.kid_assessed_variant;
		this.missing_variant_in_kid += other.missing_variant_in_kid ;
		this.nonmatching_variant_in_kid += other.nonmatching_variant_in_kid ;
		this.missing_variant_in_parents += other.missing_variant_in_parents ;
		this.mom_passed_variant += other.mom_passed_variant;
		this.dad_passed_variant += other.dad_passed_variant;

		//		this.consistent_ref += other.consistent_ref;
//		this.consistent_variant += other.consistent_variant;
//		this.inconsistent_variant += other.inconsistent_variant;
//		this.missing_variant_in_parents += other.missing_variant_in_parents;
//		this.missing_variant_in_kid += other.missing_variant_in_kid;
//		this.consistent_variant_passed += other.consistent_variant_passed;
//		this.non_biallelic_variant += other.non_biallelic_variant;
//		this.unclassified_events += other.unclassified_events;
		return this;
	}
	
	public int totalVariants() { return trio.consistent_variant + trio.inconsistent_variant + trio.non_biallelic_variant; }
	
	public String toString() {
		StringBuilder b = new StringBuilder();
		
		b.append(String.format("%ncovered in trio: %d%n", trio.covered ) );

		b.append(String.format("assessed in trio: %d (%3.2f%% covered)%n", 
				trio.assessed, Utils.percentage(trio.assessed,trio.covered )) );
		
		b.append(String.format("   reference in all samples: %d (%3.2f%% assessed)%n", 
				trio.ref, Utils.percentage(trio.ref,trio.assessed )) );
		
		b.append(String.format("   variant sites: %d (%3.2f%% assessed, or 1 per %3.2f kB)%n", 
				totalVariants(), Utils.percentage(totalVariants(), trio.assessed), ((double)trio.assessed/totalVariants())/1000.0 
		));
		
		b.append(String.format("      consistent variants: %d (%3.2f%% variants)%n", 
				trio.consistent_variant, Utils.percentage(trio.consistent_variant,totalVariants()) 
		));

//		b.append(String.format("         passed (in daughter and parent(s)): %d%n         lost (in parent(s) but not in daughter): %d%n",
//    		   consistent_variant_passed, consistent_variant - consistent_variant_passed));
       
       b.append(String.format("      multiallelic variant: %d (%3.2f%% variants)%n",
				trio.non_biallelic_variant, Utils.percentage(trio.non_biallelic_variant, totalVariants())
       ));

       b.append(String.format("      inconsistent variant: %d (%3.2f%% variants)%n",
				trio.inconsistent_variant, Utils.percentage(trio.inconsistent_variant, totalVariants())
      ));

       b.append(String.format("         missing from daughter: %d (%3.2f%% inconsistent variants)%n",
				missing_variant_in_kid, Utils.percentage(missing_variant_in_kid, trio.inconsistent_variant)
      ));

       b.append(String.format("         missing from both parents: %d (%3.2f%% inconsistent variants)%n",
				missing_variant_in_parents, Utils.percentage(missing_variant_in_parents, trio.inconsistent_variant)		
       ));

       b.append(String.format("         non-matching in daughter: %d (%3.2f%% inconsistent variants)%n",
				nonmatching_variant_in_kid, Utils.percentage(nonmatching_variant_in_kid, trio.inconsistent_variant)		
       ));
       
		b.append("per trio individual:\n");
		b.append("   mother:\n");
		b.append(mom.toString());
		b.append("   father:\n");
		b.append(dad.toString());
		b.append("   daughter:\n");
		b.append(kid.toString());
		
		return b.toString();
	}
}
