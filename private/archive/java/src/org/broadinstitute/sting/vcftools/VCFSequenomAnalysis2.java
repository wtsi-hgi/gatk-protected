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

package org.broadinstitute.sting.tools.vcf;
import org.broad.tribble.vcf.VCFGenotypeEncoding;
import org.broad.tribble.vcf.VCFGenotypeRecord;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Argument;

import java.io.*;
import java.util.*;
import java.lang.*;

import net.sf.picard.util.Interval;
import org.broadinstitute.sting.utils.genotype.vcf.VCFReader;


class VCFSequenomAnalysis2 extends CommandLineProgram 
{
		@Argument(fullName = "sequenom", shortName = "sequenom", doc = "file to open", required = true) public String filename1;
		@Argument(fullName = "sequencing", shortName = "sequencing", doc = "file to open", required = true) public String filename2;
		@Argument(fullName = "out", shortName = "out", doc = "file to write results to", required = false) public String output_filename = "/dev/stdout";
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = true;
		@Argument(fullName = "verbose", shortName = "verbose", doc = "print extremely detailed stats", required = false) public Boolean verbose = false;
		@Argument(fullName = "qual_threshold", shortName = "qual_threshold", doc = "minimum genotype quality to consider", required = false) public long qual_threshold = 1;


		@Override
		protected int execute() 
		{
			//System.out.println("Loading " + filename + "...");
		
			PrintStream output = null;
			try
			{
				output = new PrintStream(new FileOutputStream(output_filename));
			}
			catch (Exception e)
			{
				throw new RuntimeException(e);
			}

			output.printf("PROBE interval flag ref alt n_total_sequenom failure_rate_sequenom n_alt_sequencing n_alt_sequenom p_alt_sequencing p_alt_sequenom HWE_sequencing_chi HWE_sequenom_chi HWE_sequencing_p HWE_sequenom_p is_singleton_in_sequencing singleton_matched_in_sequenom n_sequencing_hets n_sequenom_hets num_hets_also_het_in_sequenom num_hets_dropped_in_sequenom\n");

			VCFReader reader1;
			VCFReader reader2;

			if (autocorrect) 
			{ 
				reader1 = new VCFReader(new File(filename1),new VCFHomogenizer());
				reader2 = new VCFReader(new File(filename2),new VCFHomogenizer());  
			}
			else 
			{ 
				reader1 = new VCFReader(new File(filename1)); 
				reader2 = new VCFReader(new File(filename2)); 
			}
		
			VCFHeader header1 = reader1.getHeader();
			VCFHeader header2 = reader2.getHeader();

			VCFRecord record1 = reader1.next();
			VCFRecord record2 = reader2.next();

			int[] sequenom_aaf_counts = new int[1024];
			int[] sequencing_aaf_counts = new int[1024];
			int max_aaf = 0;

			while(true)
			{
				if (record1 == null) { break; }
				if (record2 == null) { break; }

				Interval interval1 = VCFTool.getIntervalFromRecord(record1);
				Interval interval2 = VCFTool.getIntervalFromRecord(record2);

				int comparison = interval1.compareTo(interval2);
				
				if (comparison == 0)
				{
					// records match! compute concordance.
					
					// (unless one of them is "filtered")
					if (record1.isFiltered() || record2.isFiltered())
					{
						record1 = reader1.next();
						record2 = reader2.next();
						continue;
					}

					String[] sample_names = record2.getSampleNames();
					
					char ref = record1.getReference().charAt(0);
					char alt = VCFTool.getAlt(record2);
					
					int n_total_sequenom = VCFTool.Compute_n_total(record1, sample_names);
					int n_total_sequencing = VCFTool.Compute_n_total(record2, sample_names);
					double failure_rate_sequenom = VCFTool.Compute_failure_rate(record1);

					int n_alt_sequenom = VCFTool.Compute_n_alt(record1, sample_names);
					int n_alt_sequencing = VCFTool.Compute_n_alt(record2, sample_names);

					double p_alt_sequenom = (double)n_alt_sequenom / (double)n_total_sequenom;
					double p_alt_sequencing = (double)n_alt_sequencing / (double)n_total_sequencing;

					int n_het_sequenom = VCFTool.Compute_n_het(record1, sample_names);
					int n_het_sequencing = VCFTool.Compute_n_het(record2, sample_names);

					sequenom_aaf_counts[n_alt_sequenom] += 1;
					sequencing_aaf_counts[n_alt_sequencing] += 1;
					if (n_alt_sequenom > max_aaf) { max_aaf = n_alt_sequenom; }
					if (n_alt_sequencing > max_aaf) { max_aaf = n_alt_sequencing; }

					double HWE_sequenom = VCFTool.Compute_HWE(record1, sample_names);
					double HWE_sequencing = VCFTool.Compute_HWE(record2, sample_names);

					boolean isPolymorphic_sequenom   = (n_alt_sequenom   > 0) ? true : false;
					boolean isPolymorphic_sequencing = (n_alt_sequencing > 0) ? true : false;

					int is_singleton_in_sequencing = 0;
					int singleton_matched_in_sequenom = 0;
					
					if ((n_alt_sequencing == 1) && (VCFTool.Compute_n_alt(record1) > 0))
					{
						is_singleton_in_sequencing = 1;
						singleton_matched_in_sequenom = CheckSingletonMatch(record2, record1);
					}

					int[] het_match_ans = ComputeHetMatches(record2, record1);
					int num_hets_also_het_in_sequenom = het_match_ans[0];
				   	int num_hets_dropped_in_sequenom = het_match_ans[1];

					String flag = null;
					if (isPolymorphic_sequenom) { flag = "TP"; }
					else { flag = "FP"; }

					output.printf("PROBE %s %s %c %c %d %f %d %d %f %f %f %f %f %f %d %d %d %d %d %d\n",
										interval1,
										flag,
										ref,
										alt,
										n_total_sequenom,
										failure_rate_sequenom,
										n_alt_sequencing,
										n_alt_sequenom,
										p_alt_sequencing,
										p_alt_sequenom,
										HWE_sequencing,
										HWE_sequenom,
										VCFTool.P_from_Chi(HWE_sequencing),
										VCFTool.P_from_Chi(HWE_sequenom),
										is_singleton_in_sequencing,
										singleton_matched_in_sequenom,
										n_het_sequencing,
										n_het_sequenom,
										num_hets_also_het_in_sequenom,
				   						num_hets_dropped_in_sequenom);

					record1 = reader1.next();
					record2 = reader2.next();
				}
				else if (comparison > 0)
				{
					// interval1 is later than interval2.
					//System.err.printf("Skipping (2): %s\n", VCFTool.getIntervalFromRecord(record2));	
					record2 = reader2.next();
				}
				else if (comparison < 0)
				{
					// interval2 is later than interval1.
					//System.err.printf("Skipping (1): %s\n", VCFTool.getIntervalFromRecord(record1));	
					record1 = reader1.next();
				}

			}

			for (int i = 0; i < max_aaf; i++)
			{
				output.printf("AAF %d %d %d\n", i, sequenom_aaf_counts[i], sequencing_aaf_counts[i]);	
			}


			output.flush();
			output.close();


			return 0;
		}

		int CheckSingletonMatch(VCFRecord sequencing, VCFRecord sequenom)
		{
			String singleton_name = "";

			// first, check sequencing
			String[] sample_names = sequencing.getSampleNames();
			List<VCFGenotypeRecord> genotypes = sequencing.getVCFGenotypeRecords();
			int n_ref = 0;
			int n_alt = 0;
			for (int i = 0; i < sample_names.length; i++)
			{
				VCFGenotypeRecord rec = genotypes.get(i);
				List<VCFGenotypeEncoding> alleles = rec.getAlleles();
				String g = "";
				for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
				char[] c = g.toCharArray();
				Arrays.sort(c);
				g = new String(c);
				if (g.equals("..")) { continue; }
				if (g.charAt(0) == sequencing.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; singleton_name = sample_names[i]; }
				if (g.charAt(1) == sequencing.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; singleton_name = sample_names[i]; }
			}
			if (n_alt != 1) { throw new RuntimeException(); }
			if (singleton_name.equals("")) { throw new RuntimeException(); }

			// now, check sequenom
			sample_names = sequenom.getSampleNames();
			genotypes = sequenom.getVCFGenotypeRecords();
			n_ref = 0;
			n_alt = 0;
			for (int i = 0; i < sample_names.length; i++)
			{
				if (sample_names[i].equals(singleton_name))
				{
					VCFGenotypeRecord rec = genotypes.get(i);
					List<VCFGenotypeEncoding> alleles = rec.getAlleles();
					String g = "";
					for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
					char[] c = g.toCharArray();
					Arrays.sort(c);
					g = new String(c);
					if (g.equals("..")) { continue; }
					if (g.charAt(0) == sequenom.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; singleton_name = sample_names[i]; }
					if (g.charAt(1) == sequenom.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; singleton_name = sample_names[i]; }
					break;
				}
			}
			if (n_alt > 0) { return 1; }
			else if (n_ref != 0) { return 0; }
			else { return -1; }
		}

		int[] ComputeHetMatches(VCFRecord sequencing, VCFRecord sequenom)
		{
			// first, check sequencing
			String[] sample_names = sequencing.getSampleNames();
			List<VCFGenotypeRecord> genotypes = sequencing.getVCFGenotypeRecords();
			ArrayList<String> het_samples = new ArrayList<String>();
			for (int i = 0; i < sample_names.length; i++)
			{
				int n_ref = 0;
				int n_alt = 0;

				VCFGenotypeRecord rec = genotypes.get(i);
				List<VCFGenotypeEncoding> alleles = rec.getAlleles();
				String g = "";
				for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
				char[] c = g.toCharArray();
				Arrays.sort(c);
				g = new String(c);
				if (g.equals("..")) { continue; }
				if (g.charAt(0) == sequencing.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; }
				if (g.charAt(1) == sequencing.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; }
				if (n_alt == 1) { het_samples.add(sample_names[i]); }
			}

			// now, check sequenom
			sample_names = sequenom.getSampleNames();
			genotypes = sequenom.getVCFGenotypeRecords();
			int matched_hets = 0;
			int dropped_hets = 0;
			int mismatched_hets = 0;
			int num_hets = het_samples.size();
			for (int i = 0; i < sample_names.length; i++)
			{
				if (het_samples.contains(sample_names[i]))
				{
					het_samples.remove(sample_names[i]);

					int n_ref = 0;
					int n_alt = 0;
					VCFGenotypeRecord rec = genotypes.get(i);
					List<VCFGenotypeEncoding> alleles = rec.getAlleles();
					String g = "";
					for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
					char[] c = g.toCharArray();
					Arrays.sort(c);
					g = new String(c);
					if (g.equals("..")) { dropped_hets += 1; continue; }
					if (g.charAt(0) == sequenom.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; }
					if (g.charAt(1) == sequenom.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; }
					if (n_alt == 1) { matched_hets += 1; }	
					else { mismatched_hets += 1; }
				}
			}

			if ((matched_hets + dropped_hets + mismatched_hets) != num_hets)
			{
				String warning = String.format("WARNING: %d + %d + %d != %d ", 
															matched_hets,
															dropped_hets,
															mismatched_hets,
															num_hets);
				for (int i = 0; i < het_samples.size(); i++) { warning += het_samples.get(i) + " "; }
				System.out.println(warning);
			}

			int[] ans = new int[2];
			ans[0] = matched_hets;
			ans[1] = dropped_hets;
			return ans;
		}
}
