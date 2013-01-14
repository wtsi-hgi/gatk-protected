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
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.VCFReader;


import java.io.*;
import java.util.*;

//import org.apache.commons.math.optimization.*;
//import org.apache.commons.math.optimization.direct.*;
//import org.apache.commons.math.analysis.MultivariateRealFunction;

// Program for frequency-specific  VCF-files.


/**
 * @author jmaguire
 */


class VCFOptimize extends CommandLineProgram 
{
		@Argument(fullName = "vcf", shortName = "vcf", doc = "file to open", required = true) public String filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "target_TsTv", shortName = "target_TsTv", doc = "Minimum acceptable TsTv", required=false) public double target_TsTv = 2.07;
		@Argument(fullName = "output", shortName = "output", doc = "file to write cuts to", required = true) public String output_filename;
		@Argument(fullName = "min_calls", shortName = "min_calls", doc = "Minimum signifigant number of calls", required=false) public int min_calls = 100;
		@Argument(fullName = "num_breaks", shortName = "num_breaks", doc = "Number of breaks to search over", required=false) public int num_breaks = 100;

		// Debugging arguments:
		@Argument(fullName = "n_records", shortName = "n_records", doc = "Number of records to load (debugging)", required=false) public int n_records_to_process = Integer.MAX_VALUE;
		@Argument(fullName = "verbose", shortName = "verbose", doc = "print detailed debugging info.", required=false) public boolean verbose = false;

		class OptimizationRecord
		{
			public boolean transition;
			public int freq;
			public double[] features;
			
			public OptimizationRecord(boolean transition, int freq, double[] features)
			{
				this.transition = transition;
				this.freq       = freq;
				this.features   = features.clone();
			}

			public OptimizationRecord clone()
			{
				return new OptimizationRecord(transition, freq, features.clone());
			}
		}

		private OptimizationRecord pack(VCFHeader header, VCFRecord input)
		{
			Map<String,String> info = input.getInfoValues();
			
			if (! info.containsKey("AC")) 
			{ 
				throw new RuntimeException("AC not present in record: \n" + input.toStringEncoding(header));
			}
			if (! info.containsKey("DP")) 
			{ 
				throw new RuntimeException("DP not present in record: \n" + input.toStringEncoding(header));
			}
			if (! info.containsKey("SB")) 
			{ 
				throw new RuntimeException("SB not present in record: \n" + input.toStringEncoding(header));
			}

			boolean transition = VCFTool.isTransition(input);
			int freq           = Integer.parseInt(input.getInfoValues().get("AC"));
			double LOD         = input.getQual();
			double depth       = Double.parseDouble(input.getInfoValues().get("DP"));
			double SLOD        = Double.parseDouble(input.getInfoValues().get("SB"));

			double[] features = new double[2];
			features[0] = LOD;
			features[1] = -1*SLOD;
	
			return new OptimizationRecord(transition, freq, features);
		}

			// This is the objective function we're searching in.
			// if (tstv>=min) { return #snps; } else { return -inf; }
			public double tstv(double[] point, OptimizationRecord[] records)
			{
				double transitions   = 0;
				double transversions = 0;
				double total         = 0;
				for (int i = 0; i < records.length; i++)
				{
					int j = 0;
					for (j = 0; j < point.length; j++)
					{
//						if (records == null) { System.out.printf("records==null\n"); }
//						if (records[i] == null) { System.out.printf("records[%d]==null\n", i); }
//						if (records[i].features == null) { System.out.printf("records[%d].features==null\n", i); }

						if (records[i].features[j] < point[j]) { break; }
					}
					if (j == point.length)
					{
						if (records[i].transition == true) { transitions += 1; }
						else { transversions += 1; }
						total += 1;
					}
				}

				double tstv = transitions / transversions;
				return tstv;
			}


			// This is the objective function we're searching in.
			// if (tstv>=min) { return #snps; } else { return -inf; }
			public double num_calls(double[] point, OptimizationRecord[] records)
			{
				double total = 0;
				for (int i = 0; i < records.length; i++)
				{
					int j = 0;
					for (j = 0; j < point.length; j++)
					{
//						if (records == null) { System.out.printf("records==null\n"); }
//						if (records[i] == null) { System.out.printf("records[%d]==null\n", i); }
//						if (records[i].features == null) { System.out.printf("records[%d].features==null\n", i); }

						if (records[i].features[j] < point[j]) { break; }
					}
					if (j == point.length)
					{
						total += 1;
					}
				}

				return total;
			}


			public class Cut
			{
				public double lod;
				public double slod;
				public int freq;

				public Cut(double lod, double slod)
				{
					this.lod = lod;
					this.slod = slod;
					this.freq = -1;
				}

				public Cut(double lod, double slod, int freq)
				{
					this.lod = lod;
					this.slod = slod;
					this.freq = freq;
				}

				public Cut(String record)
				{
					String[] tokens = record.split("\\s+");
					this.freq = Integer.parseInt(tokens[0]);					
					this.lod  = Double.parseDouble(tokens[1]);					
					this.slod = Double.parseDouble(tokens[2]);					
				}

				public String toString()
				{
					return String.format("%d %f %f", freq, lod, slod);
				}
			}


		// Just a simple grid search.
		private Cut optimize(OptimizationRecord[] records, double min_TsTv, int freq)
		{


			double[] lods  = new double[records.length];
			double[] slods = new double[records.length];
			for (int i = 0; i < lods.length; i++)
			{
				lods[i]  = records[i].features[0];
				slods[i] = records[i].features[1];
			}

			Arrays.sort(lods);
			Arrays.sort(slods);

			double[] lod_breaks  = new double[num_breaks];
			double[] slod_breaks = new double[num_breaks];
			int bin_size = 1 + (records.length / num_breaks);

			//System.out.printf("BREAKS i j lod slod\n");
			int j = 0;
			for (int i = 0; i < records.length; i += bin_size)
			{
				lod_breaks[j] = lods[i];
				slod_breaks[j] = slods[i];
				j += 1;
				//System.out.printf("BREAKS %d %d %f %f\n", i, j, lods[i], slods[i]);
			}
			//System.out.printf("\n");

			double best_lod       = lod_breaks[0];
			double best_slod      = slod_breaks[0];

			int best_lod_idx       = 0;
			int best_slod_idx      = 0;

			double[] point = new double[2];
			point[0] = best_lod;
			point[1] = best_slod;

			double best_tstv      = tstv(point, records);
			double best_num_calls = num_calls(point, records);
			boolean flag = false;

			//for (double lod = 0; lod < 8000; lod += 10)
			for (int lod_idx = 0; lod_idx < num_breaks; lod_idx += 1)
			{
				double lod = lod_breaks[lod_idx];
				//for (double slod = -4000; slod < 1000; slod += 10)
				for (int slod_idx = 0; slod_idx < num_breaks; slod_idx += 1)
				{
					double slod = slod_breaks[slod_idx];

					point = new double[2];
					point[0] = lod;
					point[1] = slod;
					double tstv      = tstv(point, records);
					double num_calls = num_calls(point, records);
					
					if (num_calls < min_calls) { continue; }

					if ((tstv >= min_TsTv) && (num_calls > best_num_calls)) 
					{ 
						best_lod=lod; 
						best_slod=slod; 
						best_tstv=tstv; 
						best_num_calls=num_calls; 
						best_lod_idx = lod_idx;
						best_slod_idx = slod_idx;
						flag=true;
					}
					else if ((tstv >= best_tstv) && (!flag)) 
					{ 
						best_lod=lod; 
						best_slod=slod; 
						best_tstv=tstv; 
						best_num_calls=num_calls;
						best_lod_idx = lod_idx;
						best_slod_idx = slod_idx;
					}


					if (verbose)
					{
						System.out.printf("DEBUG: %d | %d %d | %f %f %f %f | %f %f %f %f\n", 
									freq, 
									lod_idx, slod_idx,
									lod, slod, num_calls, tstv,
									best_lod, best_slod, best_num_calls, best_tstv);
					}
				}
			}
			
			//System.out.printf("Found optimum: lod=%f slod=%f num_calls=%f tstv=%f\n", best_lod, best_slod, best_num_calls, best_tstv);
			System.out.printf("%d %d %d %f %f %f %f\n", freq, best_lod_idx, best_slod_idx, best_lod, best_slod, best_num_calls, best_tstv);

			return new Cut(best_lod, best_slod);
		}

		@Override
		protected int execute() 
		{
			System.out.println("Loading " + filename + "...");
		
			VCFReader reader = null;

			if (autocorrect) { reader = new VCFReader(new File(filename),new VCFHomogenizer()); }
			else { reader = new VCFReader(new File(filename)); }

			PrintWriter output = null;
			try
			{
				output = new PrintWriter(new FileWriter(output_filename));
			}
			catch (Exception e)
			{
				throw new RuntimeException(e); 
			}

			VCFHeader header = reader.getHeader();

			HashMap<Integer,ArrayList<OptimizationRecord>> records = new HashMap<Integer,ArrayList<OptimizationRecord>>();

			Date start_time = new Date();
			int n_records_processed = 0;
			int max_freq = 0;
			while(reader.hasNext())
			{
				VCFRecord record = reader.next();
		
				OptimizationRecord optimization_record = pack(header, record);

				if (optimization_record.freq > max_freq) { max_freq = optimization_record.freq; }

				if (! records.containsKey(optimization_record.freq)) { records.put(optimization_record.freq, new ArrayList<OptimizationRecord>()); }
				records.get(optimization_record.freq).add(optimization_record.clone());

				n_records_processed += 1;

				if (n_records_processed == n_records_to_process) { break; }
			}
			System.out.printf("Loaded %d records\n", n_records_processed);

			//for (int freq = 1; freq <= 5; freq += 1)
			for (int freq = 1; freq <= max_freq; freq += 1)
			{
				if (records.get(freq) == null) { System.out.printf("Skipping AAF %d (no calls)\n", freq); continue; }
				System.out.printf("\nOptimizing AAF %d...\n", freq);

				OptimizationRecord[] fnord = new OptimizationRecord[records.get(freq).size()];
				Cut cut = optimize(records.get(freq).toArray(fnord), target_TsTv, freq);
				cut.freq = freq;

				output.println(cut);
			}
			output.flush();
			output.close();
			
			return 0;
		}
}
