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

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.indels.Matrix;

import java.util.*;
import java.io.*;

// Beta iterative multi-sample caller
// j.maguire 6-11-2009

public class MultiSampleCallerAccuracyTest extends MultiSampleCaller
{
    @Argument(required=false, shortName="lod_threshold", doc="") public double LOD_THRESHOLD = 1e-6;
    @Argument(required=true, shortName="stats_output", doc="") public String STATS_OUTPUT;

	Matrix<Integer> n_variants;
	Matrix<Integer> n_found;

	PrintStream stats_output;

    public void initialize() 
	{
		this.DISCOVERY_OUTPUT = "/dev/null";
		this.INDIVIDUAL_OUTPUT = "/dev/null";

		super.initialize();

		n_variants = new Matrix<Integer>(sample_names.size()*2, sample_names.size()*2);
		n_found    = new Matrix<Integer>(sample_names.size()*2, sample_names.size()*2);

		for (int i = 0; i < sample_names.size()*2; i++)
		{
			for (int j = 0; j < sample_names.size()*2; j++)
			{
				n_variants.set(i,j,0);
				n_found.set(i,j,0);
			}
		}

		try
		{
			stats_output = new PrintStream(STATS_OUTPUT);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}

	}

    public MultiSampleCallResult map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) 
	{
        HapMapGenotypeROD hapmap = (HapMapGenotypeROD)tracker.lookup("hapmap", null);

		// Collect all the variants and the normals.
		ArrayList<String> variant_samples = new ArrayList<String>();
		ArrayList<String> reference_samples = new ArrayList<String>();

		int n_ref_chromosomes = 0;
		int n_alt_chromosomes = 0;

		String reference_genotype = String.format("%c%c", ref, ref);
		for (int i = 0; i < sample_names.size(); i++)
		{
			String true_genotype = hapmap.get(sample_names.get(i));
			if (true_genotype == null) { continue; }

			if (true_genotype.equals(reference_genotype)) { reference_samples.add(sample_names.get(i)); }
			else { variant_samples.add(sample_names.get(i)); }

			if (true_genotype.equals(reference_genotype)) { n_ref_chromosomes += 1; }
			else if (true_genotype.contains(String.format("%c",ref))) { n_ref_chromosomes += 1; n_alt_chromosomes += 1; }
			else { n_alt_chromosomes += 2; }

		}

			// Put together a context.
			ArrayList<String> working_samples = new ArrayList<String>();
			working_samples.addAll(variant_samples);
			working_samples.addAll(reference_samples);
			AlignmentContext working_context = filterAlignmentContextBySamples(context, working_samples);

			// Call.
			MultiSampleCallResult call_result = super.map(tracker, ref, working_context);
			EM_Result em_result = call_result.em_result;

			// Compute Statistics.
			if (n_variants == null) { System.out.printf("n_variants is null\n"); }
			if (n_found == null) { System.out.printf("n_found is null\n"); }
			n_variants.set(n_ref_chromosomes, n_alt_chromosomes, n_variants.get(n_ref_chromosomes, n_alt_chromosomes)+1);
			if ((call_result.lod > LOD_THRESHOLD) && (n_alt_chromosomes >= 1))
			{
				n_found.set(n_ref_chromosomes, n_alt_chromosomes, n_found.get(n_ref_chromosomes, n_alt_chromosomes)+1);
			}

		return null;
	}

	private void PrintStats()
	{
		stats_output.printf("n_reference_chromosomes n_variant_chromosomes n_sites n_found fraction_found\n");
		for (int i = 0; i < sample_names.size()*2; i++)
		{
			for (int j = 0; j < sample_names.size()*2; j++)
			{
				int N = (int)n_variants.get(i,j);
				int found = (int)n_found.get(i,j);

				if (N == 0) { continue; }
				if (found == 0) { continue; }

				double fraction_found = 100.0 * (double)found / (double)N;
				n_variants.set(i,j,0);
				n_found.set(i,j,0);
				stats_output.printf("%d %d %d %d %f\n", 
										i,
										j,
										N,
										found,
										fraction_found);
			}
		}
	}

    public void onTraversalDone(String sum) 
	{
		PrintStats();
		stats_output.flush();
		stats_output.close();
		out.println("MultiSampleCallerAccuracyTest done.");
		return;
	}

    public String reduceInit() 
	{
		return super.reduceInit();
	}

    public String reduce(MultiSampleCallResult record, String sum) 
	{
		return super.reduce(record, sum);
	}

	// END Walker Interface Functions
	/////////


	/////////
	// BEGIN Utility Functions

	// Filter a locus context by sample IDs
	//   (pulls out only reads from the specified samples, and returns them in one context).
    private AlignmentContext filterAlignmentContextBySamples(AlignmentContext context, List<String> sample_names)
    {
		HashSet<String> index = new HashSet<String>();
		for (int i = 0; i < sample_names.size(); i++)
		{
			index.add(sample_names.get(i));
		}

		ArrayList<SAMRecord> reads = new ArrayList();
		ArrayList<Integer> offsets = new ArrayList();

        for (int i = 0; i < context.getReads().size(); i++)
        {
            SAMRecord read = context.getReads().get(i);
            Integer offset = context.getOffsets().get(i);
            String RG = (String)(read.getAttribute("RG"));

            assert(header != null);
            assert(header.getReadGroup(RG) != null);

            String sample = header.getReadGroup(RG).getSample();
			if (SAMPLE_NAME_REGEX != null) { sample = sample.replaceAll(SAMPLE_NAME_REGEX, "$1"); }

			if (index.contains(sample))
			{
            	reads.add(read); 
            	offsets.add(offset); 
			}
        }

		return new AlignmentContext(context.getLocation(), reads, offsets);
    }

	// END Utility Functions
	/////////

}
