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

package org.broadinstitute.sting.queue.qscripts.qc

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.interval.IntervalSetRule

class DiagnoseIntervals extends QScript {
  // Required arguments
  @Argument(shortName = "I", required = true, doc = "List of BAM files") var bamList: List[File] = _
  @Argument(shortName = "L", required = true, doc = "Intervals file") var intervalsFile: List[File] = _
  @Argument(shortName = "R", required = false, doc = "Reference sequence") var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")

  // Script arguments
  @Argument(shortName = "sc", required = false, doc = "Scatter count") var jobs: Int = 50
  @Argument(shortName = "isr", required = false, doc = "interval set rule") var intervalSetRule = IntervalSetRule.UNION
  @Argument(shortName = "ds", required = false, doc = "downsampling fractions") var downsamplingFraction = 1.0
  @Argument(shortName = "mem", required = false, doc = "memory for the jvm") var memoryLimit = 4

  // Walker optional arguments
  @Argument(shortName = "BQ", doc = "The minimum Base Quality that is considered for calls", required = false) var minimumBaseQuality: Int = 20
  @Argument(shortName = "MQ", doc = "The minimum read mapping quality considered for calls", required = false) var minimumMappingQuality: Int = 20
  @Argument(shortName = "min", doc = "The minimum allowable coverage, used for calling LOW_COVERAGE", required = false) var minimumCoverage: Int = 5
  @Argument(shortName = "max", doc = "The maximum allowable coverage, used for calling EXCESSIVE_COVERAGE", required = false) var maximumCoverage: Int = 700
  @Argument(shortName = "ins", doc = "The maximum allowed distance between a read and its mate", required = false) var maximumInsertSize: Int = 500
  @Argument(shortName = "stV", doc = "The needed proportion of samples containing a call for the interval to adopt the call ", required = false) var votePercentageThreshold: Double = 0.50
  @Argument(shortName = "stBM", doc = "The proportion of the loci needed for calling BAD_MATE", required = false) var badMateStatusThreshold: Double = 0.50
  @Argument(shortName = "stC", doc = "The proportion of the loci needed for calling LOW_COVERAGE and COVERAGE_GAPS", required = false) var coverageStatusThreshold: Double = 0.20
  @Argument(shortName = "stXC", doc = "The proportion of the loci needed for calling EXCESSIVE_COVERAGE", required = false) var excessiveCoverageThreshold: Double = 0.20
  @Argument(shortName = "stQ", doc = "The proportion of the loci needed for calling POOR_QUALITY", required = false) var qualityStatusThreshold: Double = 0.50
  @Argument(shortName = "m", required = false, doc = "missing intervals") var missingIntervals:File = _

  def script {

    val dt_output = swapExt(bamList(0), ".interval_list", ".diagnose.vcf")

    val dt = new DiagnoseTargets()
    // Required arguments
    dt.input_file = bamList
    dt.intervalsString = intervalsFile
    dt.reference_sequence = referenceFile

    // Script arguments
    dt.scatterCount = jobs
    dt.interval_set_rule = intervalSetRule
    dt.dfrac = downsamplingFraction
    dt.memoryLimit = memoryLimit
    dt.out = dt_output

    // Walker optional arguments
    dt.minimum_base_quality = minimumBaseQuality
    dt.minimum_mapping_quality = minimumMappingQuality
    dt.minimum_coverage = minimumCoverage
    dt.maximum_coverage = maximumCoverage
    dt.maximum_insert_size = maximumInsertSize
    dt.voting_status_threshold = votePercentageThreshold
    dt.bad_mate_status_threshold = badMateStatusThreshold
    dt.coverage_status_threshold = coverageStatusThreshold
    dt.excessive_coverage_status_threshold = excessiveCoverageThreshold
    dt.quality_status_threshold = qualityStatusThreshold
    dt.missing = missingIntervals

    add(dt)
  }

}


