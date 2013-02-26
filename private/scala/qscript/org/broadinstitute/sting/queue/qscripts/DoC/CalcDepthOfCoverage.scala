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

package org.broadinstitute.sting.queue.qscripts.DoC

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.util.VCF_BAM_utilities
import org.broadinstitute.sting.queue.util.DoC._
import org.broadinstitute.sting.gatk.walkers.coverage.CoverageUtils

class CalcDepthOfCoverage extends QScript {
  qscript =>

  @Input(doc = "bam input, as .bam or as a list of files", shortName = "I", required = true)
  var bams: File = _

  @Input(doc = "gatk jar file", shortName = "J", required = true)
  var gatkJarFile: File = _

  @Input(doc = "xhmm executable file", shortName = "xhmmExec", required = true)
  var xhmmExec: File = _

  @Input(shortName = "R", doc = "ref", required = true)
  var referenceFile: File = _

  @Input(shortName = "L", doc = "Intervals", required = false)
  var intervals: File = _

  @Argument(doc = "Expand the intervals by this number of bases (on each side)", shortName = "expandIntervals", required = false)
  var expandIntervals = 0

  @Argument(doc = "level of parallelism for BAM DoC.   By default is set to 0 [no scattering].", shortName = "scatter", required = false)
  var scatterCountInput = 0

  @Argument(doc = "Samples to run together for DoC.   By default is set to 1 [one job per sample].", shortName = "samplesPerJob", required = false)
  var samplesPerJob = 1

  @Output(doc = "Base name for files to output", shortName = "o", required = true)
  var outputBase: File = _

  @Hidden
  @Argument(doc = "How should overlapping reads from the same fragment be handled?", shortName = "countType", required = false)
  //
  // TODO: change this to be the default once UnifiedGenotyper output is annotated with *fragment-based depth* (and this is used for filtering):
  //var countType = CoverageUtils.CountPileupType.COUNT_FRAGMENTS_REQUIRE_SAME_BASE
  //
  var countType = CoverageUtils.CountPileupType.COUNT_READS

  @Argument(doc = "Maximum depth (before GATK down-sampling kicks in...)", shortName = "MAX_DEPTH", required = false)
  var MAX_DEPTH = 20000

  @Hidden
  @Argument(doc = "Number of read-depth bins", shortName = "NUM_BINS", required = false)
  var NUM_BINS = 200

  @Hidden
  @Argument(doc = "Starting value of read-depth bins", shortName = "START_BIN", required = false)
  var START_BIN = 1

  @Argument(doc = "Minimum read mapping quality", shortName = "MMQ", required = false)
  var minMappingQuality = 0

  @Argument(doc = "Minimum base quality to be counted in depth", shortName = "MBQ", required = false)
  var minBaseQuality = 0

  @Argument(doc = "Memory (in GB) required for storing the whole matrix in memory", shortName = "wholeMatrixMemory", required = false)
  var wholeMatrixMemory = -1

  @Argument(doc = "Memory (in GB) required for merging the big matrix", shortName = "largeMemory", required = false)
  var largeMemory = -1

  @Argument(shortName = "sampleIDsMap", doc = "File mapping BAM sample IDs to desired sample IDs", required = false)
  var sampleIDsMap: String = ""

  @Argument(shortName = "sampleIDsMapFromColumn", doc = "Column number of OLD sample IDs to map", required = false)
  var sampleIDsMapFromColumn = 1

  @Argument(shortName = "sampleIDsMapToColumn", doc = "Column number of NEW sample IDs to map", required = false)
  var sampleIDsMapToColumn = 2

  @Argument(shortName = "longJobQueue", doc = "Job queue to run the 'long-running' commands", required = false)
  var longJobQueue: String = ""

  @Argument(shortName = "minCoverageCalcs", doc = "Coverage thresholds against which to measure each sample", required = false)
  var minCoverageCalcs = List(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60).toArray


  val PREPARED_TARGS_SUFFIX: String = ".merged.interval_list"

  val DOC_OUTPUT_SUFFIX: String = ".DoC.txt"
  val BASE_DOC_OUTPUT_SUFFIX: String = ".per_base.DoC.txt"


  trait LargeMemoryLimit extends CommandLineFunction {
    if (largeMemory < 0) {
      this.memoryLimit = 6
    }
    else {
      this.memoryLimit = largeMemory
    }
  }

  trait WholeMatrixMemoryLimit extends CommandLineFunction {
    // Since loading ALL of the data can take significant memory:
    if (wholeMatrixMemory < 0) {
      this.memoryLimit = 24
    }
    else {
      this.memoryLimit = wholeMatrixMemory
    }
  }

  trait LongRunTime extends CommandLineFunction {
    if (longJobQueue != "")
      this.jobQueue = longJobQueue
  }

  class extendFlanks(inIntervals: File, outIntervals: File) extends WriteFlankingIntervalsFunction {
      this.inputIntervals = inIntervals
      this.outputIntervals = outIntervals
      this.flankSize = qscript.expandIntervals
      this.reference = qscript.referenceFile
  }

  def script = {
    val prepTargetsInitial = new PrepareTargets(List(qscript.intervals), outputBase.getPath + PREPARED_TARGS_SUFFIX, xhmmExec, referenceFile)
    add(prepTargetsInitial)

    var prepTargets = prepTargetsInitial
    if (qscript.intervals != null && qscript.expandIntervals > 0) {
      val flankIntervalsFile = new File(outputBase.getPath + ".flank.interval_list")
      val extendCommand = new extendFlanks(prepTargetsInitial.out, flankIntervalsFile)
      add(extendCommand)

      prepTargets = new PrepareTargets(List(qscript.intervals, extendCommand.outputIntervals), outputBase.getPath + ".withFlanking" + PREPARED_TARGS_SUFFIX, xhmmExec, referenceFile)
      add(prepTargets)
    }

    trait CommandLineGATKArgs extends CommandLineGATK {
      this.intervals :+= prepTargets.out
      this.jarFile = qscript.gatkJarFile
      this.reference_sequence = qscript.referenceFile
      this.logging_level = "INFO"
    }

    val sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = VCF_BAM_utilities.getMapOfBAMsForSample(VCF_BAM_utilities.parseBAMsInput(bams))
    val samples: List[String] = sampleToBams.keys.toList
    Console.out.printf("Samples are %s%n", samples)

    val groups: List[Group] = buildDoCgroups(samples, sampleToBams, samplesPerJob, outputBase)
    var docs: List[DoCwithDepthOutputAtEachBase] = List[DoCwithDepthOutputAtEachBase]()
    for (group <- groups) {
      Console.out.printf("Group is %s%n", group)
      docs ::= new DoCwithDepthOutputAtEachBase(group.bams, group.DoC_output, countType, MAX_DEPTH, minMappingQuality, minBaseQuality, scatterCountInput, START_BIN, NUM_BINS, minCoverageCalcs) with CommandLineGATKArgs
    }
    addAll(docs)

    for (minCoverageCalc <- minCoverageCalcs) {
      val mergeDepths = new MergeGATKdepths(docs.map(u => u.intervalSampleOut), outputBase.getPath + ".min_" + minCoverageCalc + DOC_OUTPUT_SUFFIX, "_%_above_" + minCoverageCalc, xhmmExec, sampleIDsMap, sampleIDsMapFromColumn, sampleIDsMapToColumn, None, false) with WholeMatrixMemoryLimit
      add(mergeDepths)
    }

    // Want 0 precision for base counts:
    val mergeBaseDepths = new MergeGATKdepths(docs.map(u => u.outPrefix), outputBase.getPath + BASE_DOC_OUTPUT_SUFFIX, "Depth_for_", xhmmExec, sampleIDsMap, sampleIDsMapFromColumn, sampleIDsMapToColumn, Some(0), true) with LargeMemoryLimit with LongRunTime
    add(mergeBaseDepths)
  }
}
