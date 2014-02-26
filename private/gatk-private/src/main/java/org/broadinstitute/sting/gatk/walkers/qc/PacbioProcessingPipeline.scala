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

package org.broadinstitute.sting.gatk.walkers.qc

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.picard.{ReorderSam, SortSam, AddOrReplaceReadGroups}
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 4/20/11
 * Time: 16:29 PM
 */


class PacbioProcessingPipeline extends QScript {

  @Input(doc="input FASTA/FASTQ/BAM file - or list of FASTA/FASTQ/BAM files. ", shortName="i", required=true)
  var input: File = _

  @Input(doc="Reference fasta file", shortName="R", required=true)
  var reference: File = _

  @Input(doc="dbsnp VCF file to use ", shortName="D", required=true)
  var dbSNP: File = _

  @Argument(doc="Number of jobs to scatter/gather. Default: 0." , shortName = "sg", required=false)
  var threads: Int = 0

  @Argument(doc="Sample Name to fill in the Read Group information (only necessary if using fasta/fastq)" , shortName = "sn", required=false)
  var sample: String = "NA"

  @Input(doc="The path to the binary of bwa to align fasta/fastq files", fullName="path_to_bwa", shortName="bwa", required=false)
  var bwaPath: File = _

  @Argument(doc="Input is a BLASR generated BAM file", shortName = "blasr", fullName="blasr_bam", required=false)
  var BLASR_BAM: Boolean = false

  @Hidden
  @Argument(doc="The default base qualities to use before recalibration. Default is Q20 (should be good for every dataset)." , shortName = "dbq", required=false)
  var dbq: Int = 20

  @Hidden
  @Argument(shortName="bwastring", required=false)
  var bwastring: String = ""

  @Hidden
  @Argument(shortName = "test", fullName = "test_mode", required = false)
  var testMode: Boolean = false

  val queueLogDir: String = ".qlog/"

  def script() {

    val fileList: Seq[File] = QScriptUtils.createSeqFromFile(input)

    for (file: File <- fileList) {

      var USE_BWA: Boolean = false
      var resetQuals: Boolean = true

      if (file.endsWith(".fasta") || file.endsWith(".fq") || file.endsWith(".fastq")) {
        if (bwaPath == null) {
          throw new UserException("You provided a fasta/fastq file but didn't provide the path for BWA");
        }
        USE_BWA = true
        if (file.endsWith(".fq") || file.endsWith(".fastq"))
          resetQuals = false
      }

      // FASTA -> BAM steps
      val alignedSam: File = file.getName + ".aligned.sam"
      val sortedBam: File = swapExt(alignedSam, ".sam", ".bam")
      val rgBam: File = swapExt(file, ".bam", ".rg.bam")

      val bamBase = if (USE_BWA) {rgBam} else {file}

      // BAM Steps
      val mqBAM: File = swapExt(bamBase, ".bam", ".mq.bam")
      val recalFile1: File = swapExt(bamBase, ".bam", ".recal1.table")
      val recalFile2: File = swapExt(bamBase, ".bam", ".recal2.table")
      val recalBam: File   = swapExt(bamBase, ".bam", ".recal.bam")
      val path1: String    = recalBam + ".before"
      val path2: String    = recalBam + ".after"

      if (USE_BWA) {
        add(align(file, alignedSam),
            sortSam(alignedSam, sortedBam),
            addReadGroup(sortedBam, rgBam, sample))
      }

      else if (BLASR_BAM) {
        val reorderedBAM = swapExt(bamBase, ".bam", ".reordered.bam")
        add(reorder(bamBase, reorderedBAM),
            changeMQ(reorderedBAM, mqBAM))
      }

      val bam = if (BLASR_BAM) {mqBAM} else {bamBase}

      add(cov(bam, recalFile1, resetQuals),
          recal(bam, recalFile1, recalBam),
          cov(recalBam, recalFile2, false))
    }
  }


  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 4
    this.isIntermediate = true
  }

  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = reference
  }


  case class align(@Input inFastq: File, @Output outSam: File) extends ExternalCommonArgs {
    def commandLine = bwaPath + " bwasw -b5 -q2 -r1 -z20 -t16 " + reference + " " + inFastq + " > " + outSam
    this.memoryLimit = 8
    this.analysisName = queueLogDir + outSam + ".bwa_sam_se"
    this.jobName = queueLogDir + outSam + ".bwa_sam_se"
  }

  case class sortSam (inSam: File, outBam: File) extends SortSam with ExternalCommonArgs {
    this.input = List(inSam)
    this.output = outBam
    this.memoryLimit = 8
    this.sortOrder = SortOrder.coordinate
    this.analysisName = queueLogDir + outBam + ".sortSam"
    this.jobName = queueLogDir + outBam + ".sortSam"
  }

  case class reorder (inSam: File, outSam: File) extends ReorderSam with ExternalCommonArgs {
    this.input = List(inSam)
    this.output = outSam
    this.sortReference = reference
  }

  case class changeMQ(inBam: File, outBam: File) extends PrintReads with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.out = outBam
    this.read_filter :+= "ReassignMappingQuality"
  }

  case class addReadGroup (inBam: File, outBam: File, sample: String) extends AddOrReplaceReadGroups with ExternalCommonArgs {
    @Output(doc="output bai file") var bai = swapExt(outBam, ".bam", ".bai")
    this.input = List(inBam)
    this.output = outBam
    this.RGID = "1"
    this.RGCN = "BI"
    this.RGPL = "PacBio_RS"
    this.RGSM = sample
    this.RGLB = "default_library"
    this.RGPU = "default_pu"
    this.analysisName = queueLogDir + outBam + ".rg"
    this.jobName = queueLogDir + outBam + ".rg"
  }

  case class cov (inBam: File, outRecalFile: File, resetQuals: Boolean) extends BaseRecalibrator with CommandLineGATKArgs {
    if (resetQuals) 
      this.DBQ = dbq
    this.knownSites :+= dbSNP
    this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "ContextCovariate")
    this.input_file :+= inBam
    this.disable_indel_quals = true
    this.out = outRecalFile
    this.analysisName = queueLogDir + outRecalFile + ".covariates"
    this.jobName = queueLogDir + outRecalFile + ".covariates"
    this.scatterCount = threads
    this.read_filter :+= "BadCigar"
  }

  case class recal (inBam: File, inRecalFile: File, outBam: File) extends PrintReads with CommandLineGATKArgs {
    this.DBQ = dbq
    this.input_file :+= inBam
    this.BQSR = inRecalFile
    this.out = outBam
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
    this.read_filter :+= "BadCigar"
    this.scatterCount = threads
  }
}
