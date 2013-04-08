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

package org.broadinstitute.sting.queue.qscripts.DataProcessing

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode

import net.sf.samtools.SAMFileHeader.SortOrder

import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.utils.NGSPlatform
import org.broadinstitute.sting.queue.extensions.picard._
import io.Source
import org.broadinstitute.sting.gatk.walkers.genotyper.{UnifiedGenotyperEngine, GenotypeLikelihoodsCalculationModel}
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.gatk.downsampling.DownsampleType

class AncientDNAProcessingPipeline extends QScript {
  qscript =>

  /** ***************************************************************************
    * Required Parameters
    * ***************************************************************************/

  @Input(doc = "a table with all the necessary information to process the data", fullName = "metadata", shortName = "m", required = false, exclusiveOf = "individual")
  var metaData: File = _

  /** ******************************************************************************
    * Additional Parameters that the pipeline should have pre-defined in the image
    * ******************************************************************************/

  @Argument(doc="Reference fasta file", fullName="reference", shortName="R", required=false)
  var reference: File = new File("/groups/reich/reference-genomes/human_hg19/human_g1k_v37/human_g1k_v37.fasta")

  @Argument(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=false)
  var dbSNP: Seq[File] = Seq(new File("/groups/reich/sw/gatk/GenomeAnalysisTK-2.3-9/bundle/dbsnp_137.b37.vcf"))

  @Argument(doc = "The path to the binary of bwa", fullName = "path_to_bwa", shortName = "bwa", required = false)
  var bwaPath: File = new File("/groups/reich/sw/bwa-0.6.1/bwa")

  @Argument(doc = "The path to the binary of samtools ", fullName = "path_to_samtools", shortName = "samtools", required = false)
  var samtoolsPath: File = new File("/groups/reich/sw/samtools-0.1.18/samtools")

  @Argument(doc = "extra VCF files to use as reference indels for Indel Realignment", fullName = "extra_indels", shortName = "indels", required = false)
  var indelSites: Seq[File] = Seq()

  @Argument(doc = "short job queue for LSF", fullName = "queue", shortName = "queue", required = false)
  var queue: String = "short"

  @Argument(doc = "long job queue for LSF", fullName = "longQueue", shortName = "longQueue", required = false)
  var longQueue: String = "long"

  @Argument(doc = "long job queue for LSF", fullName = "miniQueue", shortName = "miniQueue", required = false)
  var miniQueue: String = "mini"

  @Argument(doc = "project for LSF", fullName = "project", shortName = "project", required = false)
  var project: String = "default"

  @Argument(doc = "tmp dir", fullName = "tmpDir", shortName = "tmpDir", required = false)
  var tmpDir: String = "/scratch/gd73/tmp/"

  @Argument(doc = "runname", fullName = "runname", shortName = "runname", required = false)
  var runname: String = "default"

  @Argument(doc = "Start from a BAM and don't do mapping", fullName = "startFromBAM", shortName = "startFromBAM", required = false)
  var startFromBAM: Boolean = false
  @Argument(doc = "Input bam in case startFromBAM is set", fullName = "inputBAM", shortName = "inputBAM", required = false)
  var inputBAM:File = _

  @Argument(doc = "If set, outputs will be split by chromosome", fullName = "splitByContig", shortName = "splitByContig", required = false)
  var splitByContig: Boolean = false


  @Hidden
  @Argument(doc = "Skip BWA aln", fullName = "skipAln", shortName = "skipAln", required = false)
  var skipAln: Boolean = false

  @Input(doc = "Interval file with targets used in exome capture (used for QC metrics)", fullName = "targets", shortName = "targets", required = false)
  var targets: File = _

  @Input(doc = "Interval file with baits used in exome capture (used for QC metrics)", fullName = "baits", shortName = "baits", required = false)
  var baits: File = _

  @Hidden
  @Argument(doc = "Define the default platform for Count Covariates -- useful for techdev purposes only.", fullName = "platform", shortName = "platform", required = false)
  var platform: String = "illumina"

  @Hidden
  @Argument(doc = "Define the default platform for Count Covariates -- useful for techdev purposes only.", fullName = "center", shortName = "center", required = false)
  var center: String = "BI"

  @Argument(doc = "Adaptor sequence", fullName = "adaptor", shortName = "adaptor", required = false)
  var adaptorSequences: Seq[String] = _

  @Argument(doc = "remove Adaptors and merge read pairs", fullName = "mergePairs", shortName = "merge", required = false)
  var mergeReadPairs: Boolean = false


  /** **************************************************************************
    * Hidden Parameters
    * ***************************************************************************/
  @Hidden
  @Argument(doc = "Use BWASW instead of BWA aln", fullName = "use_bwa_sw", shortName = "bwasw", required = false)
  var useBWAsw: Boolean = false

  @Hidden
  @Argument(doc = "Collect Picard QC metrics", fullName = "skipQC", shortName = "skipqc", required = false)
  var skipQCMetrics: Boolean = false

  @Hidden
  @Argument(doc = "Number of threads jobs should use when possible", fullName = "numThreads", shortName = "nct", required = false)
  var numThreads: Int = 1

  @Hidden
  @Argument(doc = "Number of threads jobs should use when possible", fullName = "bwathreads", shortName = "bwathreads", required = false)
  var bwaThreads: Int = 1

  @Hidden
  @Argument(doc = "Default memory limit per job", fullName = "mem_limit", shortName = "mem", required = false)
  var memLimit: Int = 2

  @Hidden
  @Argument(doc = "How many ways to scatter/gather", fullName = "scatter_gather", shortName = "sg", required = false)
  var nContigs: Int = 0

  @Hidden
  @Argument(doc = "Define the default platform for Count Covariates -- useful for techdev purposes only.", fullName = "default_platform", shortName = "dp", required = false)
  var defaultPlatform: String = ""

  @Hidden
  @Argument(doc = "Run the pipeline in test mode only", fullName = "test_mode", shortName = "test", required = false)
  var testMode: Boolean = false

  @Hidden
  @Argument(doc = "Run the pipeline in quick mode only", fullName = "quick", shortName = "quick", required = false)
  var quick: Boolean = false

  @Hidden
  @Argument(doc = "Base path for FASTQs", fullName = "baseFastqPath", shortName = "baseFastqPath", required = false)
  val baseFastqPath: String = ""

  @Hidden
  @Argument(doc = "Run single sample germline calling in resulting bam", fullName = "doSingleSampleCalling", shortName = "call", required = false)
  var doSingleSampleCalling: Boolean = false

  @Hidden
  @Argument(doc = "Do post-recalibration to get BQSR statistics", fullName = "doPostRecal", shortName = "postRecal", required = false)
  var doPostRecal: Boolean = false

  @Hidden
  @Argument(doc = "BWA Parameteres", fullName = "bwa_parameters", shortName = "bp", required = false)
  val bwaParameters: String = " -o 2 -n 0.01 -l 16500 " // " -q 5 -l 32 -k 2 -o 1 "

  @Hidden
  @Argument(doc = "Base path for Picard executables", fullName = "picardBase", shortName = "picardBase", required = false)
  val picardBase: String = "/groups/reich/sw/picard_1.56/picard-tools-1.56/"

  @Hidden
  @Argument(doc = "Base path for SeqPrep executables", fullName = "seqPrepCMD", shortName = "seqprep", required = false)
  val seqPrepCMD: String = "/humgen/gsa-hpprojects/dev/delangel/SeqPrep/SeqPrep"

  @Hidden
  @Argument(doc = "Do post-recalibration to get BQSR statistics", fullName = "doSeqPrep", shortName = "doSeqPrep", required = false)
  var doSeqPrep: Boolean = false

  val cleaningExtension: String = ".clean.bam"
  val headerVersion: String = "#FILE1,FILE2,RG,SAMPLE,LIBRARY,PLATFORM,PLATFORM_UNIT,CENTER,DESCRIPTION,DATE_SEQUENCED"


  /** **************************************************************************
    * Main script
    * ***************************************************************************/

  val contigs:List[String] = List("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT")
  def script() {


    val lanes: Seq[MetaInfo] =
    {
      var lanesFromFile = Seq[MetaInfo]()
      var id = 1
      for (line: String <- Source.fromFile(metaData).getLines()) {
        if (line.startsWith("#")) {
          checkMetaDataHeader(line)
        }
        else {
          lanesFromFile :+= new MetaInfo(id, line.split(","))
          id = id + 1
        }
      }
      lanesFromFile
    }
    val commonName = lanes.last.commonName
    var sampleBAM = new File(tmpDir+ commonName + ".bam")
    var laneBAMs = Seq[File]()

    if (startFromBAM) {
      laneBAMs :+= inputBAM
      if(!inputBAM.endsWith("bam.list")) sampleBAM = inputBAM

    } else {
       for (lane <- lanes) {
        laneBAMs +:=  performAlignment(lane)
      }
//      add(joinBAMsFast(laneBAMs, sampleBAM,""))
    }


    clean(laneBAMs, sampleBAM)

    val preRecalFile = swapExt(tmpDir, sampleBAM, ".bam",".pre_recal.table")
    val postRecalFile = swapExt(tmpDir, sampleBAM, ".bam", ".post_recal.table")
    val recalBAM = swapExt(tmpDir, sampleBAM, ".bam", ".clean.dedup.recal.bam")
    val outVCF = swapExt(tmpDir, recalBAM, ".bam", ".vcf")
    if (splitByContig) {
      var chrBAMs = Seq[File]()
      for (chr <- contigs) {
        val cleanBAM = swapExt(tmpDir, sampleBAM, ".bam", ".chr"+chr+cleaningExtension)
        val dedupBAM = swapExt(tmpDir, sampleBAM, ".bam", ".chr"+chr+ ".clean.dedup.bam")
        val duplicateMetricsFile = swapExt(tmpDir, sampleBAM, ".bam", ".chr"+chr+ ".duplicateMetrics")

        add(dedup(cleanBAM, dedupBAM, duplicateMetricsFile))
        chrBAMs :+= dedupBAM
      }
      recalibrate(chrBAMs, preRecalFile, postRecalFile, recalBAM)
    } else {
      val cleanBAM = swapExt(tmpDir, sampleBAM, ".bam", cleaningExtension)
      val dedupBAM = swapExt(tmpDir, sampleBAM, ".bam", ".clean.dedup.bam")
      val recalBAM = swapExt(tmpDir, sampleBAM, ".bam", ".clean.dedup.recal.bam")
      val duplicateMetricsFile = swapExt(tmpDir, sampleBAM, ".bam", ".duplicateMetrics")

      add(dedup(cleanBAM, dedupBAM, duplicateMetricsFile))

      recalibrate(Seq(dedupBAM), preRecalFile, postRecalFile, recalBAM)
    }

    if (!qscript.skipQCMetrics) {
      if (qscript.targets != null && qscript.baits != null) {
        add(calculateHSMetrics(recalBAM, swapExt(recalBAM, ".bam", ".hs_metrics")))
      }
      // collect QC metrics based on full BAM
      val outGcBiasMetrics = swapExt(recalBAM, ".bam", ".gc_metrics")
      val outMultipleMetrics = swapExt(recalBAM, ".bam", ".multipleMetrics")

      add(calculateGCMetrics(recalBAM, outGcBiasMetrics))
      add(calculateMultipleMetrics(recalBAM, outMultipleMetrics))

    }


    add(call(recalBAM, outVCF))

  }

  /** **************************************************************************
    * Helper classes and methods
    * ***************************************************************************/

  private class MetaInfo( val id: Int,
                          val file1: File,
                          val file2: File,
                          val readGroup: String,
                          val sample: String,
                          val library: String,
                          val platform: NGSPlatform,
                          val platformUnit: String,
                          val center: String,
                          val description: String,
                          val dateSequenced: String
                          ) {


    def this(idP: Int, headerArray: Array[String]) =
      this(
        idP,
        new File(headerArray(0)),
        if (headerArray(1).isEmpty) {
          null
        } else {
          new File(headerArray(1))
        },
        headerArray(2),     //rg
        headerArray(3),     //sm
        headerArray(4),     //lib
        NGSPlatform.fromReadGroupPL(headerArray(5)),        //pl
        headerArray(6),     //pu
        headerArray(7),     //cn
        headerArray(8),     //desc
        headerArray(9)      //date
      )

    def bamFileName = project+"."+runname + "." + sample +"." +id + ".bam"
    def commonName = project+"."+runname + "." + sample
    def readGroupString = "@RG\tID:%s\tCN:%s\tDS:%s\tDT:%s\tLB:%s\tPL:%s\tPU:%s\tSM:%s".format(readGroup, center, description, dateSequenced, library, platform, platformUnit, sample)
  }

  def checkMetaDataHeader(header: String) {
    assert(header == headerVersion,
      String.format("Your header doesn't match the header this version of the pipeline is expecting.\n\tYour header: %s\n\t Our header: %s\n", header, headerVersion))
  }


  /**
   * BWA alignment for the lane (pair ended or not)
   *
   * @return an aligned bam file for the lane
   */
  def performAlignment(metaInfo: MetaInfo): File = {
    val alnBAM: File = new File(tmpDir+metaInfo.bamFileName)
    val alnSAM: File = swapExt(tmpDir,alnBAM,".bam",".aligned.sam")

    // check if we need to merge read pairs and strip adaptor sequence first
    if (mergeReadPairs) {
      // merge pairs, then split between merged and non-merged reads
      val mergeOut: File = swapExt(tmpDir, alnBAM, "bam","mergeOut.bam")
      val saiFile1: File = swapExt(tmpDir, mergeOut, ".bam",".1.sai")
      val saiFile2: File = swapExt(tmpDir, mergeOut, ".bam",".2.sai")
      if (metaInfo.file1.endsWith("bam") ) {
        add(mergePairs(metaInfo.file1, mergeOut))
      } else {
        // want to read merge pairs but input file is in fastq format:
        // first convert to BAM and index
        val unmergedBAM: File = swapExt(tmpDir,alnBAM,".bam", ".unmerged.bam")
        val unmergedBAMIdx: File = swapExt(tmpDir, unmergedBAM, ".bam",".bam.bai")
        add(fastqToSam(metaInfo.file1,metaInfo.file2,unmergedBAM))
        add(indexBAM(unmergedBAM,unmergedBAMIdx))
        add(mergePairs(unmergedBAM,mergeOut))
      }

      // BWA ALN
      if (!skipAln){
        add(bwa_aln_pe(mergeOut,saiFile1, 1))
        add(bwa_aln_pe(mergeOut,saiFile2, 2))
      }
      //BWA SAMPE
      add(bwa_sam_pe_bam(mergeOut, saiFile1, saiFile2, alnSAM, metaInfo.readGroupString))

    } else {
      val saiFile1: File = swapExt(tmpDir,alnBAM,".bam",".1.sai")
      val saiFile2: File = swapExt(tmpDir,alnBAM,".bam",".2.sai")
      // no read pair merging/trimming
      if (metaInfo.file1.endsWith("bam")){
        // no read pair merging, bam file as input:
        if (!skipAln) {
          add(bwa_aln_pe(metaInfo.file1, saiFile1,1))
          add(bwa_aln_pe(metaInfo.file1, saiFile2,2))
        }
        add(bwa_sam_pe_bam(metaInfo.file1, saiFile1, saiFile2, alnSAM, metaInfo.readGroupString))
      } else {
        // fastq input:
        if (!skipAln) {
          add(bwa_aln(metaInfo.file1, saiFile1))
          if (metaInfo.file2 != null) {
            add(bwa_aln(metaInfo.file2, saiFile2))
            // no second pair: do single-ended alignment
            add(bwa_sam_se(metaInfo.file1, saiFile1, alnSAM, metaInfo.readGroupString))

        }
        add(bwa_sam_pe(metaInfo.file1, metaInfo.file2, saiFile1, saiFile2, alnSAM, metaInfo.readGroupString))
        }
      }

    }

    add(sortSam(alnSAM, alnBAM, SortOrder.coordinate))

    alnBAM
  }

  def clean(inputFiles: Seq[File], baseName: File) {
    if (splitByContig) {
      for (chr <- contigs) {
        val targetIntervals = swapExt(tmpDir,baseName, ".bam", ".chr"+chr+ ".cleaning.interval_list")
        add(target(inputFiles,baseName, targetIntervals, chr))
        add(indel(inputFiles,baseName, targetIntervals, chr))

      }
    }
    else {
      val targetIntervals = swapExt(tmpDir,baseName, ".bam", ".cleaning.interval_list")
      add(target(inputFiles,baseName, targetIntervals,""), indel(inputFiles,baseName, targetIntervals,""))
    }
  }

  def recalibrate(dedupBAMs: Seq[File], preRecalFile: File, postRecalFile: File, recalBAM: File) {
    add(bqsr(dedupBAMs, preRecalFile),
      apply_bqsr(dedupBAMs, preRecalFile, recalBAM))

    if (qscript.doPostRecal) {
      add(bqsr(Seq(recalBAM), postRecalFile))
    }
  }


  /** **************************************************************************
    * Classes (GATK Walkers)
    * ***************************************************************************/


  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = qscript.memLimit
    this.isIntermediate = true
    this.jobQueue = queue
  }

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
  }

  trait SAMargs extends PicardBamFunction with ExternalCommonArgs {
    this.maxRecordsInRam = 100000
  }

  case class target(inBAMs: Seq[File],baseName: File, outIntervals: File, chr: String) extends RealignerTargetCreator with CommandLineGATKArgs {
    this.input_file = inBAMs
    this.out = outIntervals
    this.mismatchFraction = 0.0
    this.known ++= qscript.dbSNP
    if (indelSites != null)
      this.known ++= qscript.indelSites
    this.analysisName = outIntervals + ".target"
    this.jobName = outIntervals + ".target"
    if (splitByContig) {
      this.scatterCount = 0
      this.intervalsString :+= chr
    }
    else {
      this.scatterCount = nContigs
      if (qscript.targets != null)
        this.intervals :+= qscript.targets
    }
    this.isIntermediate = false
  }

  case class indel(inBAMs: Seq[File],baseName: File, tIntervals: File, chr: String) extends IndelRealigner with CommandLineGATKArgs {
    var out1: File = swapExt(tmpDir, baseName, ".bam", cleaningExtension)

    if (splitByContig) {
      out1 = swapExt(tmpDir, baseName, ".bam", ".chr"+chr+cleaningExtension)
      this.intervalsString :+= chr
      this.scatterCount = 0
    }
    else {
      this.scatterCount = nContigs
      if (qscript.targets != null)
        this.intervals :+= qscript.targets
    }

    this.input_file = inBAMs
    this.targetIntervals = tIntervals


    this.o = out1

    this.known ++= qscript.dbSNP
    if (qscript.indelSites != null)
      this.known ++= qscript.indelSites
    this.consensusDeterminationModel = ConsensusDeterminationModel.USE_READS

    this.noPGTag = qscript.testMode
    this.analysisName = baseName.toString + "clean"
    this.jobName = baseName.toString + ".clean"
//    this.isIntermediate = false

  }

  case class bqsr(inBAMs: Seq[File], outRecalFile: File) extends BaseRecalibrator with CommandLineGATKArgs {
    this.knownSites ++= qscript.dbSNP
    this.covariate ++= Seq("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "ContextCovariate")
    this.input_file = inBAMs
    this.disable_indel_quals = true
    this.out = outRecalFile
    if (!defaultPlatform.isEmpty) this.default_platform = defaultPlatform
    this.scatterCount = nContigs
    this.analysisName = outRecalFile + ".covariates"
    this.jobName = outRecalFile + ".covariates"
    this.memoryLimit = Some(4) // needs 4 GB to store big tables in memory
    this.nct = Some(qscript.numThreads)
  }

  case class apply_bqsr(inBAMs: Seq[File], inRecalFile: File, outBAM: File) extends PrintReads with CommandLineGATKArgs {
    this.input_file = inBAMs
    this.BQSR = inRecalFile
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.out = outBAM
    this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = outBAM + ".recalibration"
    this.jobName = outBAM + ".recalibration"
    this.nct = Some(qscript.numThreads)
  }

  case class reduce(inBAM: File, outBAM: File) extends ReduceReads with CommandLineGATKArgs {
    this.input_file :+= inBAM
    this.out = outBAM
    this.isIntermediate = false
    this.analysisName = outBAM + ".reduce"
    this.jobName = outBAM + ".reduce"
    this.memoryLimit = Some(4)

  }

  case class call(inBAM: File, outVCF: File) extends UnifiedGenotyper with CommandLineGATKArgs {
    this.input_file :+= inBAM
    this.out = outVCF
    this.isIntermediate = false
    this.analysisName = outVCF + ".singleSampleCalling"
    this.jobName = outVCF + ".singleSampleCalling"
    this.dbsnp = qscript.dbSNP(0)
    this.downsample_to_coverage = 600
    this.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.scatterCount = nContigs
    this.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES

    if (qscript.targets != null)
      this.intervals :+= qscript.targets
  }

  case class samToFastQ(inBAM: File, out1: File, out2: File) extends SamToFastq with ExternalCommonArgs {
    this.input :+= inBAM
    this.fastq = out1
    this.secondEndFastQ = out2

  }

  case class fastqToSam(inFQ1: File, inFQ2: File, out: File) extends FastqToSam with ExternalCommonArgs {
    this.fastq=inFQ1
    this.secondEndFastQ=inFQ2
    this.bam=out
    //this.memoryLimit = 1
    this.analysisName = out + ".fq2Bam"
    this.jobName = out + ".fq2Bam"

  }
  case class mergePairs(inBAM: File, outBAM: File) extends ReadAdaptorTrimmer with CommandLineGATKArgs {
    this.input_file :+= inBAM
    this.out = outBAM
    this.analysisName = outBAM + ".mergePairs"
    this.jobName = outBAM + ".mergePairs"

  }
  case class joinBAMsFast(inBAMs: Seq[File], outBAM: File,contig:String) extends PrintReads with CommandLineGATKArgs {
    this.input_file = inBAMs
    this.o = outBAM

    this.dt = DownsampleType.NONE
    this.baq = CalculationMode.OFF
    this.analysisName = outBAM + contig+".joinBAMs"
    this.jobName = outBAM + contig+ ".joinBAMs"
    this.isIntermediate = false
    this.nct = Some(qscript.numThreads)
   // this.scatterCount = nContigs
  }

  case class joinBAMs(inBAMs: Seq[File], outBAM: File,contig:String) extends MergeSamFiles with ExternalCommonArgs {
    this.input = inBAMs
    this.output = outBAM
    this.analysisName = outBAM + ".joinBAMs"
    this.jobName = outBAM + ".joinBAMs"
   // this.compressionLevel = Some(0)
  }

  /** **************************************************************************
    * Classes (non-GATK programs)
    * ***************************************************************************/


  case class dedup(inBAM: File, outBAM: File, metricsFile: File) extends MarkDuplicates with ExternalCommonArgs {
    this.input :+= inBAM
    this.output = outBAM
    this.metrics = metricsFile
    //this.memoryLimit = 4
    this.analysisName = outBAM + ".dedup"
    this.jobName = outBAM + ".dedup"
    this.assumeSorted = Some(true)
//    this.isIntermediate = false
  }

  case class calculateHSMetrics(inBAM: File, outFile: File) extends CalculateHsMetrics with ExternalCommonArgs {
    @Output(doc = "Metrics output", required = false) var ouths: File = outFile
    this.isIntermediate = false
    this.reference = qscript.reference
    this.input :+= inBAM
    this.output = outFile
    if (qscript.targets != null)
      this.targets = qscript.targets
    this.baits = qscript.baits
    this.analysisName = outFile + ".hsMetrics"
    this.jobName = outFile + ".hsMetrics"
    this.jarFile = new File(qscript.picardBase + "CalculateHsMetrics.jar")
    // todo - do we want to compute per-read group HS metrics?

  }

  case class calculateGCMetrics(inBAM: File, outFile: File) extends CollectGcBiasMetrics with ExternalCommonArgs {
    @Output(doc = "Metrics output", required = false) var outgc: File = outFile
    this.reference = qscript.reference
    this.isIntermediate = false
    this.input :+= inBAM
    this.output = outFile
    this.analysisName = inBAM + ".gcMetrics"
    this.jobName = inBAM + ".gcMetrics"
    this.jarFile = new File(qscript.picardBase + "CollectGcBiasMetrics.jar")
  }

  case class calculateMultipleMetrics(inBAM: File, outFile: File) extends CollectMultipleMetrics with ExternalCommonArgs {
    @Output(doc = "Metrics output", required = false) var outmm: File = outFile
    this.reference = qscript.reference
    this.input :+= inBAM
    this.isIntermediate = false
    this.output = outFile
    this.analysisName = inBAM + ".multipleMetrics"
    this.jobName = inBAM + ".multipleMetrics"
    this.jarFile = new File(qscript.picardBase + "CollectMultipleMetrics.jar")
  }

  case class sortSam(inSam: File, outBAM: File, sortOrderP: SortOrder) extends SortSam with ExternalCommonArgs {
    this.input :+= inSam
    this.output = outBAM
    this.sortOrder = sortOrderP
    this.analysisName = outBAM + ".sortSam"
    this.jobName = outBAM + ".sortSam"
//    this.isIntermediate = false

    //    this.compressionLevel = Some(0)
  }

  case class validate(inBAM: File, outLog: File) extends ValidateSamFile with ExternalCommonArgs {
    this.input :+= inBAM
    this.output = outLog
    this.REFERENCE_SEQUENCE = qscript.reference
    this.isIntermediate = false
    this.analysisName = outLog + ".validate"
    this.jobName = outLog + ".validate"
  }

  case class revert(inBAM: File, outBAM: File, removeAlignmentInfo: Boolean) extends RevertSam with ExternalCommonArgs {
    this.output = outBAM
    this.input :+= inBAM
    this.removeAlignmentInformation = removeAlignmentInfo
    this.sortOrder = if (removeAlignmentInfo) {
      SortOrder.queryname
    } else {
      SortOrder.coordinate
    }
    this.analysisName = outBAM + "revert"
    this.jobName = outBAM + ".revert"
  }

  case class convertToFastQ(inBAM: File, outFQ: File) extends SamToFastq with ExternalCommonArgs {
    this.input :+= inBAM
    this.fastq = outFQ
    this.analysisName = outFQ + "convert_to_fastq"
    this.jobName = outFQ + ".convert_to_fastq"
  }

  case class sam2bam(inSAM: File,outBAM: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "sam file") var sam = inSAM
    @Output(doc = "output bam file") var bam = outBAM

    def commandLine = samtoolsPath + " view -bS " + " " + sam+ " > " + bam

    this.memoryLimit = 2
    this.analysisName = outBAM + ".sam2bam"
    this.jobName = outBAM + ".sam2bam"

  }

  case class indexBAM(inBAM: File, outFile: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "sam file") var bam = inBAM
    @Output(doc = "output bam file") var bai = outFile

    def commandLine = samtoolsPath + " index " + " " + bam

    this.memoryLimit = 2
    this.analysisName = inBAM + ".index"
    this.jobName = inBAM + ".index"
    this.jobQueue = miniQueue

  }

  case class FilterPairs(inBAM: File, outBAM: File, filter: Boolean) extends PrintReads with CommandLineGATKArgs with PairedRead {
    this.input_file :+= inBAM
    this.out= outBAM
    this.dt = DownsampleType.NONE
    this.baq = CalculationMode.OFF
    this.analysisName = outBAM+".filterPairedReads"
    this.jobName = outBAM + ".filterPairedReads"
//    this.isIntermediate = false
    this.nct = Some(qscript.numThreads)

   // this.rf :+= "PairedReadFilter"
    this.filter_paired = filter

  }
  case class seqPrep(inFQ1: File, inFQ2: File, outMerged: File, outUnmerged1: File, outUnmerged2: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "fq file 1 to be merged") var fq1 = inFQ1
    @Input(doc = "fq file 2 to be merged") var fq2 = inFQ2
    @Output(doc = "merged fq file") var merged = outMerged
    @Output(doc = "unmerged file 1") var unmerged1 = outUnmerged1
    @Output(doc = "unmerged file 1") var unmerged2 = outUnmerged2

    this.memoryLimit = 2
    this.analysisName = outMerged + ".seqprep"
    this.jobName = outMerged + ".seqprep"
    this.isIntermediate = false

    var cmd:String = seqPrepCMD + " -r " + fq1 + " -f " + fq2 + " -1 " + unmerged1 + " -2 " + unmerged2 + " -s "+merged
    if (!adaptorSequences.isEmpty) {
      cmd += " -A " + adaptorSequences(0) +  " -B " + adaptorSequences(1)
    }
    def commandLine = cmd

  }
  case class bwa_sam_se(inBAM: File, inSai: File, outBAM: File, readGroupString: String) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "bam file to be aligned") var bam = inBAM
    @Input(doc = "bwa alignment index file") var sai = inSai
    @Output(doc = "output aligned bam file") var alignedBam = outBAM

    def commandLine = bwaPath + " samse " + reference + " " + sai + " " + bam + " -r \"" + readGroupString + "\" > " + alignedBam

    this.memoryLimit = 6
    this.analysisName = outBAM + ".bwa_sam_se"
    this.jobName = outBAM + ".bwa_sam_se"
//    this.isIntermediate = false
  }

  case class bwa_sam_pe(inFile1: File, inFile2: File, inSai1: File, inSai2: File, outBAM: File, readGroupString: String) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "First file to be aligned") var first = inFile1
    @Input(doc = "Second file to be aligned") var second = inFile2
    @Input(doc = "bwa alignment index file for 1st mating pair") var sai1 = inSai1
    @Input(doc = "bwa alignment index file for 2nd mating pair") var sai2 = inSai2
    @Output(doc = "output aligned bam file") var alignedBam = outBAM

    def commandLine = bwaPath + " sampe " + reference + " " + sai1 + " " + sai2 + " " + first + " " + second + " -r \"" + readGroupString + "\" > " + alignedBam

    this.memoryLimit = 4
    this.analysisName = outBAM + ".bwa_sam_pe"
    this.jobName = outBAM + ".bwa_sam_pe"
//    this.isIntermediate = false
  }

  case class bwa_sam_pe_bam(inBam: File, inSai1: File, inSai2: File, outBAM: File, readGroupString: String) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Input(doc = "bwa alignment index file for 1st mating pair") var sai1 = inSai1
    @Input(doc = "bwa alignment index file for 2nd mating pair") var sai2 = inSai2
    @Output(doc = "output aligned bam file") var alignedBam = outBAM

    def commandLine = bwaPath + " sampe " + reference + " " + sai1 + " " + sai2 + " " + bam + " " + bam + " -r \"" + readGroupString + "\" > " + alignedBam
    this.memoryLimit = 4
    this.analysisName = outBAM + ".bwa_sam_pe"
    this.jobName = outBAM + ".bwa_sam_pe"
//    this.isIntermediate = false
  }

  case class bwa_sw(inFastQ: File, outBAM: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "fastq file to be aligned") var fq = inFastQ
    @Output(doc = "output bam file") var bam = outBAM

    def commandLine = bwaPath + " bwasw -t " + numThreads + " " + reference + " " + fq + " > " + bam

    this.analysisName = outBAM + ".bwasw"
    this.jobName = outBAM + ".bwasw"
  }


  case class bwa_aln(inBAM: File, outSai: File) extends CommandLineFunction {
    @Input(doc = "bam file to be aligned") var bam = inBAM
    @Output(doc = "output sai file") var sai = outSai

    var ex:String = ""
    if (inBAM.endsWith("bam") ) ex = " -b "

    def commandLine = bwaPath + " aln -t " + bwaThreads + bwaParameters + reference + " " + bam + ex + " > " + sai
    this.analysisName = outSai + ".bwa_aln_se"
    this.jobName = outSai + ".bwa_aln_se"
    this.memoryLimit = Some(6)
    this.isIntermediate = false
    this.jobQueue = longQueue

  }
  case class bwa_aln_pe (inBam: File, outSai1: File, index: Int) extends CommandLineFunction {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Output(doc="output sai file for 1st mating pair") var sai = outSai1
    def commandLine = bwaPath + " aln -t " + bwaThreads + bwaParameters + reference + " -b" + index + " " + bam + " > " + sai
    this.analysisName = outSai1 + ".bwa_aln_pe1"
    this.jobName = outSai1 + ".bwa_aln_pe1"
    this.memoryLimit = Some(8)
    this.isIntermediate = false
    this.jobQueue = longQueue
  }

  case class writeList(inBAMs: Seq[File], outBAMList: File) extends ListWriterFunction {
    this.inputFiles = inBAMs
    this.listFile = outBAMList
    this.analysisName = outBAMList + ".bamList"
    this.jobName = outBAMList + ".bamList"
  }


}

