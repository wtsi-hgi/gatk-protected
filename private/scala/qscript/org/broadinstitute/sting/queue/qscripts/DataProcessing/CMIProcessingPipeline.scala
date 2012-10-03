/**
 * Created with IntelliJ IDEA.
 * User: carneiro
 * Date: 9/25/12
 * Time: 12:04 PM
 */

package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.picard._
import org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode

import net.sf.samtools.SAMFileHeader.SortOrder

import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.commandline.Hidden
import io.Source
import org.broadinstitute.sting.utils.NGSPlatform
import collection.mutable

class CMIProcessingPipeline extends QScript {
  qscript =>

  /*****************************************************************************
    * Required Parameters
    ****************************************************************************/

  @Input(doc="a table with all the necessary information to process the data", fullName="metadata", shortName="m", required=true)
  var metaData: File = _

  /********************************************************************************
    * Additional Parameters that the pipeline should have pre-defined in the image
    *******************************************************************************/

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=true)
  var reference: File = _

  @Input(doc="DBSNP or known callset to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=true)
  var dbSNP: Seq[File] = Seq()

  @Input(doc="The path to the binary of bwa (usually BAM files have already been mapped - but if you want to remap this is the option)", fullName="path_to_bwa", shortName="bwa", required=false)
  var bwaPath: File = _

  @Input(doc="extra VCF files to use as reference indels for Indel Realignment", fullName="extra_indels", shortName="indels", required=false)
  var indels: Seq[File] = Seq()

  /****************************************************************************
    * Hidden Parameters
    ****************************************************************************/
  @Hidden
  @Argument(doc="Use BWASW instead of BWA aln", fullName="use_bwa_sw", shortName="bwasw", required=false)
  var useBWAsw: Boolean = false

  @Hidden
  @Argument(doc="Number of threads BWA should use", fullName="bwa_threads", shortName="bt", required=false)
  var bwaThreads: Int = 1

  @Hidden
  @Argument(doc="Number of threads BWA should use", fullName="mem_limit", shortName="mem", required=false)
  var memLimit: Int = 4

  @Hidden
  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = 0

  @Hidden
  @Argument(doc="Define the default platform for Count Covariates -- useful for techdev purposes only.", fullName="default_platform", shortName="dp", required=false)
  var defaultPlatform: String = ""

  @Hidden
  @Argument(doc="Run the pipeline in test mode only", fullName = "test_mode", shortName = "test", required=false)
  var testMode: Boolean = false

  @Hidden
  @Argument(doc="BWA Parameteres", fullName = "bwa_parameters", shortName = "bp", required=false)
  val bwaParameters: String = " -q 5 -l 32 -k 2 -t 4 -o 1 "

  val cleaningExtension: String = ".clean.bam"
  val headerVersion: String = "#FILE1,FILE2,INDIVIDUAL,SAMPLE,LIBRARY,SEQUENCING,TUMOR,PLATFORM,PLATFORM_UNIT,CENTER,DESCRIPTION,DATE_SEQUENCED"

  /****************************************************************************
    * Main script
    ****************************************************************************/

  def script() {

    // todo -- preprocess metadata
    var lanes = Seq[MetaInfo]()
    var id = 1
    for (line: String <- Source.fromFile(metaData).getLines()) {
      if (line.startsWith("#")) {
        checkMetaDataHeader(line)
      }
      else {
        lanes :+= new MetaInfo(id, line.split(","))
        id = id + 1
      }
    }

    // todo -- align fastQ's from lanes object
    val samples = new mutable.HashMap[String, Seq[File]]()
    for (meta <- lanes) {
      val bamFile = performAlignment(meta)
      if (samples.contains(meta.sample)) {
        samples.put(meta.sample, samples(meta.sample) ++ Seq(bamFile))
      }
      else {
        samples.put(meta.sample, Seq(bamFile))
      }
    }
    // todo -- add optional bam de-multiplexing and re-multiplexing pipeline

    var allBAMs = Seq[File]()
    for ((sample,bams) <- samples) {
      val sampleBAM = new File(sample + ".bam")
      allBAMs +:= sampleBAM
      add(joinBAMs(bams, sampleBAM))
    }

    clean(allBAMs)
    for (bam <- allBAMs) {
      val cleanBAM      = swapExt(bam, ".bam", cleaningExtension)
      val dedupBAM      = swapExt(bam, ".bam", ".clean.dedup.bam")
      val recalBAM      = swapExt(bam, ".bam", ".clean.dedup.recal.bam")
      val reducedBAM    = swapExt(bam, ".bam", ".clean.dedup.recal.reduced.bam")
      val metricsFile   = swapExt(bam, ".bam", ".metrics")
      val preRecalFile  = swapExt(bam, ".bam", ".pre_recal.table")
      val postRecalFile = swapExt(bam, ".bam", ".post_recal.table")

      add(dedup(cleanBAM, dedupBAM, metricsFile))
      recalibrate(dedupBAM, preRecalFile, postRecalFile, recalBAM)
      add(reduce(recalBAM, reducedBAM))
    }
  }

  /****************************************************************************
    * Helper classes and methods
    ****************************************************************************/

  private class MetaInfo (
    val id: Int,
    val file1: File,
    val file2: File,
    val individual: String,
    val sample: String,
    val library: String,
    val sequencing: String,
    val tumor: Int,
    val platform: NGSPlatform,
    val platformUnit: String,
    val center: String,
    val description: String,
    val dateSequenced: String
  ) {
    def this(idP: Int, headerArray: Array[String]) =
      this (
        idP,
        new File(headerArray(0)),
        if (headerArray(1).isEmpty) {null} else {new File(headerArray(1))},
        headerArray(2),
        headerArray(3),
        headerArray(4),
        headerArray(5),
        headerArray(6).toInt,
        NGSPlatform.fromReadGroupPL(headerArray(7)),
        headerArray(8),
        headerArray(9),
        headerArray(10),
        headerArray(11)
      )
    def bamFileName = individual + "." + sample + "." + library + "." + sequencing + "." + id + "." + tumor + ".bam"
    def readGroupString = "@RG\tID:%d\tCN:%s\tDS:%s\tDT:%s\tLB:%s\tPL:%s\tPU:%s\tSM:%s".format(id, center, description, dateSequenced, library, platform, platformUnit, sample)
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
    val saiFile1: File = new File(metaInfo.file1 + ".1.sai")
    val saiFile2: File = new File(metaInfo.file2 + ".2.sai")
    val alnSAM: File   = new File(metaInfo.file1 + ".sam")
    val alnBAM: File   = new File(metaInfo.bamFileName)

    add(bwa(" ", metaInfo.file1, saiFile1))

    if (!metaInfo.file2.isEmpty) {
      add(bwa(" ", metaInfo.file2, saiFile2),
          bwa_sam_pe(metaInfo.file1, metaInfo.file2, saiFile1, saiFile2, alnSAM, metaInfo.readGroupString))
    }
    else {
          add(bwa_sam_se(metaInfo.file1, saiFile1, alnSAM, metaInfo.readGroupString))
    }
    add(sortSam(alnSAM, alnBAM, SortOrder.coordinate))

    alnBAM
  }

  def clean(allBAMs: Seq[File]) {
    val bam: File = allBAMs(0)
    val targetIntervals = swapExt(bam, ".bam", ".cleaning_intervals")
    add(target(allBAMs, targetIntervals), indel(allBAMs, targetIntervals))
  }

  def recalibrate(dedupBAM: File, preRecalFile: File, postRecalFile: File, recalBAM: File) {
    add(bqsr(dedupBAM, preRecalFile),
        apply_bqsr(dedupBAM, preRecalFile, recalBAM),
        bqsr(recalBAM, postRecalFile)) 
  }


  /****************************************************************************
    * Classes (GATK Walkers)
    ****************************************************************************/



  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = qscript.memLimit
    this.isIntermediate = true
  }

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
  }

  trait SAMargs extends PicardBamFunction with ExternalCommonArgs {
    this.maxRecordsInRam = 100000
  }

  case class target (inBAMs: Seq[File], outIntervals: File) extends RealignerTargetCreator with CommandLineGATKArgs {
    this.input_file = inBAMs
    this.out = outIntervals
    this.mismatchFraction = 0.0
    this.known ++= qscript.dbSNP
    if (indels != null)
      this.known ++= qscript.indels
    this.scatterCount = nContigs
    this.analysisName = outIntervals + ".target"
    this.jobName = outIntervals + ".target"
  }

  case class indel (inBAMs: Seq[File], tIntervals: File) extends IndelRealigner with CommandLineGATKArgs {
    // TODO -- THIS IS A WORKAROUND FOR QUEUE'S LIMITATION OF TRACKING LISTS OF FILES (implementation limited to 5 files)
    @Output(doc="first cleaned bam file", required=true) var out1: File = swapExt(inBAMs(0), ".bam", cleaningExtension)
    @Output(doc="first cleaned bam file", required=true) var ind1: File = swapExt(out1, ".bam", ".bai")
    @Output(doc="first cleaned bam file", required=false) var out2: File = if (inBAMs.length >= 2) {swapExt(inBAMs(1), ".bam", cleaningExtension)} else {null}
    @Output(doc="first cleaned bam file", required=false) var ind2: File = if (inBAMs.length >= 2) {swapExt(out2, ".bam", ".bai")} else {null}
    @Output(doc="first cleaned bam file", required=false) var out3: File = if (inBAMs.length >= 3) {swapExt(inBAMs(2), ".bam", cleaningExtension)} else {null}
    @Output(doc="first cleaned bam file", required=false) var ind3: File = if (inBAMs.length >= 2) {swapExt(out3, ".bam", ".bai")} else {null}
    @Output(doc="first cleaned bam file", required=false) var out4: File = if (inBAMs.length >= 4) {swapExt(inBAMs(3), ".bam", cleaningExtension)} else {null}
    @Output(doc="first cleaned bam file", required=false) var ind4: File = if (inBAMs.length >= 2) {swapExt(out4, ".bam", ".bai")} else {null}
    @Output(doc="first cleaned bam file", required=false) var out5: File = if (inBAMs.length >= 5) {swapExt(inBAMs(4), ".bam", cleaningExtension)} else {null}
    @Output(doc="first cleaned bam file", required=false) var ind5: File = if (inBAMs.length >= 2) {swapExt(out5, ".bam", ".bai")} else {null}
    this.input_file = inBAMs
    this.targetIntervals = tIntervals
    this.nWayOut = cleaningExtension
    this.known ++= qscript.dbSNP
    if (qscript.indels != null)
      this.known ++= qscript.indels
    this.consensusDeterminationModel = ConsensusDeterminationModel.USE_READS
    this.compress = 0
    this.noPGTag = qscript.testMode
    this.scatterCount = nContigs
    this.analysisName = inBAMs(0).toString + "clean"
    this.jobName = inBAMs(0).toString + ".clean"
  }

  case class bqsr (inBAM: File, outRecalFile: File) extends BaseRecalibrator with CommandLineGATKArgs {
    this.knownSites ++= qscript.dbSNP
    this.covariate ++= Seq("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "ContextCovariate")
    this.input_file :+= inBAM
    this.disable_indel_quals = true
    this.out = outRecalFile
    if (!defaultPlatform.isEmpty) this.default_platform = defaultPlatform
    this.scatterCount = nContigs
    this.analysisName = outRecalFile + ".covariates"
    this.jobName = outRecalFile + ".covariates"
  }

  case class apply_bqsr (inBAM: File, inRecalFile: File, outBAM: File) extends PrintReads with CommandLineGATKArgs {
    this.input_file :+= inBAM
    this.BQSR = inRecalFile
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.out = outBAM
    this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = outBAM + ".recalibration"
    this.jobName = outBAM + ".recalibration"
  }
  
  case class reduce (inBAM: File, outBAM: File) extends ReduceReads with CommandLineGATKArgs {
    this.input_file :+= inBAM
    this.out = outBAM
    this.isIntermediate = false
    this.analysisName = outBAM + ".recalibration"
    this.jobName = outBAM + ".recalibration"
  }



  /****************************************************************************
    * Classes (non-GATK programs)
    ****************************************************************************/


  case class dedup (inBAM: File, outBAM: File, metricsFile: File) extends MarkDuplicates with ExternalCommonArgs {
    this.input :+= inBAM
    this.output = outBAM
    this.metrics = metricsFile
    this.memoryLimit = 16
    this.analysisName = outBAM + ".dedup"
    this.jobName = outBAM + ".dedup"
  }

  case class joinBAMs (inBAMs: Seq[File], outBAM: File) extends MergeSamFiles with ExternalCommonArgs {
    this.input = inBAMs
    this.output = outBAM
    this.analysisName = outBAM + ".joinBAMs"
    this.jobName = outBAM + ".joinBAMs"
  }

  case class sortSam (inSam: File, outBAM: File, sortOrderP: SortOrder) extends SortSam with ExternalCommonArgs {
    this.input :+= inSam
    this.output = outBAM
    this.sortOrder = sortOrderP
    this.analysisName = outBAM + ".sortSam"
    this.jobName = outBAM + ".sortSam"
  }

  case class validate (inBAM: File, outLog: File) extends ValidateSamFile with ExternalCommonArgs {
    this.input :+= inBAM
    this.output = outLog
    this.REFERENCE_SEQUENCE = qscript.reference
    this.isIntermediate = false
    this.analysisName = outLog + ".validate"
    this.jobName = outLog + ".validate"
  }

  case class revert (inBAM: File, outBAM: File, removeAlignmentInfo: Boolean) extends RevertSam with ExternalCommonArgs {
    this.output = outBAM
    this.input :+= inBAM
    this.removeAlignmentInformation = removeAlignmentInfo
    this.sortOrder = if (removeAlignmentInfo) {SortOrder.queryname} else {SortOrder.coordinate}
    this.analysisName = outBAM + "revert"
    this.jobName = outBAM + ".revert"
  }

  case class convertToFastQ (inBAM: File, outFQ: File) extends SamToFastq with ExternalCommonArgs {
    this.input :+= inBAM
    this.fastq = outFQ
    this.analysisName = outFQ + "convert_to_fastq"
    this.jobName = outFQ + ".convert_to_fastq"
  }

  case class bwa_sam_se (inBAM: File, inSai: File, outBAM: File, readGroupString: String) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBAM
    @Input(doc="bwa alignment index file") var sai = inSai
    @Output(doc="output aligned bam file") var alignedBam = outBAM
    def commandLine = bwaPath + " samse " + reference + " " + sai + " " + bam + " -r \"" + readGroupString + "\" > " + alignedBam
    this.memoryLimit = 6
    this.analysisName = outBAM + ".bwa_sam_se"
    this.jobName = outBAM + ".bwa_sam_se"
  }

  case class bwa_sam_pe (inFile1: File, inFile2: File, inSai1: File, inSai2:File, outBAM: File, readGroupString: String) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="First file to be aligned") var first = inFile1
    @Input(doc="Second file to be aligned") var second = inFile2
    @Input(doc="bwa alignment index file for 1st mating pair") var sai1 = inSai1
    @Input(doc="bwa alignment index file for 2nd mating pair") var sai2 = inSai2
    @Output(doc="output aligned bam file") var alignedBam = outBAM
    def commandLine = bwaPath + " sampe " + reference + " " + sai1 + " " + sai2 + " " + first + " " + second + " -r \"" + readGroupString + "\" > " + alignedBam
    this.memoryLimit = 6
    this.analysisName = outBAM + ".bwa_sam_pe"
    this.jobName = outBAM + ".bwa_sam_pe"
  }

  case class bwa_sw (inFastQ: File, outBAM: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="fastq file to be aligned") var fq = inFastQ
    @Output(doc="output bam file") var bam = outBAM
    def commandLine = bwaPath + " bwasw -t " + bwaThreads + " " + reference + " " + fq + " > " + bam
    this.analysisName = outBAM + ".bwasw"
    this.jobName = outBAM + ".bwasw"
  }
  
  
  case class bwa (inputParms: String, inBAM: File, outSai: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBAM
    @Output(doc="output sai file") var sai = outSai
    def commandLine = bwaPath + " aln -t " + bwaThreads + bwaParameters + reference + inputParms + bam + " > " + sai
    this.analysisName = outSai + ".bwa_aln_se"
    this.jobName = outSai + ".bwa_aln_se"
  }

  case class writeList(inBAMs: Seq[File], outBAMList: File) extends ListWriterFunction {
    this.inputFiles = inBAMs
    this.listFile = outBAMList
    this.analysisName = outBAMList + ".bamList"
    this.jobName = outBAMList + ".bamList"
  }
}

