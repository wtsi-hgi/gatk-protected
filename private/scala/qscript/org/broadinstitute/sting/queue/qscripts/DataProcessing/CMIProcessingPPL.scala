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

import collection.JavaConversions._
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileHeader.SortOrder

import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException

class CMIProcessingPPL extends QScript {
  qscript =>

  /****************************************************************************
    * Required Parameters
    ****************************************************************************/

  @Input(doc="Single ended or first of pair ended fastQ list", fullName="normal_fastq_1", shortName="nf1", required=false, exclusiveOf = "normalBAMs")
  var normalFastQ1: File = _

  @Input(doc="Single or pair ended BAM list", fullName="bam_list", shortName="bl", required=false, exclusiveOf = "normalFastQ1")
  var normalBAMs: File = _

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=true)
  var reference: File = _

  @Input(doc="DBSNP or known callset to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=true)
  var dbSNP: Seq[File] = Seq()

  /****************************************************************************
    * Optional Input Parameters
    ****************************************************************************/

  @Input(doc="Second of pair ended fastQ", fullName="normal_fastq_2", shortName="nf2", required=false)
  var normalFastQ2: File = _

  @Input(doc="Single ended or first of pair ended TUMOR fastQ list", fullName="tumor_fastq_1", shortName="tf1", required=false)
  var tumorFastQ1: File = _

  @Input(doc="Second of pair ended TUMOR fastQ list", fullName="tumor_fastq_2", shortName="tf2", required=false)
  var tumorFastQ2: File = _

  @Input(doc="Single or pair ended BAM list", fullName="bam", shortName="bam", required=true)
  var tumorBAMs: File = _

  @Argument(doc="Perform pair ended analysis", fullName="pair_ended", shortName="pe", required=false)
  var pairEndedAnalysis: Boolean = false


  /****************************************************************************
    * Additional Parameters
    ****************************************************************************/

  @Input(doc="extra VCF files to use as reference indels for Indel Realignment", fullName="extra_indels", shortName="indels", required=false)
  var indels: Seq[File] = Seq()

  @Input(doc="The path to the binary of bwa (usually BAM files have already been mapped - but if you want to remap this is the option)", fullName="path_to_bwa", shortName="bwa", required=false)
  var bwaPath: File = _

  @Argument(doc="Use BWASW instead of BWA aln", fullName="use_bwa_sw", shortName="bwasw", required=false)
  var useBWAsw: Boolean = false

  @Argument(doc="Number of threads BWA should use", fullName="bwa_threads", shortName="bt", required=false)
  var bwaThreads: Int = 1

  @Argument(doc="Perform validation on the BAM files", fullName="validation", shortName="vs", required=false)
  var validation: Boolean = false


  /****************************************************************************
    * Hidden Parameters
    ****************************************************************************/
  @Hidden
  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = 0

  @Hidden
  @Argument(doc="Define the default platform for Count Covariates -- useful for techdev purposes only.", fullName="default_platform", shortName="dp", required=false)
  var defaultPlatform: String = ""

  @Hidden
  @Argument(doc="Run the pipeline in test mode only", fullName = "test_mode", shortName = "test", required=false)
  var testMode: Boolean = false


  /****************************************************************************
    * Global Variables
    ****************************************************************************/
  val bwaParameters: String = " -q 5 -l 32 -k 2 -t 4 -o 1 "


  /****************************************************************************
    * Main script
    ****************************************************************************/

  def script() {

    if (normalFastQ1 == null && normalBAMs == null) {
      throw new ReviewedStingException("You need to provide either a BAM list (-bl) or a fastQ list (-nf1) for this analysis")
    }

    val useBAMs = normalBAMs != null
    val normals1: Seq[File] = if (useBAMs) {QScriptUtils.createSeqFromFile(normalBAMs)} else {QScriptUtils.createSeqFromFile(normalFastQ1)}
    val normals2: Seq[File] = if (pairEndedAnalysis && !useBAMs) {QScriptUtils.createSeqFromFile(normalFastQ2)} else {Seq()}
    val tumors1: Seq[File] = if (useBAMs) {QScriptUtils.createSeqFromFile(tumorBAMs)} else {QScriptUtils.createSeqFromFile(tumorFastQ1)}
    val tumors2: Seq[File] = if (pairEndedAnalysis && !useBAMs) {QScriptUtils.createSeqFromFile(tumorFastQ2)} else {Seq()}

    if (pairEndedAnalysis && !useBAMs) {
      assert(normals1.length == normals2.length, "Different number of first and second of pairs in the lists provided " + normals1.length + " != " + normals2.length)
    }

    val alignedNormals = performAlignment(normals1, normals2, useBAMs)
    val alignedTumors = performAlignment(tumors1, tumors2, useBAMs)
    val allBAMs = alignedNormals ++ alignedTumors


    // BAM files generated by the pipeline
    val normal       = alignedNormals(0)
    val normalClean  = swapExt(normal, ".bam", ".clean.bam")
    val normalDedup  = swapExt(normal, ".bam", ".clean.dedup.bam")
    val normalRecal  = swapExt(normal, ".bam", ".clean.dedup.recal.bam")

    // Accessory files
    val targetIntervals = swapExt(normal, ".bam", ".intervals")
    val metricsFile     = swapExt(normal, ".bam", ".metrics")
    val preRecalFile    = swapExt(normal, ".bam", ".pre_recal.table")
    val postRecalFile   = swapExt(normal, ".bam", ".post_recal.table")


    add(target(allBAMs, targetIntervals),
        clean(allBAMs, targetIntervals, normalClean),
        dedup(normalClean, normalDedup, metricsFile),
        cov(normalDedup, preRecalFile),
        recal(normalDedup, preRecalFile, normalRecal),
        cov(normalRecal, postRecalFile))

    if (!alignedTumors.isEmpty) {
      val tumor        = alignedTumors(0)
      val tumorClean   = swapExt(tumor, ".bam", ".clean.bam")
      val tumorDedup   = swapExt(tumor, ".bam", ".clean.dedup.bam")
      val tumorRecal   = swapExt(tumor, ".bam", ".clean.dedup.recal.bam")

      val tumorTargetIntervals = swapExt(normal, ".bam", ".intervals")
      val tumorMetricsFile     = swapExt(normal, ".bam", ".metrics")
      val tumorPreRecalFile    = swapExt(normal, ".bam", ".pre_recal.table")
      val tumorPostRecalFile   = swapExt(normal, ".bam", ".post_recal.table")


    add(dedup(tumorClean, tumorDedup, tumorMetricsFile),
        cov(tumorDedup, tumorPreRecalFile),
        recal(tumorDedup, tumorPreRecalFile, tumorRecal),
        cov(tumorRecal, tumorPostRecalFile))
    }
  }






  /****************************************************************************
    * Helper classes and methods
    ****************************************************************************/

  class ReadGroup (val id: String,
                   val lb: String,
                   val pl: String,
                   val pu: String,
                   val sm: String,
                   val cn: String,
                   val ds: String)
  {}

  // Utility function to merge all bam files of similar samples. Generates one BAM file per sample.
  // It uses the sample information on the header of the input BAM files.
  //
  // Because the realignment only happens after these scripts are executed, in case you are using
  // bwa realignment, this function will operate over the original bam files and output over the
  // (to be realigned) bam files.
  def createSampleFiles(bamFiles: Seq[File], realignedBamFiles: Seq[File]): Map[String, Seq[File]] = {

    // Creating a table with SAMPLE information from each input BAM file
    val sampleTable = scala.collection.mutable.Map.empty[String, Seq[File]]
    val realignedIterator = realignedBamFiles.iterator
    for (bam <- bamFiles) {
      val rBam = realignedIterator.next()  // advance to next element in the realignedBam list so they're in sync.

      val samReader = new SAMFileReader(bam)
      val header = samReader.getFileHeader
      val readGroups = header.getReadGroups

      // only allow one sample per file. Bam files with multiple samples would require pre-processing of the file
      // with PrintReads to separate the samples. Tell user to do it himself!
      assert(!QScriptUtils.hasMultipleSamples(readGroups), "The pipeline requires that only one sample is present in a BAM file. Please separate the samples in " + bam)

      // Fill out the sample table with the readgroups in this file
      for (rg <- readGroups) {
        val sample = rg.getSample
        if (!sampleTable.contains(sample))
          sampleTable(sample) = Seq(rBam)
        else if ( !sampleTable(sample).contains(rBam))
          sampleTable(sample) :+= rBam
      }
    }
    sampleTable.toMap
  }

  // Rebuilds the Read Group string to give BWA
  def addReadGroups(inBam: File, outBam: File, samReader: SAMFileReader) {
    val readGroups = samReader.getFileHeader.getReadGroups
    var index: Int = readGroups.length
    for (rg <- readGroups) {
      val intermediateInBam: File = if (index == readGroups.length) { inBam } else { swapExt(outBam, ".bam", index+1 + "-rg.bam") }
      val intermediateOutBam: File = if (index > 1) {swapExt(outBam, ".bam", index + "-rg.bam") } else { outBam}
      val readGroup = new ReadGroup(rg.getReadGroupId, rg.getLibrary, rg.getPlatform, rg.getPlatformUnit, rg.getSample, rg.getSequencingCenter, rg.getDescription)
      add(addReadGroup(intermediateInBam, intermediateOutBam, readGroup))
      index = index - 1
    }
  }


  /**
   *
   */
  def generateBWAInputParameters(isBAM: Boolean, isSecondOfPair: Boolean):String = {
    var parms: String = ""
    if (isBAM) {
      parms = if (pairEndedAnalysis) {if (isSecondOfPair) {" -b2 "} else {" -b1 "}} else {" -b "}
    }
    parms
  }

  /**
   * BWA alignment
   *
   * @return
   */
  def performAlignment(firstOfPairs: Seq[File], secondOfPairs: Seq[File], useBAMs: Boolean): Seq[File] = {
    var result: Seq[File] = Seq()

    for (i <- 0 until firstOfPairs.length) {
      val first = firstOfPairs(i)
      val splitted = first.toString.split(".")
      val extension = "." + splitted(splitted.length - 1)
      val saiFile1: File = swapExt(first, extension, ".1.sai")
      val saiFile2: File = swapExt(first, extension, ".2.sai")
      val alignedSAM: File = swapExt(first, extension, ".sam")
      val alignedBAM: File = swapExt(alignedSAM, ".sam", ".aln.bam")
      val rgBAM : File = swapExt(alignedBAM, ".bam", ".rg.bam")

      // BAM pipeline
      if (useBAMs) {
        val revertedBAM = revertBAM(first, true)
        add(bwa(generateBWAInputParameters(useBAMs, false), revertedBAM, saiFile1))
        if (pairEndedAnalysis) {
          add(bwa(generateBWAInputParameters(true, true), revertedBAM, saiFile2),
              bwa_sam_pe(revertedBAM, revertedBAM, saiFile1, saiFile2, alignedSAM))
        }
        else {
          add(bwa_sam_se(revertedBAM, saiFile1, alignedSAM))
        }
      }

      // FASTQ pipeline
      else {
        add(bwa(generateBWAInputParameters(false, false), first, saiFile1))
        if (pairEndedAnalysis) {
          val second = secondOfPairs(i)
          add(bwa(generateBWAInputParameters(false, true), second, saiFile2),
              bwa_sam_pe(first, second, saiFile1, saiFile2, alignedSAM))
        }
        else {
          add(bwa_sam_se(first, saiFile1, alignedSAM))
        }
      }

      add(sortSam(alignedSAM, alignedBAM, SortOrder.coordinate))
      addReadGroups(alignedBAM, rgBAM, new SAMFileReader(alignedSAM))
      result :+= rgBAM
    }

    result
  }

  def revertBams(bams: Seq[File], removeAlignmentInformation: Boolean): Seq[File] = {
    var revertedBAMList: Seq[File] = Seq()
    for (bam <- bams)
      revertedBAMList :+= revertBAM(bam, removeAlignmentInformation)
    revertedBAMList
  }

  def revertBAM(bam: File, removeAlignmentInformation: Boolean): File = {
    val revertedBAM = swapExt(bam, ".bam", ".reverted.bam")
    add(revert(bam, revertedBAM, removeAlignmentInformation))
    revertedBAM
  }




  /****************************************************************************
    * Classes (GATK Walkers)
    ****************************************************************************/



  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 4
    this.isIntermediate = true
  }

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
  }

  trait SAMargs extends PicardBamFunction with ExternalCommonArgs {
    this.maxRecordsInRam = 100000
  }

  case class target (inBams: Seq[File], outIntervals: File) extends RealignerTargetCreator with CommandLineGATKArgs {
    this.input_file = inBams
    this.out = outIntervals
    this.mismatchFraction = 0.0
    this.known ++= qscript.dbSNP
    if (indels != null)
      this.known ++= qscript.indels
    this.scatterCount = nContigs
    this.analysisName = outIntervals + ".target"
    this.jobName = outIntervals + ".target"
  }

  case class clean (inBams: Seq[File], tIntervals: File, outBam: File) extends IndelRealigner with CommandLineGATKArgs {
    this.input_file = inBams
    this.targetIntervals = tIntervals
    this.out = outBam
    this.known ++= qscript.dbSNP
    if (qscript.indels != null)
      this.known ++= qscript.indels
    this.consensusDeterminationModel = ConsensusDeterminationModel.USE_READS
    this.compress = 0
    this.noPGTag = qscript.testMode
    this.scatterCount = nContigs
    this.analysisName = outBam + ".clean"
    this.jobName = outBam + ".clean"
  }

  case class cov (inBam: File, outRecalFile: File) extends BaseRecalibrator with CommandLineGATKArgs {
    this.knownSites ++= qscript.dbSNP
    this.covariate ++= Seq("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "ContextCovariate")
    this.input_file :+= inBam
    this.disable_indel_quals = true
    this.out = outRecalFile
    if (!defaultPlatform.isEmpty) this.default_platform = defaultPlatform
    this.scatterCount = nContigs
    this.analysisName = outRecalFile + ".covariates"
    this.jobName = outRecalFile + ".covariates"
  }

  case class recal (inBam: File, inRecalFile: File, outBam: File) extends PrintReads with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.BQSR = inRecalFile
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.out = outBam
    this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = outBam + ".recalibration"
    this.jobName = outBam + ".recalibration"
  }



  /****************************************************************************
    * Classes (non-GATK programs)
    ****************************************************************************/


  case class dedup (inBam: File, outBam: File, metricsFile: File) extends MarkDuplicates with ExternalCommonArgs {
    this.input :+= inBam
    this.output = outBam
    this.metrics = metricsFile
    this.memoryLimit = 16
    this.analysisName = outBam + ".dedup"
    this.jobName = outBam + ".dedup"
  }

  case class joinBams (inBams: Seq[File], outBam: File) extends MergeSamFiles with ExternalCommonArgs {
    this.input = inBams
    this.output = outBam
    this.analysisName = outBam + ".joinBams"
    this.jobName = outBam + ".joinBams"
  }

  case class sortSam (inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam with ExternalCommonArgs {
    this.input :+= inSam
    this.output = outBam
    this.sortOrder = sortOrderP
    this.analysisName = outBam + ".sortSam"
    this.jobName = outBam + ".sortSam"
  }

  case class validate (inBam: File, outLog: File) extends ValidateSamFile with ExternalCommonArgs {
    this.input :+= inBam
    this.output = outLog
    this.REFERENCE_SEQUENCE = qscript.reference
    this.isIntermediate = false
    this.analysisName = outLog + ".validate"
    this.jobName = outLog + ".validate"
  }


  case class addReadGroup (inBam: File, outBam: File, readGroup: ReadGroup) extends AddOrReplaceReadGroups with ExternalCommonArgs {
    this.input :+= inBam
    this.output = outBam
    this.RGID = readGroup.id
    this.RGCN = readGroup.cn
    this.RGDS = readGroup.ds
    this.RGLB = readGroup.lb
    this.RGPL = readGroup.pl
    this.RGPU = readGroup.pu
    this.RGSM = readGroup.sm
    this.analysisName = outBam + ".rg"
    this.jobName = outBam + ".rg"
  }

  case class revert (inBam: File, outBam: File, removeAlignmentInfo: Boolean) extends RevertSam with ExternalCommonArgs {
    this.output = outBam
    this.input :+= inBam
    this.removeAlignmentInformation = removeAlignmentInfo
    this.sortOrder = if (removeAlignmentInfo) {SortOrder.queryname} else {SortOrder.coordinate}
    this.analysisName = outBam + "revert"
    this.jobName = outBam + ".revert"
  }

  case class convertToFastQ (inBam: File, outFQ: File) extends SamToFastq with ExternalCommonArgs {
    this.input :+= inBam
    this.fastq = outFQ
    this.analysisName = outFQ + "convert_to_fastq"
    this.jobName = outFQ + ".convert_to_fastq"
  }

  case class bwa_aln_se (inBam: File, outSai: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Output(doc="output sai file") var sai = outSai
    def commandLine = bwaPath + " aln -t " + bwaThreads + " -q 5 " + reference + " -b " + bam + " > " + sai
    this.analysisName = outSai + ".bwa_aln_se"
    this.jobName = outSai + ".bwa_aln_se"
  }

  case class bwa_aln_pe (inBam: File, outSai1: File, index: Int) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Output(doc="output sai file for 1st mating pair") var sai = outSai1
    def commandLine = bwaPath + " aln -t " + bwaThreads + " -q 5 " + reference + " -b" + index + " " + bam + " > " + sai
    this.analysisName = outSai1 + ".bwa_aln_pe1"
    this.jobName = outSai1 + ".bwa_aln_pe1"
  }

  case class bwa_sam_se (inBam: File, inSai: File, outBam: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Input(doc="bwa alignment index file") var sai = inSai
    @Output(doc="output aligned bam file") var alignedBam = outBam
    def commandLine = bwaPath + " samse " + reference + " " + sai + " " + bam + " > " + alignedBam
    this.memoryLimit = 6
    this.analysisName = outBam + ".bwa_sam_se"
    this.jobName = outBam + ".bwa_sam_se"
  }

  case class bwa_sam_pe (inFile1: File, inFile2: File, inSai1: File, inSai2:File, outBam: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var first = inFile1
    @Input(doc="bam file to be aligned") var second = inFile2
    @Input(doc="bwa alignment index file for 1st mating pair") var sai1 = inSai1
    @Input(doc="bwa alignment index file for 2nd mating pair") var sai2 = inSai2
    @Output(doc="output aligned bam file") var alignedBam = outBam
    def commandLine = bwaPath + " sampe " + reference + " " + sai1 + " " + sai2 + " " + first + " " + second + " > " + alignedBam
    this.memoryLimit = 6
    this.analysisName = outBam + ".bwa_sam_pe"
    this.jobName = outBam + ".bwa_sam_pe"
  }

  case class bwa_sw (inFastQ: File, outBam: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="fastq file to be aligned") var fq = inFastQ
    @Output(doc="output bam file") var bam = outBam
    def commandLine = bwaPath + " bwasw -t " + bwaThreads + " " + reference + " " + fq + " > " + bam
    this.analysisName = outBam + ".bwasw"
    this.jobName = outBam + ".bwasw"
  }
  
  
  case class bwa (inputParms: String, inBam: File, outSai: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Output(doc="output sai file") var sai = outSai
    def commandLine = bwaPath + " aln -t " + bwaThreads + bwaParameters + reference + inputParms + bam + " > " + sai
    this.analysisName = outSai + ".bwa_aln_se"
    this.jobName = outSai + ".bwa_aln_se"
  }

  case class writeList(inBams: Seq[File], outBamList: File) extends ListWriterFunction {
    this.inputFiles = inBams
    this.listFile = outBamList
    this.analysisName = outBamList + ".bamList"
    this.jobName = outBamList + ".bamList"
  }
}

