/**
 * Created with IntelliJ IDEA.
 * User: carneiro
 * Date: 9/25/12
 * Time: 12:04 PM
 */

package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode

import net.sf.samtools.SAMFileHeader.SortOrder

import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.commandline.Hidden
import io.Source
import org.broadinstitute.sting.utils.NGSPlatform
import collection.mutable
import org.broadinstitute.sting.queue.extensions.picard._
import org.broadinstitute.sting.queue.extensions.cancer.MuTect

class CMIBAMProcessingPipeline extends QScript {
  qscript =>

  /*****************************************************************************
    * Required Parameters
    ****************************************************************************/

  @Input(doc="a table with all the necessary information to process the data", fullName="metadata", shortName="m", required=true)
  var metaData: File = _

  /********************************************************************************
    * Additional Parameters that the pipeline should have pre-defined in the image
    *******************************************************************************/

  // TODO: The metadata file should be decomposed into separate lists. Think SelectVariant's --filter/--filterName pairings.
  @Input(doc="argument that allows Queue to see and download files", fullName="file1", required=false)
  var file1: Seq[File] = Nil
  @Input(doc="even more files that should be downloaded", fullName="file2", required=false)
  var file2: Seq[File] = Nil

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=false)
  var reference: File = new File("/refdata/human_g1k_v37_decoy.fasta")

  @Input(doc="DBSNP or known callset to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=false)
  var dbSNP: Seq[File] = Seq(new File("/refdata/dbsnp_135.b37.vcf"))

  @Input(doc="The path to the binary of bwa (usually BAM files have already been mapped - but if you want to remap this is the option)", fullName="path_to_bwa", shortName="bwa", required=false)
  var bwaPath: File = "/opt/bwa/bin/bwa"

  @Input(doc="extra VCF files to use as reference indels for Indel Realignment", fullName="extra_indels", shortName="indels", required=false)
  var indelSites: Seq[File] = Seq()

  @Input(doc="Interval file with targets used in exome capture (used for QC metrics)", fullName="targets", shortName="targets", required=false)
  var targets: File = new File("/refdata/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list")

  @Input(doc="Interval file with baits used in exome capture (used for QC metrics)", fullName="baits", shortName="baits", required=false)
  var baits: File = new File("/refdata/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.interval_list")


  // todo hotfix: paste mutect parameters here
  @Argument(doc="Run mutect calling in resulting bam", fullName = "doMutect", shortName = "doMutect", required=false)
  var doMutect: Boolean = false

  @Input(doc="COSMIC sites to use (must be in VCF format)", fullName="cosmic", shortName="C", required=false)
  var cosmic: Seq[File] = Seq(new File("/refdata/hg19_cosmic_v54_120711.vcf"))

  @Input(doc="The path to the binary of MuTect", fullName="mutect_jar", shortName="mj", required=false)
  var mutectJar: File = new File("/opt/cmi-cancer/mutect/muTect.jar")

  @Input(doc="Panel Of Normals or known artifact sites to use (must be in VCF format)", fullName="panel_of_normals", shortName="pon", required=false)
  var pon: Seq[File] = Seq()

  /****************************************************************************
   * Output files, to be passed in by messaging service
   ****************************************************************************/

  @Output(doc="Processed unreduced normal BAM", fullName="unreducedNormalBAM", shortName="unb", required=true)  // Using full name, so json field is mixed case "unfilteredVcf" or "uv"
  var unreducedNormalBAM: File = _

  @Output(doc="Processed reduced normal BAMs", fullName="reducedNormalBAM", shortName="rnb", required=true)
  var reducedNormalBAM: File = _

  @Output(doc="Processed unreduced normal BAM Index", fullName="unreducedNormalBAMIndex", shortName="unbi", required=true)  // Using full name, so json field is mixed case "unfilteredVcf" or "uv"
  var unreducedNormalBAMIndex: File = _

  @Output(doc="Processed reduced normal BAM Index", fullName="reducedNormalBAMIndex", shortName="rnbi", required=true)
  var reducedNormalBAMIndex: File = _

  @Output(doc="Processed unreduced tumor BAM", fullName="unreducedTumorBAM", shortName="utb", required=false)  // Using full name, so json field is mixed case "unfilteredVcf" or "uv"
  var unreducedTumorBAM: File = _

  @Output(doc="Processed reduced tumor BAM", fullName="reducedTumorBAM", shortName="rtb", required=false)
  var reducedTumorBAM: File = _

  @Output(doc="Processed unreduced tumor BAM Index", fullName="unreducedTumorBAMIndex", shortName="utbi", required=false)  // Using full name, so json field is mixed case "unfilteredVcf" or "uv"
  var unreducedTumorBAMIndex: File = _

  @Output(doc="Processed reduced tumor BAM Index", fullName="reducedTumorBAMIndex", shortName="rtbi", required=false)
  var reducedTumorBAMIndex: File = _


  // in case single sample calls are requested
  @Output(doc="Processed single sample VCF", fullName="singleSampleVCF", shortName="ssvcf", required=false)  // Using full name, so json field is mixed case "unfilteredVcf" or "uv"
  var singleSampleVCF: File = _

  @Output(doc="Processed single sample VCF index", fullName="singleSampleVCFIndex", shortName="ssvcfi", required=false)
  var singleSampleVCFIndex: File = _

  // in case single sample calls are requested
/*
  @Output(doc="Processed single sample tumor point mutations", fullName="tumorPointMutationVCF", shortName="tvcf", required=false)  // Using full name, so json field is mixed case "unfilteredVcf" or "uv"
  var tumorPointMutationVCF: File = _

  @Output(doc="Processed single sample tumor point mutation index", fullName="tumorPointMutationVCFIndex", shortName="tvcfi", required=false)
  var tumorPointMutationVCFIndex: File = _

  @Output(doc="Tumor raw coverage file", fullName="tumorRawCoverage", shortName="trc", required=false)  // Using full name, so json field is mixed case "unfilteredVcf" or "uv"
  var tumorRawCoverage: File = _
   */

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
  @Argument(doc="Default memory limit per job", fullName="mem_limit", shortName="mem", required=false)
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
  @Argument(doc="Run single sample germline calling in resulting bam", fullName = "doSingleSampleCalling", shortName = "call", required=false)
  var doSingleSampleCalling: Boolean = false

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

    // todo -- parse all the sample level metadata into a single structure (not per lane)
    val tumorInfo = new mutable.HashMap[String, Int]()

    for (meta <- lanes) {
      val bamFile = performAlignment(meta)
      if (samples.contains(meta.sample)) {
        samples.put(meta.sample, samples(meta.sample) ++ Seq(bamFile))
      }
      else {
        samples.put(meta.sample, Seq(bamFile))
      }

      if (tumorInfo.contains(meta.sample)) {
        // check that the type of this sample hasn't changed
        assert(meta.tumor == tumorInfo(meta.sample),
          String.format("Tumor type for sample %s is internally inconsistent within metadata file.  Found %s and %s\n", meta.sample, meta.tumor.toString, tumorInfo(meta.sample).toString))
      } else {
        tumorInfo.put(meta.sample, meta.tumor)
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
      val duplicateMetricsFile   = swapExt(bam, ".bam", ".duplicateMetrics")
      val preRecalFile  = swapExt(bam, ".bam", ".pre_recal.table")
      val postRecalFile = swapExt(bam, ".bam", ".post_recal.table")
      val outVCF        = swapExt(reducedBAM,".bam",".vcf")

      add(dedup(cleanBAM, dedupBAM, duplicateMetricsFile))
      if (qscript.targets != null && qscript.baits != null) {
        add(calculateHSMetrics(dedupBAM,swapExt(dedupBAM,".bam",".hs_metrics")))
      }
      recalibrate(dedupBAM, preRecalFile, postRecalFile, recalBAM)
      add(reduce(recalBAM, reducedBAM))

      // collect QC metrics based on full BAM
      val outGcBiasMetrics = swapExt(recalBAM,".bam",".gc_metrics")
      val outMultipleMetrics = swapExt(recalBAM,".bam",".multipleMetrics")
      add(calculateGCMetrics(recalBAM, outGcBiasMetrics))
      add(calculateMultipleMetrics(recalBAM, outMultipleMetrics))

      val normalName = tumorInfo.filter( P => P._2 == 0).keysIterator.next()
      val tumorName = tumorInfo.filter( P => P._2 == 1).keysIterator.next()

      // todo hotfix: there is not a way right now to look up the BAM for a sample?
      val normalBam = normalName + ".clean.dedup.recal.bam"
      val tumorBam = tumorName + ".clean.dedup.recal.bam"
      qscript.unreducedNormalBAM = new File(normalBam)
      qscript.unreducedTumorBAM = new File(tumorBam)
      qscript.unreducedNormalBAMIndex = swapExt(qscript.unreducedNormalBAM,".bam",".bai")
      qscript.unreducedTumorBAMIndex = swapExt(qscript.unreducedTumorBAMIndex,".bam",".bai")
      qscript.reducedNormalBAM = swapExt(normalBam,".bam",".reduced.bam")
      qscript.reducedTumorBAM = swapExt(tumorBam,".bam",".reduced.bam")
      qscript.reducedNormalBAMIndex = swapExt(qscript.reducedNormalBAM,".bam",".bai")
      qscript.reducedTumorBAMIndex = swapExt(qscript.reducedTumorBAMIndex,".bam",".bai")

      // add single sample vcf germline calling
      if (qscript.doSingleSampleCalling) {
        add(call(reducedBAM, outVCF))
        // report output parameters
        qscript.singleSampleVCF = outVCF
        qscript.singleSampleVCFIndex = outVCF + ".idx"
      }
    }
    // todo hotfix part 3: add cancer-specific calling also at the end
    if (qscript.doMutect) {

      val tumorFractionContamination:Float = 0.01f
      val outPrefix = tumorName + "-vs-" + normalName
      val rawMutations = outPrefix + ".call_stats.txt"
      val rawVcf = outPrefix + ".vcf"
      val rawCoverage = outPrefix + ".wig.txt"
      print("MUTECT Tumor BAM is: " + tumorBam.getName)
      add(mutect(tumorName, tumorBam, normalName, normalBam, tumorFractionContamination, rawMutations, rawVcf, rawCoverage))

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
    val targetIntervals = swapExt(bam, ".bam", ".cleaning.interval_list")
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
    if (indelSites != null)
      this.known ++= qscript.indelSites
    this.scatterCount = nContigs
    this.analysisName = outIntervals + ".target"
    this.jobName = outIntervals + ".target"
    this.intervals :+= qscript.targets // todo hotfix - only run in given targets, otherwise we're just wasting computation
  }

  case class indel (inBAMs: Seq[File], tIntervals: File) extends IndelRealigner with CommandLineGATKArgs {
    // TODO -- THIS IS A WORKAROUND FOR QUEUE'S LIMITATION OF TRACKING LISTS OF FILES (implementation limited to 5 files)
    @Output(doc="first cleaned bam file", required=true) var out1: File = swapExt(inBAMs(0), ".bam", cleaningExtension)
    @Output(doc="first cleaned bam file", required=true) var ind1: File = swapExt(out1, ".bam", ".bai")
    @Output(doc="first cleaned bam file", required=false) var out2: File = if (inBAMs.length >= 2) {swapExt(inBAMs(1), ".bam", cleaningExtension)} else {null}
    @Output(doc="first cleaned bam file", required=false) var ind2: File = if (inBAMs.length >= 2) {swapExt(out2, ".bam", ".bai")} else {null}
    @Output(doc="first cleaned bam file", required=false) var out3: File = if (inBAMs.length >= 3) {swapExt(inBAMs(2), ".bam", cleaningExtension)} else {null}
    @Output(doc="first cleaned bam file", required=false) var ind3: File = if (inBAMs.length >= 3) {swapExt(out3, ".bam", ".bai")} else {null}
    @Output(doc="first cleaned bam file", required=false) var out4: File = if (inBAMs.length >= 4) {swapExt(inBAMs(3), ".bam", cleaningExtension)} else {null}
    @Output(doc="first cleaned bam file", required=false) var ind4: File = if (inBAMs.length >= 4) {swapExt(out4, ".bam", ".bai")} else {null}
    @Output(doc="first cleaned bam file", required=false) var out5: File = if (inBAMs.length >= 5) {swapExt(inBAMs(4), ".bam", cleaningExtension)} else {null}
    @Output(doc="first cleaned bam file", required=false) var ind5: File = if (inBAMs.length >= 5) {swapExt(out5, ".bam", ".bai")} else {null}
    this.input_file = inBAMs
    this.targetIntervals = tIntervals

    // FIXME - nWayOut doesn't seem to work, for now really only support single sample BAMs


    if (inBAMs.length == 1) {
      this.o = out1
    } else {
      this.nWayOut = cleaningExtension
    }
    this.known ++= qscript.dbSNP
    if (qscript.indelSites != null)
      this.known ++= qscript.indelSites
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
    this.analysisName = outBAM + ".reduce"
    this.jobName = outBAM + ".reduce"
    this.intervals :+= qscript.targets // TODO hotfix - just so we don't waste 50 hours to completion on a 10kb BAM

  }

  case class call (inBAM: File, outVCF: File) extends UnifiedGenotyper with CommandLineGATKArgs {
    this.input_file :+= inBAM
    this.out = outVCF
    this.isIntermediate = false
    this.analysisName = outVCF + ".singleSampleCalling"
    this.jobName = outVCF + ".singleSampleCalling"
    //this.dbsnp = qscript.dbSNP.
    this.downsample_to_coverage = 600
    this.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.scatterCount = nContigs

     // for now have the caller just traverse the whole genome even though it will waste 99.999% of its time
    this.intervals :+= qscript.targets // TODO hotfix - just so we don't waste 50 hours to completion on a 10kb BAM
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
    this.assumeSorted = Some(true)
  }

  case class calculateHSMetrics (inBAM:File, outFile: File) extends CalculateHsMetrics with ExternalCommonArgs {
    this.reference = qscript.reference
    this.input :+= inBAM
    this.output = outFile
    this.targets = qscript.targets
    this.baits = qscript.baits
    this.analysisName = outFile + ".hsMetrics"
    this.jobName = outFile + ".hsMetrics"
    // todo - do we want to compute per-read group HS metrics?

  }

  case class calculateGCMetrics (inBAM:File, outFile: File) extends CollectGcBiasMetrics with ExternalCommonArgs {
    this.reference = qscript.reference
    this.input :+= inBAM
    this.output = outFile
    this.analysisName = outFile + ".gcMetrics"
    this.jobName = outFile + ".gcMetrics"
  }

  case class calculateMultipleMetrics (inBAM:File, outFile: File) extends CollectMultipleMetrics with ExternalCommonArgs {
    this.reference = qscript.reference
    this.input :+= inBAM
    this.output = outFile
    this.analysisName = outFile + ".metrics"
    this.jobName = outFile + ".metrics"
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
    this.memoryLimit = 1
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

  // todo hotfix - pasted from cancer pipeline

  case class mutect (tumorName : String, tumorBam : File, normalName : String, normalBam : File, tumorFractionContamination : Float, outMutations : File,  outVcf : File, outCoverage : File) extends MuTect with CommandLineGATKArgs {
    this.isIntermediate = false
    this.scatterCount = 1
    this.memoryLimit = 4
    this.jarFile = qscript.mutectJar
    this.intervals :+= qscript.targets

    this.dbsnp = qscript.dbSNP
    this.cosmic = qscript.cosmic
    this.normal_panel = qscript.pon

    this.only_passing_calls = true
    this.enable_extended_output = true
    this.downsample_to_coverage = 1000 // TODO: how deep should this be?
    this.fraction_contamination = Some(tumorFractionContamination)

    this.input_file :+= new TaggedFile(tumorBam, "tumor")
    this.input_file :+= new TaggedFile(normalBam, "normal")

    this.out = outMutations
    this.coverage_file = outCoverage
    this.vcf = outVcf


    this.analysisName = tumorBam.toString + ".mutect"
    this.jobName = this.analysisName

    print ("MuTect CMD:" + this.commandLine)
  }

}

