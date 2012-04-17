/*
 * Copyright (c) 2012, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.interval.IntervalUtils
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model

class WholeGenomePipeline extends QScript {
  @Input(doc="Bam file list", shortName = "I", required=true)
  var bamList: File = _

  @Input(doc="If set will skip indel realignment", shortName = "skipClean", required=false)
  var skipClean = false

  @Input(doc="Exclude intervals list", shortName = "XL", required=false)
  var excludeIntervals: Seq[File] = Nil

  @Argument(doc="path to tmp space for storing intermediate bam files", shortName="cleanerTmpDir", required=true)
  var cleanerTmpDir: String = _

  @Argument(doc="Flag for running the whole genome (wg) or chromosome 20 (chr20). The default is chr20.", shortName="runType", required=false)
  var runType = "chr20"

  @Argument(doc="Chunk size. Defaults to 250,000", shortName="chunk", required=false)
  var chunkSize = 250000

  @Argument(doc="Standard memory limit. Defaults to 4g", shortName="pipeMem", required=false)
  var pipelineMemoryLimit = 4

  @Argument(doc="Memory limit for VQSR. Defaults to 32g", shortName="vqsrMem", required=false)
  var vqsrMemoryLimit = 32

  @Argument(doc="Memory limit for variant calling and indel realignment. Defaults to 4g", shortName="callMem", required=false)
  var callMemoryLimit = 4

  @Argument(doc="Memory reservation for variant calling and indel realignment. Defaults to caller memory limit", shortName="callRes", required=false)
  var callMemoryReservation: Option[Double] = None

  def script() {
    // TODO: examine dontRequestMultipleCores as a global variable versus a local
    // -nit uses a background thread to avoid 1024 ulimit, but the thread 1-10% of a CPU.
    this.qSettings.dontRequestMultipleCores = true

    var runIntervals = Seq.empty[String]
    var intervals = Traversable.empty[Interval]

    var reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
    val resources = "/humgen/gsa-pipeline/resources/b37/v4/"
    val dbsnp129 = resources + "dbsnp_135.b37.excluding_sites_after_129.vcf"
    val dbsnp135 = resources + "dbsnp_135.b37.vcf"
    val hapmap = resources + "hapmap_3.3.b37.sites.vcf"
    val k1gOmni = resources + "1000G_omni2.5.b37.sites.vcf"
    val k1gPhase1Indels = resources + "1000G_phase1.indels.b37.vcf"
    val millsK1gIndels = resources + "Mills_and_1000G_gold_standard.indels.b37.sites.vcf"

    trait CommandLineGATKArgs extends CommandLineGATK {
      this.reference_sequence = reference
      this.intervalsString = runIntervals
      this.memoryLimit = pipelineMemoryLimit
    }

    trait CallerGATKArgs extends CommandLineGATKArgs {
      this.residentRequest = callMemoryReservation
      this.memoryLimit = callMemoryLimit
      this.javaMemoryLimit = callMemoryLimit
      this.residentLimit = callMemoryLimit * 1.2
    }

    trait VQSRGATKArgs extends CommandLineGATKArgs {
      this.memoryLimit = vqsrMemoryLimit
    }

    case class Interval(chr: String, start: Long, stop: Long) {
      override def toString = "%s:%d-%d".format(chr, start, stop)
    }

    runType = runType.toLowerCase
    val contigSizes = IntervalUtils.getContigSizes(reference)
    if (runType == "wg") {
      val contigs = (1 to 22).map(_.toString) ++ Seq("X", "Y", "MT")
      intervals = contigs.map(chr => new Interval(chr, 1, contigSizes.get(chr).longValue))
      runIntervals = Nil
    } else {
      val locs = Map(
        "cent1" -> new Interval("1", 121429168, 121529168),
        "cent16" -> new Interval("16", 40844464, 40944464),
        "chr20" -> new Interval("20", 1, contigSizes.get("20").longValue),
        "chr20_100k" -> new Interval("20", 100001, 200000))

      locs.get(runType) match {
        case Some(range) =>
          intervals = Seq(range)
          runIntervals = Seq(range.toString)
        case None =>
          throw new RuntimeException("Invalid runType: " + runType + ". Must be one of: " + locs.keys.mkString(", ") + ", or wg")
      }
    }

    require(bamList != null && bamList.getName.endsWith(".bam.list"), "-I/--bamList must be specified as <projectName>.bam.list")
    val project = bamList.getName.stripSuffix(".bam.list")

    var snpChrVcfs = Seq.empty[File]
    var indelChrVcfs = Seq.empty[File]
    val glModels = Seq(Model.SNP, Model.INDEL)

    for (interval <- intervals) {
      var snpChunkVcfs = Seq.empty[File]
      var indelChunkVcfs = Seq.empty[File]

      val chr = interval.chr
      val chrStart = interval.start
      val chrStop = interval.stop
      val chrSize = interval.stop - interval.start + 1
      val chrDir = "chrs/" + chr + "/"
      val chrBase = project + ".chr" + chr

      var start = chrStart
      var chunkNumber = 1
      val chunkCount = (chrSize / chunkSize) + (if (chrSize % chunkSize == 0) 0 else 1)
      val chunkFormat = "%s_chunk_%%0%dd_of_%d".format(chrBase, chunkCount.toString.length(), chunkCount)

      while (start <= chrStop) {
        val stop = (start + chunkSize - 1) min chrStop

        val chunkBase: String = chrDir + "chunks/" + chunkFormat.format(chunkNumber)
        val tmpBase: String = cleanerTmpDir + "/" + chunkBase

        val chunkInterval = Seq("%s:%d-%d".format(chr, start, stop))

        var call_input_file = bamList

        if (!skipClean) {
          val target = new RealignerTargetCreator with CallerGATKArgs
          target.input_file :+= bamList
          target.intervalsString = chunkInterval
          target.excludeIntervals = excludeIntervals
          target.known :+= k1gPhase1Indels
          target.known :+= millsK1gIndels
          //target.num_threads = 2
          //target.num_io_threads = 2 // Avoid 1024 ulimits
          target.out = tmpBase + ".target.intervals"
          target.isIntermediate = true
          add(target)

          val clean = new IndelRealigner with CallerGATKArgs
          clean.input_file :+= bamList
          clean.intervalsString = chunkInterval
          clean.excludeIntervals = excludeIntervals
          clean.targetIntervals = target.out
          clean.knownAlleles :+= k1gPhase1Indels
          clean.knownAlleles :+= millsK1gIndels
          clean.simplifyBAM = true
          //clean.num_threads = 2
          //clean.num_io_threads = 2 // Avoid 1024 ulimits
          clean.out = tmpBase + ".cleaned.bam"
          clean.isIntermediate = true
          add(clean)

          call_input_file = clean.out
        }

        for (glModel <- glModels) {
          val call = new UnifiedGenotyper with CallerGATKArgs
          call.input_file :+= call_input_file
          call.intervalsString = chunkInterval
          call.excludeIntervals = excludeIntervals
          call.dbsnp = dbsnp135
          call.downsample_to_coverage = 50
          call.standard_min_confidence_threshold_for_calling = 4.0
          call.standard_min_confidence_threshold_for_emitting = 4.0
          call.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
          call.genotype_likelihoods_model = glModel
          call.out = tmpBase + "." + glModel.toString.toLowerCase + "s.vcf"
          call.isIntermediate = true
          add(call)

          glModel match {
            case Model.SNP =>
              snpChunkVcfs :+= call.out
            case Model.INDEL =>
              call.max_alternate_alleles = 2 // reduce processing time with fewer alts
              indelChunkVcfs :+= call.out
          }
        }

        start += chunkSize
        chunkNumber += 1
      }

      for (glModel <- glModels) {
        val chunkVcfs = glModel match {
          case Model.SNP => snpChunkVcfs
          case Model.INDEL => indelChunkVcfs
        }

        val combineChunks = new CombineVariants with CommandLineGATKArgs
        combineChunks.variant = chunkVcfs.zipWithIndex.map { case (vcf, index) => TaggedFile(vcf, "input"+index) }
        combineChunks.rod_priority_list = chunkVcfs.zipWithIndex.map { case (vcf, index) => "input"+index }.mkString(",")
        combineChunks.assumeIdenticalSamples = true
        combineChunks.out = chrDir + chrBase + "." + glModel.toString.toLowerCase + "s.unfiltered.vcf"
        add(combineChunks)

        glModel match {
          case Model.SNP =>
            snpChrVcfs :+= combineChunks.out
          case Model.INDEL =>
            indelChrVcfs :+= combineChunks.out
        }
      }
    }

    val tranche = 98.5
    val tranches = Seq(
      "100.0", "99.9", "99.5", "99.3",
      "99.0", "98.9", "98.8",
      "98.5", "98.4", "98.3", "98.2", "98.1",
      "98.0",
      "97.0",
      "95.0",
      "90.0")
    var evalVcfs = Seq.empty[File]

    for (glModel <- glModels) {
      val buildModel = new VariantRecalibrator with VQSRGATKArgs
      buildModel.trustAllPolymorphic = true
      buildModel.TStranche = tranches
      buildModel.recal_file = project + "." + glModel.toString.toLowerCase + "s.recal"
      buildModel.tranches_file = project + "." + glModel.toString.toLowerCase + "s.tranches"
      add(buildModel)

      glModel match {
        case Model.SNP =>
          buildModel.input = snpChrVcfs
          buildModel.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
          buildModel.use_annotation = Seq("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "InbreedingCoeff", "DP")
          buildModel.resource :+= TaggedFile(hapmap, "training=true,truth=true,prior=15.0")
          buildModel.resource :+= TaggedFile(k1gOmni, "training=true,prior=12.0")
          buildModel.resource :+= TaggedFile(dbsnp135, "known=true,prior=6.0")
        case Model.INDEL =>
          buildModel.input = indelChrVcfs
          buildModel.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
          buildModel.use_annotation = Seq("QD", "FS", "HaplotypeScore", "ReadPosRankSum", "InbreedingCoeff")
          buildModel.resource :+= TaggedFile(millsK1gIndels, "known=true,training=true,truth=true,prior=12.0")
      }


      for (chrVcf <- buildModel.input) {
        val applyRecalibration = new ApplyRecalibration with VQSRGATKArgs
        applyRecalibration.input :+= chrVcf
        applyRecalibration.recal_file = buildModel.recal_file
        applyRecalibration.tranches_file = buildModel.tranches_file
        applyRecalibration.ts_filter_level = tranche
        applyRecalibration.out = swapBaseExt(chrVcf, ".unfiltered.vcf", ".recalibrated.tranche_" + tranche.toString.replace(".", "_") + ".vcf")
        add(applyRecalibration)

        evalVcfs :+= applyRecalibration.out
      }
    }

    for (strats <- Seq(
      Seq("AlleleCount"),
      Seq("Sample")
    )) {
      val eval = new VariantEval with CommandLineGATKArgs
      eval.eval = evalVcfs
      eval.mergeEvals = true
      eval.dbsnp = dbsnp129
      eval.doNotUseAllStandardModules = true
      eval.evalModule = Seq("TiTvVariantEvaluator", "CountVariants", "CompOverlap")
      eval.doNotUseAllStandardStratifications = true
      eval.stratificationModule = Seq("EvalRod", "CompRod", "Novelty", "FunctionalClass") ++ strats
      eval.out = project + strats.map(_.toLowerCase).mkString(".by_", "_", "") + ".eval"
      eval

      add(eval)
    }
  }

  private def swapBaseExt(file: File, oldExt: String, newExt: String): File =
    file.getParentFile + "/" + swapExt(file, oldExt, newExt)
}
