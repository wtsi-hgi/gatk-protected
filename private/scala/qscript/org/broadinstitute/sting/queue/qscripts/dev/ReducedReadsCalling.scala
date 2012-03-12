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

import org.apache.commons.io.FilenameUtils
import org.broadinstitute.sting.pipeline.PicardAggregationUtils
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.gatk.walkers.compression.reducereads.ReduceReadsWalker.DownsampleStrategy
import io.Source._
import collection.JavaConversions._

class ReducedReadsCalling extends QScript {
  qscript =>

  @Input(doc="BAM list files. Name must be <projectName>.bam.list", shortName="I", required=false)
  var originalBamList: File = _

  @Input(doc="GATK or Picard intervals file.", shortName="L", required=false)
  var intervals: File = _

  @Input(doc="Level of parallelism for UnifiedGenotyper. By default set to 20.", shortName="varScatter", required=false)
  var variantCallerScatterCount = 20

  @Argument(doc="Pipeline memory limit. By default set to 2g.", shortName="pipeMemory", required=false)
  var pipelineMemoryLimit = 2

  @Argument(doc="Reduced reads memory limit. By default set to 6g.", shortName="rrMemory", required=false)
  var reducedMemoryLimit = 6

  @Argument(doc="Expand each target in input intervals by the specified number of bases. By default set to 50 bases.", shortName="expand", required=false)
  var expandIntervals = 50

  @Argument(doc="Subdirectory to store the reduced bams. By default set to 'reducedBams'.", shortName="bamDir", required=false)
  var bamDir = "reducedBams"

  def script() {
    val exomeIntervals = new File("/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list")
    val resources = "/humgen/gsa-pipeline/resources/b37/v3/"
    val k1gTrainingHighQuality = resources + "phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"
    val k1gTrainingTerrible = resources + "phase1.wgs.projectConsensus.v2b.recal.terrible.vcf"
    val omni = resources + "Omni25_sites_1525_samples.b37.vcf"
    val hapmap = resources + "hapmap_3.3.b37.sites.vcf"
    val dbsnp129 = resources + "dbsnp_132.b37.excluding_sites_after_129.vcf"
    val dbsnp132 = resources + "dbsnp_132.b37.vcf"

    var reference = resources + "human_g1k_v37.fasta"

    var projectName: String = null

    require(originalBamList != null && originalBamList.getName.endsWith(".bam.list"), "-I/--bamList must be specified as <projectName>.bam.list")
    require(intervals != null, "-L/--intervals must be specified")
    projectName = originalBamList.getName.stripSuffix(".bam.list")

    qSettings.runName = projectName

    val flankIntervals = projectName + ".flanks.intervals"

    if (qscript.expandIntervals > 0) {
      val writeFlanks = new WriteFlankingIntervalsFunction
      writeFlanks.reference = reference
      writeFlanks.inputIntervals = intervals
      writeFlanks.flankSize = qscript.expandIntervals
      writeFlanks.outputIntervals = flankIntervals
      writeFlanks.jobOutputFile = writeFlanks.outputIntervals + ".out"
      add(writeFlanks)
    }

    trait CommandLineGATKArgs extends CommandLineGATK {
      this.reference_sequence = reference
      this.intervals = Seq(qscript.intervals)
      this.memoryLimit = pipelineMemoryLimit
    }

    trait ExpandedIntervals extends CommandLineGATK {
      if (qscript.expandIntervals > 0)
        this.intervals :+= flankIntervals
    }

    var reducedBamList: List[File] = List()

    for (file: String <- fromFile(originalBamList).getLines()) {
      val reducedBAM = swapExt(file.getName, ".bam", ".reduced.bam")
      val reducedBAMFile = new File(bamDir + "/" + reducedBAM)
      reducedBamList :+= reducedBAMFile

      // reduce
      val reduce = new ReduceReads() with CommandLineGATKArgs with ExpandedIntervals
      reduce.memoryLimit = reducedMemoryLimit
      reduce.input_file :+= new File(file)
      reduce.out = reducedBAMFile
      add(reduce)
    }

    // TODO -- If we encounter too much overhead from having so many individual files,
    // TODO --   first merge these files together into batches (of 100?)

    val writeBamList = new ListWriterFunction
    writeBamList.inputFiles = reducedBamList
    writeBamList.listFile = projectName + ".reduced.bam.list"
    add(writeBamList)

    val call = new UnifiedGenotyper with CommandLineGATKArgs with ExpandedIntervals
    call.input_file = Seq(writeBamList.listFile)
    call.dbsnp = dbsnp132
    call.downsample_to_coverage = 600
    call.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    call.out = projectName + ".unfiltered.vcf"
    call.jobOutputFile = call.out + ".out"
    call.scatterCount = qscript.variantCallerScatterCount
    add(call)

    val selectSNPs = new SelectVariants with CommandLineGATKArgs with ExpandedIntervals
    selectSNPs.selectType :+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.SNP
    selectSNPs.variant = call.out
    selectSNPs.out = projectName + ".snps.unfiltered.vcf"
    selectSNPs.jobOutputFile = selectSNPs.out + ".out"
    add(selectSNPs)

    val selectIndels = new SelectVariants with CommandLineGATKArgs with ExpandedIntervals
    selectIndels.selectType :+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.INDEL
    selectIndels.variant = call.out
    selectIndels.out = projectName + ".indels.unfiltered.vcf"
    selectIndels.jobOutputFile = selectIndels.out + ".out"
    add(selectIndels)

    val trancheLevels = Seq(
      "100.0", "99.9", "99.5", "99.3",
      "99.0", "98.9", "98.8",
      "98.5", "98.4", "98.3", "98.2", "98.1",
      "98.0", "97.9", "97.8",
      "97.5",
      "97.0",
      "95.0",
      "90.0")

    // TODO -- This step can be really slow for big VCF files (many samples = many genotypes);
    // TODO --   maybe we should run VQSR in parallel?

    val buildSNPModel = new VariantRecalibrator with CommandLineGATKArgs with ExpandedIntervals
    buildSNPModel.input :+= selectSNPs.out
    buildSNPModel.resource :+= TaggedFile(hapmap, "training=true,truth=true,prior=15.0")
    buildSNPModel.resource :+= TaggedFile(omni, "training=true,truth=true,prior=12.0")
    buildSNPModel.resource :+= TaggedFile(k1gTrainingHighQuality, "training=true,prior=12.0")
    buildSNPModel.resource :+= TaggedFile(k1gTrainingTerrible, "bad=true,prior=2.0")
    buildSNPModel.resource :+= TaggedFile(dbsnp129, "known=true,prior=4.0")
    buildSNPModel.use_annotation = Seq("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "InbreedingCoeff")
    buildSNPModel.trustAllPolymorphic = true
    buildSNPModel.maxGaussians = 6
    buildSNPModel.stdThreshold = 14
    buildSNPModel.percentBadVariants = 0.03
    buildSNPModel.TStranche = trancheLevels
    buildSNPModel.tranches_file = projectName + ".snps.tranches"
    buildSNPModel.recal_file = projectName + ".snps.recal"
    buildSNPModel.jobOutputFile = buildSNPModel.recal_file + ".out"
    add(buildSNPModel)

    val applySNPModel = new ApplyRecalibration with CommandLineGATKArgs with ExpandedIntervals
    applySNPModel.input :+= selectSNPs.out
    applySNPModel.tranches_file = buildSNPModel.tranches_file
    applySNPModel.recal_file = buildSNPModel.recal_file
    applySNPModel.ts_filter_level = 98.5
    val filteredSNPsVcf = projectName + ".snps.recalibrated.vcf"
    applySNPModel.out = filteredSNPsVcf
    applySNPModel.jobOutputFile = applySNPModel.out + ".out"
    add(applySNPModel)

    val filterIndels = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
    filterIndels.variant = selectIndels.out
    filterIndels.filterName = Seq("Indel_FS", "Indel_QD", "Indel_ReadPosRankSum", "Indel_InbreedingCoeff")
    filterIndels.filterExpression = Seq("FS>200.0", "QD<2.0", "ReadPosRankSum<-20.0", "InbreedingCoeff<-0.8")
    filterIndels.out = projectName + ".indels.filtered.vcf"
    filterIndels.jobOutputFile = filterIndels.out + ".out"
    add(filterIndels)

    val combineSNPsIndels = new CombineVariants with CommandLineGATKArgs with ExpandedIntervals
    combineSNPsIndels.variant :+= TaggedFile(filterIndels.out, "indels")
    combineSNPsIndels.variant :+= TaggedFile(filteredSNPsVcf, "snps")
    combineSNPsIndels.filteredrecordsmergetype = org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    combineSNPsIndels.assumeIdenticalSamples = true
    combineSNPsIndels.out = projectName + ".unannotated.vcf"
    combineSNPsIndels.jobOutputFile = combineSNPsIndels.out + ".out"
    add(combineSNPsIndels)

    val eval = new VariantEval with CommandLineGATKArgs
    eval.eval :+= combineSNPsIndels.out
    eval.dbsnp = dbsnp129
    eval.doNotUseAllStandardModules = true
    eval.evalModule = Seq("TiTvVariantEvaluator", "CountVariants", "CompOverlap", "MultiallelicSummary")
    eval.doNotUseAllStandardStratifications = true
    eval.stratificationModule = Seq("EvalRod", "CompRod", "Novelty")
    eval.out = projectName + ".eval"
    eval.jobOutputFile = eval.out + ".out"
    add(eval)
  }
}
