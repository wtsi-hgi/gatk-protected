/*
 * Copyright (c) 2011, The Broad Institute
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

import org.broadinstitute.sting.pipeline.{PicardPipeline, Pipeline}
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.queue.library.ipf.intervals.ExpandIntervals
import org.broadinstitute.sting.queue.QScript
import collection.JavaConversions._

class HybridSelectionPipeline extends QScript {
  qscript =>

  private final val FILTER_VQSR = "VQSR"
  private final val FILTER_HARD = "HARD"

  @Argument(doc="the YAML file specifying inputs, interval lists, reference sequence, etc.", shortName="Y")
  var yamlFile: File = _

  @Input(doc="level of parallelism for UnifiedGenotyper. By default set to 20.", shortName="varScatter", required=false)
  var variantCallerScatterCount = 20

  @Argument(doc="memory limit for UnifiedGenotyper. By default set to 2g.", shortName="varMemory", required=false)
  var variantCallerMemory = 2

  @Argument(doc="memory limit for VariantRecalibrator and ApplyRecalibration. By default set to 8g.", shortName="vqsrMemory", required=false)
  var variantRecalibratorMemory = 8

  @Argument(doc="expand each target in input intervals by the specified number of bases. By default set to 50 bases.", shortName="expand", required=false)
  var expandIntervals = 50

  @Argument(doc="pipeline memory limit. By default set to 2g.", shortName="pipeMemory", required=false)
  var pipelineMemoryLimit = 2

  @Argument(doc="variant filter type, " + FILTER_VQSR + " or " + FILTER_HARD + ". By default hard filters are used unless the bait set is whole_exome_agilent_1.", shortName="varFilter", required=false)
  var variantFilterType: String = _

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = qscript.pipeline.getProject.getReferenceFile
    this.intervals = List(qscript.pipeline.getProject.getIntervalList)
    this.memoryLimit = pipelineMemoryLimit
  }

  val resources = "/humgen/gsa-pipeline/resources/b37/v2/"
  val exomeIntervals = "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"
  val k1gExomesBam = resources + "1000_Genomes_Whole_Exome_50_Samples.bam"
  val k1gExomesSamples = resources + "1000_Genomes_Whole_Exome_50_Samples.samples"
  val k1gTrainingHighQuality = resources + "phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"
  val k1gTrainingTerrible = resources + "phase1.wgs.projectConsensus.v2b.recal.terrible.vcf"
  val dbsnp129 = resources + "dbsnp_129.b37.vcf"
  val omni = resources + "Omni25_sites_1525_samples.b37.vcf"
  val hapmap = resources + "hapmap_3.3.b37.sites.vcf"

  def script() {
    pipeline = PicardPipeline.parse(qscript.yamlFile)

    val projectBase = qscript.pipeline.getProject.getName
    val bamType = "cleaned"

    val writeBamList = new ListWriterFunction
    writeBamList.inputFiles = qscript.pipeline.getSamples.filter(_.getBamFiles.contains(bamType)).map(_.getBamFiles.get(bamType)).toList
    writeBamList.listFile = projectBase +".bam.list"
    writeBamList.jobOutputFile = writeBamList.listFile + ".out"
    add(writeBamList)

    val isExomeIntervals = (qscript.pipeline.getProject.getIntervalList.getAbsolutePath == exomeIntervals)

    if (variantFilterType == null) {
      variantFilterType = if (isExomeIntervals) FILTER_VQSR else FILTER_HARD
    } else {
      variantFilterType = variantFilterType.toUpperCase
    }

    val useK1gExomes = (qscript.pipeline.getSamples.size() < 50 && isExomeIntervals && variantFilterType == FILTER_VQSR)

    val flankIntervals = projectBase + ".flanks.intervals"

    if (qscript.expandIntervals > 0) {
      val ei = new ExpandIntervals(
        qscript.pipeline.getProject.getIntervalList,
        1,
        qscript.expandIntervals,
        flankIntervals,
        qscript.pipeline.getProject.getReferenceFile,
        "INTERVALS",
        "INTERVALS")
      ei.jobOutputFile = ei.outList + ".out"

      add(ei)
    }

    trait ExpandedIntervals extends CommandLineGATK {
      if (qscript.expandIntervals > 0)
        this.intervals :+= flankIntervals
    }

    val call = new UnifiedGenotyper with CommandLineGATKArgs with ExpandedIntervals
    call.input_file = List(writeBamList.listFile)
    if (useK1gExomes)
      call.input_file :+= k1gExomesBam
    call.dbsnp = qscript.pipeline.getProject.getGenotypeDbsnp
    call.downsample_to_coverage = 600
    call.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    call.out = projectBase + ".unfiltered.vcf"
    call.jobOutputFile = call.out + ".out"
    call.scatterCount = qscript.variantCallerScatterCount
    call.memoryLimit = qscript.variantCallerMemory
    add(call)

    val selectSNPs = new SelectVariants with CommandLineGATKArgs with ExpandedIntervals
    selectSNPs.selectType :+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.SNP
    selectSNPs.variant = call.out
    selectSNPs.out = projectBase + ".snps.unfiltered.vcf"
    selectSNPs.jobOutputFile = selectSNPs.out + ".out"
    add(selectSNPs)

    val selectIndels = new SelectVariants with CommandLineGATKArgs with ExpandedIntervals
    selectIndels.selectType :+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.INDEL
    selectIndels.variant = call.out
    selectIndels.out = projectBase + ".indels.unfiltered.vcf"
    selectIndels.jobOutputFile = selectIndels.out + ".out"
    add(selectIndels)

    val filteredSNPsVcf = variantFilterType match {
      case FILTER_VQSR =>
        val trancheLevels = List(
          "100.0", "99.9", "99.5", "99.3",
          "99.0", "98.9", "98.8",
          "98.5", "98.4", "98.3", "98.2", "98.1",
          "98.0", "97.9", "97.8",
          "97.5",
          "97.0",
          "95.0",
          "90.0")

        val buildSNPModel = new VariantRecalibrator with CommandLineGATKArgs with ExpandedIntervals
        buildSNPModel.input :+= selectSNPs.out
        buildSNPModel.training :+= TaggedFile(hapmap, "prior=15.0")
        buildSNPModel.training :+= TaggedFile(omni, "prior=12.0")
        buildSNPModel.training :+= TaggedFile(k1gTrainingHighQuality, "prior=12.0")
        buildSNPModel.truth :+= TaggedFile(hapmap, "prior=15.0")
        buildSNPModel.truth :+= TaggedFile(omni, "prior=12.0")
        buildSNPModel.badSites :+= TaggedFile(k1gTrainingTerrible, "prior=2.0")
        buildSNPModel.known :+= TaggedFile(dbsnp129, "prior=4.0")
        buildSNPModel.use_annotation = List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "InbreedingCoeff")
        buildSNPModel.trustAllPolymorphic = true
        buildSNPModel.maxGaussians = 6
        buildSNPModel.stdThreshold = 14
        buildSNPModel.percentBadVariants = 0.03
        buildSNPModel.TStranche = trancheLevels
        buildSNPModel.tranches_file = projectBase + ".tranches"
        buildSNPModel.recal_file = projectBase + ".recal"
        buildSNPModel.jobOutputFile = buildSNPModel.recal_file + ".out"
        buildSNPModel.memoryLimit = qscript.variantRecalibratorMemory
        add(buildSNPModel)

        val applySNPModel = new ApplyRecalibration with CommandLineGATKArgs with ExpandedIntervals
        applySNPModel.input :+= selectSNPs.out
        applySNPModel.tranches_file = buildSNPModel.tranches_file
        applySNPModel.recal_file = buildSNPModel.recal_file
        applySNPModel.ts_filter_level = 98.5
        applySNPModel.out = projectBase + ".snps.recalibrated.vcf"
        applySNPModel.jobOutputFile = applySNPModel.out + ".out"
        applySNPModel.memoryLimit = qscript.variantRecalibratorMemory
        add(applySNPModel)

        applySNPModel.out

      case FILTER_HARD =>
        val filterSNPs = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
        filterSNPs.variant = selectSNPs.out
        filterSNPs.filterName = List("SNP_QD", "SNP_MQ", "SNP_HaplotypeScore", "SNP_MQRankSum", "SNP_ReadPosRankSum", "SNP_FS")
        filterSNPs.filterExpression = List("\"QD<2.0\"", "\"MQ<40.0\"", "\"HaplotypeScore>13.0\"", "\"MQRankSum<-12.5\"", "\"ReadPosRankSum<-8.0\"", "\"FS>60.0\"")
        filterSNPs.out = projectBase + ".snps.filtered.vcf"
        filterSNPs.jobOutputFile = filterSNPs.out + ".out"
        add(filterSNPs)

        filterSNPs.out

      case unknown =>
        throw new IllegalArgumentException("""If set the variantFilterType must be "%s" or "%s"""".format(FILTER_VQSR, FILTER_HARD))
    }

    val filterIndels = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
    filterIndels.variant = selectIndels.out
    filterIndels.filterName = List("Indel_FS", "Indel_QD", "Indel_ReadPosRankSum", "Indel_InbreedingCoeff")
    filterIndels.filterExpression = List("\"FS>200.0\"", "\"QD<2.0\"", "\"ReadPosRankSum<-20.0\"", "\"InbreedingCoeff<-0.8\"")
    filterIndels.out = projectBase + ".indels.filtered.vcf"
    filterIndels.jobOutputFile = filterIndels.out + ".out"
    add(filterIndels)

    val combineSNPsIndels = new CombineVariants with CommandLineGATKArgs with ExpandedIntervals
    combineSNPsIndels.variant :+= TaggedFile(filterIndels.out, "indels")
    combineSNPsIndels.variant :+= TaggedFile(filteredSNPsVcf, "snps")
    combineSNPsIndels.filteredrecordsmergetype = org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    combineSNPsIndels.assumeIdenticalSamples = true
    combineSNPsIndels.out = projectBase + ".unannotated.vcf"
    combineSNPsIndels.jobOutputFile = combineSNPsIndels.out + ".out"
    add(combineSNPsIndels)

    val snpsIndelsVcf =
      if (!useK1gExomes) {
        combineSNPsIndels.out
      } else {
        val removeK1gSamples = new SelectVariants with CommandLineGATKArgs with ExpandedIntervals
        removeK1gSamples.variant = combineSNPsIndels.out
        removeK1gSamples.exclude_sample_file :+= k1gExomesSamples
        removeK1gSamples.out = projectBase + ".selected.vcf"
        removeK1gSamples.jobOutputFile = removeK1gSamples.out + ".out"
        add(removeK1gSamples)

        removeK1gSamples.out
      }

    //
    // TODO -- David will replace the Genomic Annotator with snpEff
    //

    //val annotate = new GenomicAnnotator with CommandLineGATKArgs with ExpandedIntervals
    //annotate.rodBind :+= RodBind("variant", "VCF", snpsIndelsVcf)
    //annotate.rodBind :+= RodBind("refseq", "AnnotatorInputTable", qscript.pipeline.getProject.getRefseqTable)
    //annotate.rodToIntervalTrackName = "variant"
    //annotate.out = projectBase + ".vcf"
    //annotate.jobOutputFile = annotate.out + ".out"
    //add(annotate)

    val targetEval = new VariantEval with CommandLineGATKArgs
    //targetEval.eval :+= annotate.out
    targetEval.eval :+= snpsIndelsVcf
    targetEval.dbsnp = qscript.pipeline.getProject.getEvalDbsnp
    targetEval.doNotUseAllStandardStratifications = true
    targetEval.doNotUseAllStandardModules = true
    targetEval.evalModule = List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
    targetEval.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "FunctionalClass", "Sample")
    targetEval.out = projectBase + ".eval"
    targetEval.jobOutputFile = targetEval.out + ".out"
    add(targetEval)

    if (qscript.expandIntervals > 0) {
      val flanksEval = new VariantEval with CommandLineGATKArgs
      //flanksEval.eval :+= annotate.out
      flanksEval.eval :+= snpsIndelsVcf
      flanksEval.dbsnp = qscript.pipeline.getProject.getEvalDbsnp
      flanksEval.intervals = List(flankIntervals)
      flanksEval.doNotUseAllStandardStratifications = true
      flanksEval.doNotUseAllStandardModules = true
      flanksEval.evalModule = List("SimpleMetricsByAC", "TiTvVariantEvaluator", "CountVariants")
      flanksEval.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "FunctionalClass", "Sample")
      flanksEval.out = projectBase + ".flanks.eval"
      flanksEval.jobOutputFile = flanksEval.out + ".out"
      add(flanksEval)
    }
  }
}
