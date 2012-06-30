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
import org.broadinstitute.sting.queue.extensions.snpeff.SnpEff
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.queue.QScript
import collection.JavaConversions._

class HybridSelectionPipeline extends QScript {
  qscript =>

  private final val FILTER_VQSR = "VQSR"
  private final val FILTER_HARD = "HARD"

  @Input(doc="Tab separated squid projects and samples. Name must be <projectName>.tsv", shortName="tsv", exclusiveOf="bamList", required=false)
  var projectSampleTsv: File = _

  @Input(doc="BAM list files. Name must be <projectName>.bam.list", shortName="I", exclusiveOf="projectSampleTsv", required=false)
  var bamList: File = _

  @Input(doc="GATK or Picard intervals file.", shortName="L", exclusiveOf="projectSampleTsv", required=false)
  var intervals: File = _

  @Input(doc="Level of parallelism for UnifiedGenotyper. By default set to 20.", shortName="varScatter", required=false)
  var variantCallerScatterCount = 20

  @Argument(doc="Pipeline memory limit. By default set to 2g.", shortName="pipeMemory", required=false)
  var pipelineMemoryLimit = 2

  @Argument(doc="Variant filter type, " + FILTER_VQSR + " or " + FILTER_HARD + ". By default hard filters are used unless the bait set is whole_exome_agilent_1.", shortName="varFilter", required=false)
  var variantFilterType: String = _

  @Argument(doc="Memory limit for snpEff. By default set to 4g.", shortName="snpEffMemory", required=false)
  var snpEffMemory = 4

  @Argument(doc="Expand each target in input intervals by the specified number of bases. By default set to 50 bases.", shortName="expand", required=false)
  var expandIntervals = 50

  def script() {
    val exomeIntervals = new File("/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list")
    val resources = "/humgen/gsa-pipeline/resources/b37/v4/"
    val k1gExomesBam = resources + "1000_Genomes_Whole_Exome_50_Samples.bam"
    val k1gExomesSamples = resources + "1000_Genomes_Whole_Exome_50_Samples.samples"
    val k1gTrainingHighQuality = resources + "phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"
    val k1gTrainingTerrible = resources + "phase1.wgs.projectConsensus.v2b.recal.terrible.vcf"
    val omni = resources + "1000G_omni2.5.b37.sites.vcf"
    val hapmap = resources + "hapmap_3.3.b37.sites.vcf"
    val dbsnp129 = resources + "dbsnp_135.b37.excluding_sites_after_129.vcf"
    val dbsnp135 = resources + "dbsnp_135.b37.vcf"
    val goldStandardIndels = resources + "Mills_and_1000G_gold_standard.indels.b37.sites.vcf"

    var reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
    var useK1gExomes = false

    var projectName: String = null

    if (projectSampleTsv == null) {
      require(bamList != null && bamList.getName.endsWith(".bam.list"), "-I/--bamList must be specified as <projectName>.bam.list")
      require(intervals != null, "-L/--intervals must be specified")
      projectName = bamList.getName.stripSuffix(".bam.list")
    } else {
      projectName = FilenameUtils.getBaseName(projectSampleTsv)

      val picardSamples = PicardAggregationUtils.parseSamples(projectSampleTsv)
      val picardIntervals = PicardAggregationUtils.readAnalysisIntervals(picardSamples)

      val writeBamList = new ListWriterFunction
      writeBamList.inputFiles = PicardAggregationUtils.getSampleBams(picardSamples).toSeq
      writeBamList.listFile = projectName +".bam.list"
      add(writeBamList)

      // Use the Picard reference instead of the 1000 Genomes version
      reference = picardIntervals.getReference
      intervals = picardIntervals.getTargets
      bamList = writeBamList.listFile

      if (intervals == exomeIntervals && writeBamList.inputFiles.size < 50) {
        // Allow the HSPTest to explicitly use hard filters. Adding 50 other 1KG exomes increases test runtime.
        useK1gExomes = (variantFilterType == null || variantFilterType == FILTER_VQSR)
      }
    }

    qSettings.runName = projectName

    variantFilterType = {
      if (variantFilterType != null)
        variantFilterType.toUpperCase
      else if (intervals == exomeIntervals)
        FILTER_VQSR
      else
        FILTER_HARD
    }

    trait CommandLineGATKArgs extends CommandLineGATK {
      this.reference_sequence = reference
      this.intervals = Seq(qscript.intervals)
      this.memoryLimit = pipelineMemoryLimit
    }

    trait ExpandedIntervals extends CommandLineGATK {
      if (qscript.expandIntervals > 0)
        this.interval_padding = qscript.expandIntervals
    }

    val call = new UnifiedGenotyper with CommandLineGATKArgs with ExpandedIntervals
    call.input_file = Seq(bamList)
    if (useK1gExomes)
      call.input_file :+= k1gExomesBam
    call.dbsnp = dbsnp135
    call.downsample_to_coverage = 600
    call.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    call.out = projectName + ".unfiltered.vcf"
    call.scatterCount = qscript.variantCallerScatterCount
    add(call)

    val selectSNPs = new SelectVariants with CommandLineGATKArgs with ExpandedIntervals
    selectSNPs.selectType :+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.SNP
    selectSNPs.variant = call.out
    selectSNPs.out = projectName + ".snps.unfiltered.vcf"
    add(selectSNPs)

    val selectIndels = new SelectVariants with CommandLineGATKArgs with ExpandedIntervals
    selectIndels.selectType :+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.INDEL
    selectIndels.variant = call.out
    selectIndels.out = projectName + ".indels.unfiltered.vcf"
    add(selectIndels)

    val filteredSNPsVcf = variantFilterType match {
      case FILTER_VQSR =>
        val trancheLevels = Seq(
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
        buildSNPModel.resource :+= TaggedFile(hapmap, "training=true,truth=true,prior=15.0")
        buildSNPModel.resource :+= TaggedFile(omni, "training=true,truth=true,prior=12.0")
        buildSNPModel.resource :+= TaggedFile(k1gTrainingHighQuality, "training=true,prior=12.0")
        buildSNPModel.resource :+= TaggedFile(k1gTrainingTerrible, "bad=true,prior=2.0")
        buildSNPModel.resource :+= TaggedFile(dbsnp135, "known=true,prior=4.0")
        buildSNPModel.use_annotation = Seq("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "InbreedingCoeff")
        buildSNPModel.trustAllPolymorphic = true
        buildSNPModel.maxGaussians = 6
        buildSNPModel.TStranche = trancheLevels
        buildSNPModel.tranches_file = projectName + ".snps.tranches"
        buildSNPModel.recal_file = projectName + ".snps.recal"
        add(buildSNPModel)

        val applySNPModel = new ApplyRecalibration with CommandLineGATKArgs with ExpandedIntervals
        applySNPModel.input :+= selectSNPs.out
        applySNPModel.tranches_file = buildSNPModel.tranches_file
        applySNPModel.recal_file = buildSNPModel.recal_file
        applySNPModel.ts_filter_level = 98.5
        applySNPModel.out = projectName + ".snps.recalibrated.vcf"
        add(applySNPModel)

        applySNPModel.out

      case FILTER_HARD =>
        val filterSNPs = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
        filterSNPs.variant = selectSNPs.out
        filterSNPs.filterName = Seq("SNP_QD", "SNP_MQ", "SNP_HaplotypeScore", "SNP_MQRankSum", "SNP_ReadPosRankSum", "SNP_FS")
        filterSNPs.filterExpression = Seq("QD<2.0", "MQ<40.0", "HaplotypeScore>13.0", "MQRankSum<-12.5", "ReadPosRankSum<-8.0", "FS>60.0")
        filterSNPs.out = projectName + ".snps.filtered.vcf"
        add(filterSNPs)

        filterSNPs.out

      case unknown =>
        throw new IllegalArgumentException("""If set the variantFilterType must be "%s" or "%s"""".format(FILTER_VQSR, FILTER_HARD))
    }

    val filterIndels = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
    filterIndels.variant = selectIndels.out
    filterIndels.filterName = Seq("Indel_FS", "Indel_QD", "Indel_ReadPosRankSum", "Indel_InbreedingCoeff")
    filterIndels.filterExpression = Seq("FS>200.0", "QD<2.0", "ReadPosRankSum<-20.0", "InbreedingCoeff<-0.8")
    filterIndels.out = projectName + ".indels.filtered.vcf"
    add(filterIndels)

    val combineSNPsIndels = new CombineVariants with CommandLineGATKArgs with ExpandedIntervals
    combineSNPsIndels.variant :+= TaggedFile(filterIndels.out, "indels")
    combineSNPsIndels.variant :+= TaggedFile(filteredSNPsVcf, "snps")
    combineSNPsIndels.filteredrecordsmergetype = org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    combineSNPsIndels.assumeIdenticalSamples = true
    combineSNPsIndels.out = projectName + ".unannotated.vcf"
    add(combineSNPsIndels)

    val snpsIndelsVcf =
      if (!useK1gExomes) {
        combineSNPsIndels.out
      } else {
        val removeK1gSamples = new SelectVariants with CommandLineGATKArgs with ExpandedIntervals
        removeK1gSamples.variant = combineSNPsIndels.out
        removeK1gSamples.exclude_sample_file :+= k1gExomesSamples
        removeK1gSamples.excludeNonVariants = true
        removeK1gSamples.out = projectName + ".selected.vcf"
        add(removeK1gSamples)

        removeK1gSamples.out
      }

    val snpEff = new SnpEff
    snpEff.inVcf = snpsIndelsVcf
    snpEff.config = "/humgen/gsa-pipeline/resources/snpEff/v2_0_5/snpEff.config"
    snpEff.genomeVersion = "GRCh37.64"
    snpEff.onlyCoding = true
    snpEff.outVcf = projectName + ".snpeff.vcf"
    snpEff.memoryLimit = snpEffMemory
    add(snpEff)

    val annotate = new VariantAnnotator with CommandLineGATKArgs with ExpandedIntervals
    annotate.variant = snpsIndelsVcf
    annotate.snpEffFile = snpEff.outVcf
    annotate.annotation :+= "SnpEff"
    annotate.out = projectName + ".vcf"
    add(annotate)

    def newVariantEval = {
      val eval = new VariantEval with CommandLineGATKArgs
      eval.eval :+= annotate.out
      eval.dbsnp = dbsnp129
      eval.goldStandard = goldStandardIndels
      eval.doNotUseAllStandardModules = true
      eval.doNotUseAllStandardStratifications = true
      eval
    }

    for (strats <- Seq(
      Seq("AlleleCount"),
      Seq("Sample")
    )) {
      val stratsEval = newVariantEval
      stratsEval.evalModule = Seq("TiTvVariantEvaluator", "CountVariants", "CompOverlap", "IndelSummary", "MultiallelicSummary")
      stratsEval.stratificationModule = Seq("EvalRod", "CompRod", "Novelty", "FunctionalClass") ++ strats
      stratsEval.out = projectName + strats.map(_.toLowerCase).mkString(".by_", "_", ".eval")
      add(stratsEval)
    }

    val indelEval = newVariantEval
    indelEval.evalModule = List("IndelSummary", "IndelLengthHistogram")
    indelEval.stratificationModule = Seq("EvalRod", "CompRod", "Sample", "TandemRepeat", "OneBPIndel")
    indelEval.out = projectName + ".indel.eval"
    add(indelEval)
  }
}
