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

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.snpeff.SnpEff
import org.broadinstitute.sting.queue.function._
import org.broadinstitute.sting.queue.QScript

class LargeScaleHybridSelectionPipeline extends QScript {
  qscript =>

  @Input(doc="BAM list files. Name must be <projectName>.bam.list", shortName="I")
  var bamList: File = _

  @Argument(doc="chrs", shortName="C", fullName = "chr")
  var chrs: Seq[String] = Nil

  def script() {
    val mergeBamMemoryLimit = 16
    val callingMemoryLimit = 8
    val pipelineMemoryLimit = 4
    val VQSRMemoryLimit = 32
    val filterCalls = true
    val outputFormat = "bcf"
    val mergeBamDirRoot = "/humgen/gsa-firehose/ReduceReads_v2fixed_JointCalling/MacArthur_Chr1_RRv2fixed_JointCalling/"
    val expandIntervals = 50
    val variantCallerScatterCount = 1000

    val reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
    val resources = "/humgen/gsa-pipeline/resources/b37/v5/"
    val k1gTrainingHighQuality = resources + "phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"
    val k1gTrainingTerrible = resources + "phase1.wgs.projectConsensus.v2b.recal.terrible.vcf"
    val omni = resources + "1000G_omni2.5.b37.vcf"
    val hapmap = resources + "hapmap_3.3.b37.vcf"
    val dbsnp129 = resources + "dbsnp_137.b37.excluding_sites_after_129.vcf"
    val dbsnp137 = resources + "dbsnp_137.b37.vcf"
    val goldStandardIndels = resources + "Mills_and_1000G_gold_standard.indels.b37.vcf"
    val intervalsFile = resources + "gencode.v12_broad.agilent_merged.interval_list"
    val standardBroadIntervals = resources + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"

    require(bamList != null && bamList.getName.endsWith(".reduced.bam.list"), "-I/--bamList must be specified as <projectName>.reduced.bam.list")
    val projectName = bamList.getName.stripSuffix(".reduced.bam.list")

    qSettings.runName = projectName

    trait CommandLineGATKArgs extends CommandLineGATK with RetryMemoryLimit {
      this.reference_sequence = reference
      this.intervals = Seq(intervalsFile)
      this.interval_set_rule = org.broadinstitute.sting.utils.interval.IntervalSetRule.INTERSECTION
      this.memoryLimit = pipelineMemoryLimit
    }

    trait ExpandedIntervals extends CommandLineGATK {
      this.interval_padding = expandIntervals
    }

    var unfilteredSNPs = Seq.empty[File]
    var unfilteredIndels = Seq.empty[File]
    var recalibratedSNPs = Seq.empty[File]
    var recalibratedIndels = Seq.empty[File]
    var evalInputs = Seq.empty[File]
    var annotatedIndels = Seq.empty[File]
    var annotatedSNPs = Seq.empty[File]

    val unannotatedSitesOnly = projectName + ".unannotated.sites_only.vcf"
    val snpEffOutput = projectName + ".snpeff.vcf"
    val tranchesSNPs = projectName + ".snps.tranches"
    val recalSNPs = projectName + ".snps.recal"
    val tranchesIndels = projectName + ".indels.tranches"
    val recalIndels = projectName + ".indels.recal"

    val bams = io.Source.fromFile(bamList).getLines().toSeq
    val bamGroups = bams.grouped(1000).toSeq

    val mergedBamDir = mergeBamDirRoot +
      "%s/mergedBams/".format(projectName)

    var bamGroupNumber = 0
    var mergeBamLists = Seq.empty[File]

    val annotate = false

    for (bamGroup <- bamGroups) {
      val mergeBamList = new ListWriterFunction
      mergeBamList.inputFiles = bamGroup
      mergeBamList.listFile = mergedBamDir + "bamLists/%03d.bam.list".format(bamGroupNumber)
      add(mergeBamList)

      mergeBamLists :+= mergeBamList.listFile
      bamGroupNumber += 1
    }

    for (chr <- chrs) {

      val chrDir = "chrs/" + chr + "/"
      val chrBase = projectName + ".chr" + chr

      trait ChromosomeIntervals extends CommandLineGATKArgs with ExpandedIntervals {
        this.intervalsString :+= chr
      }

      bamGroupNumber = 0
      var chrBams = Seq.empty[File]
      for (bamGroup <- bamGroups) {
        val mergeBam = new PrintReads with ChromosomeIntervals
        mergeBam.input_file = Seq(mergeBamLists(bamGroupNumber))
        mergeBam.out = mergedBamDir + chrDir + "chr%s.%03d.bam".format(chr, bamGroupNumber)
	      mergeBam.memoryLimit = mergeBamMemoryLimit
        mergeBam.nct = 4 //todo is it make sense to use it with such an high memoryLimit?
        add(mergeBam)

        chrBams :+= mergeBam.out
        bamGroupNumber += 1
      }


      val chrMergeBamList = new ListWriterFunction
      chrMergeBamList.inputFiles = chrBams
      chrMergeBamList.listFile = mergedBamDir + chrDir + "chr%s.merged.bam.list".format(chr)
      add(chrMergeBamList)


      val calling = true
      if(calling){
        val call = new UnifiedGenotyper with ChromosomeIntervals
        call.input_file = Seq(chrMergeBamList.listFile)
        call.dbsnp = dbsnp137
        call.downsample_to_coverage = 60
        call.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
        call.out = chrDir + chrBase + ".unfiltered." + outputFormat
        call.nct = 4
        call.maxRuntimeUnits = java.util.concurrent.TimeUnit.HOURS
        call.maxRuntime = 24
        call.memoryLimit = callingMemoryLimit
        call.scatterCount = variantCallerScatterCount
        call.jobName = chrBase + ".unfiltered." + outputFormat

        add(call)


        val selectSNPs = new SelectVariants with ChromosomeIntervals
        selectSNPs.selectType :+= org.broadinstitute.variant.variantcontext.VariantContext.Type.SNP
        selectSNPs.variant = call.out
        selectSNPs.out = chrDir + chrBase  + ".snps.unfiltered." + outputFormat
        selectSNPs.nt = 16
        add(selectSNPs)

        val selectIndels = new SelectVariants with ChromosomeIntervals
        selectIndels.selectType :+= org.broadinstitute.variant.variantcontext.VariantContext.Type.INDEL
        selectIndels.variant = call.out
        selectIndels.out = chrDir + chrBase + ".indels.unfiltered." + outputFormat
        selectSNPs.nt = 16
        add(selectIndels)

        if (filterCalls) {
          val applySNPModel = new ApplyRecalibration with ChromosomeIntervals
          applySNPModel.input :+= selectSNPs.out
          applySNPModel.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
          applySNPModel.tranches_file = tranchesSNPs
          applySNPModel.recal_file = recalSNPs
          applySNPModel.ts_filter_level = 98.5
          applySNPModel.out = chrDir + chrBase + ".snps.recalibrated." + outputFormat
          applySNPModel.nt = 16
          add(applySNPModel)

          val applyIndelModel = new ApplyRecalibration with ChromosomeIntervals
          applyIndelModel.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
          applyIndelModel.input :+= selectIndels.out
          applyIndelModel.tranches_file = tranchesIndels
          applyIndelModel.recal_file = recalIndels
          applyIndelModel.ts_filter_level = 98.5
          applyIndelModel.out = chrDir + chrBase + ".indels.recalibrated." + outputFormat
          applyIndelModel.nt = 16
          add(applyIndelModel)


          if (annotate){
            val annotateSNPs = new VariantAnnotator with ChromosomeIntervals
            annotateSNPs.variant = applySNPModel.out
            annotateSNPs.snpEffFile = snpEffOutput
            annotateSNPs.annotation :+= "SnpEff"
            annotateSNPs.out = chrDir + chrBase + ".snps." + outputFormat
            add(annotateSNPs)

            val annotateIndels = new VariantAnnotator with ChromosomeIntervals
            annotateIndels.variant = applyIndelModel.out
            annotateIndels.snpEffFile = snpEffOutput
            annotateIndels.annotation :+= "SnpEff"
            annotateIndels.out = chrDir + chrBase + ".indels." + outputFormat
            add(annotateIndels)

            if (outputFormat != "vcf") {
              val bcfToVcfSNPs = new SelectVariants with ChromosomeIntervals
              bcfToVcfSNPs.variant = annotateSNPs.out
              bcfToVcfSNPs.out = chrDir + chrBase + ".SNPs.vcf"
              add(bcfToVcfSNPs)

              val bcfToVcfIndels = new SelectVariants with ChromosomeIntervals
              bcfToVcfIndels.variant = annotateIndels.out
              bcfToVcfIndels.out = chrDir + chrBase + ".indels.vcf"
              add(bcfToVcfIndels)
            }
            evalInputs = Seq(annotateSNPs.out, annotateIndels.out)
            annotatedSNPs :+= annotateIndels.out
            annotatedIndels :+= annotateIndels.out
          }
          else{
            if (outputFormat != "vcf") {
              val bcfToVcfSNPs = new SelectVariants with ChromosomeIntervals
              bcfToVcfSNPs.variant = applySNPModel.out
              bcfToVcfSNPs.out = chrDir + chrBase + ".SNPs.recalibrated.vcf"
              add(bcfToVcfSNPs)

              val bcfToVcfIndels = new SelectVariants with ChromosomeIntervals
              bcfToVcfIndels.variant = applyIndelModel.out
              bcfToVcfIndels.out = chrDir + chrBase + ".indels.recalibrated.vcf"
              add(bcfToVcfIndels)
            }
          }

          unfilteredSNPs :+= selectSNPs.out
          unfilteredIndels :+= selectIndels.out
          recalibratedSNPs :+= applySNPModel.out
          recalibratedIndels :+= applyIndelModel.out

          if (chr == "20") {
            if (!annotate)
              evalInputs = Seq(applySNPModel.out, applyIndelModel.out)

            def newVariantEval = {
              val eval = new VariantEval with CommandLineGATKArgs
              eval.intervals = Seq(standardBroadIntervals)
              eval.eval = evalInputs
              eval.mergeEvals = true
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
              stratsEval.out = chrBase + strats.map(_.toLowerCase).mkString(".by_", "_", ".eval")
              stratsEval.nt = 8
              add(stratsEval)
            }

            val indelEval = newVariantEval
            indelEval.evalModule = List("IndelSummary", "IndelLengthHistogram")
            indelEval.stratificationModule = Seq("EvalRod", "CompRod", "Sample", "TandemRepeat", "OneBPIndel")
            indelEval.out = chrBase + ".indel.eval"
            indelEval.nt = 8
            add(indelEval)


          }


        }
      }
    }

    if (filterCalls) {

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
      buildSNPModel.input = unfilteredSNPs
      buildSNPModel.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
      buildSNPModel.resource :+= TaggedFile(hapmap, "training=true,truth=true,prior=15.0")
      buildSNPModel.resource :+= TaggedFile(omni, "training=true,truth=true,prior=12.0")
      buildSNPModel.resource :+= TaggedFile(k1gTrainingHighQuality, "training=true,prior=12.0")
      buildSNPModel.resource :+= TaggedFile(k1gTrainingTerrible, "bad=true,prior=2.0")
      buildSNPModel.resource :+= TaggedFile(dbsnp137, "known=true,prior=4.0")
      buildSNPModel.use_annotation = Seq("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "InbreedingCoeff")
      buildSNPModel.trustAllPolymorphic = true
      buildSNPModel.maxGaussians = 6
      buildSNPModel.TStranche = trancheLevels
      buildSNPModel.tranches_file = projectName + ".snps.tranches"
      buildSNPModel.recal_file = projectName + ".snps.recal"
      buildSNPModel.nt = 16
      buildSNPModel.memoryLimit = VQSRMemoryLimit
      add(buildSNPModel)

      val buildIndelModel = new VariantRecalibrator with CommandLineGATKArgs with ExpandedIntervals
      buildIndelModel.input = unfilteredIndels
      buildIndelModel.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
      buildIndelModel.resource :+= TaggedFile(goldStandardIndels, "known=true,training=true,truth=true,prior=12.0")
      buildIndelModel.use_annotation = Seq("QD", "HaplotypeScore", "ReadPosRankSum", "FS", "InbreedingCoeff")
      buildIndelModel.trustAllPolymorphic = true
      buildIndelModel.maxGaussians = 4
      buildIndelModel.stdThreshold = 10
      buildIndelModel.percentBadVariants = 0.12
      buildIndelModel.TStranche = trancheLevels
      buildIndelModel.tranches_file = projectName + ".indels.tranches"
      buildIndelModel.recal_file = projectName + ".indels.recal"
      buildIndelModel.nt = 16
      buildIndelModel.memoryLimit = VQSRMemoryLimit
      add(buildIndelModel)


      if (annotate){
        val combineSNPsIndels = new CombineVariants with CommandLineGATKArgs with ExpandedIntervals
        combineSNPsIndels.variant ++= recalibratedSNPs
        combineSNPsIndels.variant ++= recalibratedIndels
        combineSNPsIndels.sites_only = true
        combineSNPsIndels.filteredrecordsmergetype = org.broadinstitute.variant.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
        combineSNPsIndels.assumeIdenticalSamples = true
        combineSNPsIndels.out = unannotatedSitesOnly
        add(combineSNPsIndels)

        // Create a sites only VCF for snpEff

        val snpEff = new SnpEff
        snpEff.inVcf = combineSNPsIndels.out   //same as unannotatedSitesOnly
        snpEff.config = "/humgen/gsa-pipeline/resources/snpEff/v2_0_5/snpEff.config"
        snpEff.genomeVersion = "GRCh37.64"
        snpEff.onlyCoding = true
        snpEff.outVcf = snpEffOutput
        snpEff.memoryLimit = pipelineMemoryLimit
        add(snpEff)

        evalInputs = annotatedIndels ++ annotatedSNPs
      }

      if (!annotate)
        evalInputs = recalibratedSNPs ++ recalibratedIndels

      def newVariantEval = {
        val eval = new VariantEval with CommandLineGATKArgs
        eval.intervals = Seq(standardBroadIntervals)
        eval.eval = evalInputs
        eval.mergeEvals = true
        eval.dbsnp = dbsnp129
        eval.goldStandard = goldStandardIndels
        eval.doNotUseAllStandardModules = true
        eval.doNotUseAllStandardStratifications = true
        eval.memoryLimit = 32
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
        stratsEval.nt = 8
        add(stratsEval)
      }

      val indelEval = newVariantEval
      indelEval.evalModule = List("IndelSummary", "IndelLengthHistogram")
      indelEval.stratificationModule = Seq("EvalRod", "CompRod", "Sample", "TandemRepeat", "OneBPIndel")
      indelEval.out = projectName + ".indelQC.eval"
      indelEval.nt = 8
      add(indelEval)





    }
  }
}
