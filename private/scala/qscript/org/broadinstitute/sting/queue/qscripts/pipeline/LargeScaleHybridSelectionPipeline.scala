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

package org.broadinstitute.sting.queue.qscripts.pipeline

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
            annotatedSNPs :+= annotateSNPs.out
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
