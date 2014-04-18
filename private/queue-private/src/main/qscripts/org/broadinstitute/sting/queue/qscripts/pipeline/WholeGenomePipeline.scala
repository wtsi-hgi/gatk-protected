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

import collection.immutable.ListMap
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
    //this.qSettings.dontRequestMultipleCores = true

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

    case class Interval(description: String, chr: String, start: Long, stop: Long) {
      def intervalString = "%s:%d-%d".format(chr, start, stop)
    }

    runType = runType.toLowerCase
    val contigSizes = IntervalUtils.getContigSizes(reference)
    var callIntervals = Seq.empty[Interval]
    var vqsrIntervals = Seq.empty[Interval]
    var evalIntervals = Seq.empty[Interval]
    if (runType == "wg") {
      val contigs = (1 to 22).map(_.toString) ++ Seq("X", "Y", "MT")
      callIntervals = contigs.map(chr => new Interval("chr" + chr, chr, 1, contigSizes.get(chr).longValue))
      vqsrIntervals = Nil
      evalIntervals = callIntervals.find(_.description == "chr20").toSeq
    } else {
      val locs = Seq(
        new Interval("cent1", "1", 121429168, 121529168),
        new Interval("cent16", "16", 40844464, 40944464),
        new Interval("chr20", "20", 1, contigSizes.get("20").longValue),
        new Interval("chr20_100k", "20", 100001, 200000))

      locs.find(_.description == runType) match {
        case Some(interval) =>
          callIntervals = Seq(interval)
          vqsrIntervals = Seq(interval)
          evalIntervals = Seq(interval)
        case None =>
          throw new RuntimeException("Invalid runType: " + runType + ". Must be one of: " + locs.map(_.description).mkString(", ") + ", or wg")
      }
    }

    require(bamList != null && bamList.getName.endsWith(".bam.list"), "-I/--bamList must be specified as <projectName>.bam.list")
    val project = bamList.getName.stripSuffix(".bam.list")

    var snpChrVcfs = ListMap.empty[Interval, File]
    var indelChrVcfs = ListMap.empty[Interval, File]
    val glModels = Seq(Model.SNP, Model.INDEL)

    for (interval <- callIntervals) {
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
          call.standard_min_confidence_threshold_for_calling = 10.0
          call.standard_min_confidence_threshold_for_emitting = 10.0
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
        combineChunks.intervalsString = Seq(interval.intervalString)
        combineChunks.variant = chunkVcfs.zipWithIndex.map { case (vcf, index) => TaggedFile(vcf, "input"+index) }
        combineChunks.rod_priority_list = chunkVcfs.zipWithIndex.map { case (vcf, index) => "input"+index }.mkString(",")
        combineChunks.suppressCommandLineHeader = true
        combineChunks.assumeIdenticalSamples = true
        combineChunks.out = chrDir + chrBase + "." + glModel.toString.toLowerCase + "s.unfiltered.vcf"
        add(combineChunks)

        glModel match {
          case Model.SNP =>
            snpChrVcfs += interval -> combineChunks.out
          case Model.INDEL =>
            indelChrVcfs += interval -> combineChunks.out
        }
      }
    }

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
      buildModel.intervalsString = vqsrIntervals.map(_.intervalString)
      buildModel.trustAllPolymorphic = true
      buildModel.TStranche = tranches
      buildModel.recal_file = project + "." + glModel.toString.toLowerCase + "s.recal.vcf"
      buildModel.tranches_file = project + "." + glModel.toString.toLowerCase + "s.tranches"
      add(buildModel)

      var tranche = 98.5
      var chrVcfs = ListMap.empty[Interval, File]
      glModel match {
        case Model.SNP =>
          chrVcfs = snpChrVcfs
          //tranche = ...
          buildModel.input = snpChrVcfs.values.toSeq
          buildModel.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
          buildModel.use_annotation = Seq("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "InbreedingCoeff", "DP")
          buildModel.resource :+= TaggedFile(hapmap, "training=true,truth=true,prior=15.0")
          buildModel.resource :+= TaggedFile(k1gOmni, "training=true,prior=12.0")
          buildModel.resource :+= TaggedFile(dbsnp135, "known=true,prior=6.0")
        case Model.INDEL =>
          chrVcfs = indelChrVcfs
          //tranche = ...
          buildModel.input = indelChrVcfs.values.toSeq
          buildModel.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
          buildModel.use_annotation = Seq("QD", "FS", "HaplotypeScore", "ReadPosRankSum", "InbreedingCoeff")
          buildModel.resource :+= TaggedFile(millsK1gIndels, "known=true,training=true,truth=true,prior=12.0")
      }

      for ((interval, chrVcf) <- chrVcfs) {
        val applyRecalibration = new ApplyRecalibration with CommandLineGATKArgs
        applyRecalibration.intervalsString = Seq(interval.intervalString)
        applyRecalibration.input :+= chrVcf
        applyRecalibration.mode = buildModel.mode
        applyRecalibration.recal_file = buildModel.recal_file
        applyRecalibration.tranches_file = buildModel.tranches_file
        applyRecalibration.ts_filter_level = tranche
        applyRecalibration.out = swapBaseExt(chrVcf, ".unfiltered.vcf", ".recalibrated.tranche_" + tranche.toString.replace(".", "_") + ".vcf")
        add(applyRecalibration)

        if (evalIntervals.contains(interval))
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
      eval.intervalsString = evalIntervals.map(_.intervalString)
      eval.dbsnp = dbsnp129
      eval.doNotUseAllStandardModules = true
      eval.evalModule = Seq("TiTvVariantEvaluator", "CountVariants", "CompOverlap")
      eval.doNotUseAllStandardStratifications = true
      eval.stratificationModule = Seq("EvalRod", "CompRod", "Novelty", "FunctionalClass") ++ strats
      eval.out = project + evalIntervals.map(_.description).mkString(".", "_", "") + strats.map(_.toLowerCase).mkString(".by_", "_", "") + ".eval"
      add(eval)
    }
  }

  private def swapBaseExt(file: File, oldExt: String, newExt: String): File =
    file.getParentFile + "/" + swapExt(file, oldExt, newExt)
}
