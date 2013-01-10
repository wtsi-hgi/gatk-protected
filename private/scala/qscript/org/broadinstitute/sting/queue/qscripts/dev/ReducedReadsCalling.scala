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

package org.broadinstitute.sting.queue.qscripts.dev

import org.apache.commons.io.FilenameUtils
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.snpeff.SnpEff
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.queue.QScript
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
      // TODO -- fix Queue so that ReadFilters work correctly
      // val reduce = new ReduceReads() with CommandLineGATKArgs with ExpandedIntervals with BadMate
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

    // TODO -- 600 is too high when there are many samples; find a better value
    // call.downsample_to_coverage = 600
    call.downsample_to_coverage = 60

    val CALL_TOGETHER = false
    var snpsFile: String = null
    var indelsFile: String = null

    // to call both at once
    if ( CALL_TOGETHER ) {
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
      snpsFile = selectSNPs.out
      add(selectSNPs)

      val selectIndels = new SelectVariants with CommandLineGATKArgs with ExpandedIntervals
      selectIndels.selectType :+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.INDEL
      selectIndels.variant = call.out
      selectIndels.out = projectName + ".indels.unfiltered.vcf"
      selectIndels.jobOutputFile = selectIndels.out + ".out"
      indelsFile = selectIndels.out
      add(selectIndels)
    } else {
      call.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
      call.out = projectName + ".snps.unfiltered.vcf"
      call.jobOutputFile = call.out + ".out"
      call.scatterCount = qscript.variantCallerScatterCount
      snpsFile = call.out
      add(call)

      val callIndels = new UnifiedGenotyper with CommandLineGATKArgs with ExpandedIntervals
      callIndels.input_file = call.input_file
      callIndels.dbsnp = call.dbsnp
      callIndels.downsample_to_coverage = call.downsample_to_coverage
      callIndels.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
      callIndels.max_alternate_alleles = 2
      callIndels.out = projectName + ".indels.unfiltered.vcf"
      callIndels.jobOutputFile = callIndels.out + ".out"
      callIndels.scatterCount = qscript.variantCallerScatterCount
      indelsFile = callIndels.out
      add(callIndels)
    }


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
    buildSNPModel.input :+= snpsFile
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
    applySNPModel.input :+= snpsFile
    applySNPModel.tranches_file = buildSNPModel.tranches_file
    applySNPModel.recal_file = buildSNPModel.recal_file
    applySNPModel.ts_filter_level = 98.5
    val filteredSNPsVcf = projectName + ".snps.recalibrated.vcf"
    applySNPModel.out = filteredSNPsVcf
    applySNPModel.jobOutputFile = applySNPModel.out + ".out"
    add(applySNPModel)

    val filterIndels = new VariantFiltration with CommandLineGATKArgs with ExpandedIntervals
    filterIndels.variant = indelsFile
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

    val eval = new VariantEval with CommandLineGATKArgs with ExpandedIntervals
    eval.eval :+= combineSNPsIndels.out
    eval.dbsnp = dbsnp129
    eval.doNotUseAllStandardModules = true
    eval.evalModule = Seq("TiTvVariantEvaluator", "CountVariants", "CompOverlap", "MultiallelicSummary")
    eval.doNotUseAllStandardStratifications = true
    eval.stratificationModule = Seq("EvalRod", "CompRod", "Novelty")
    eval.out = projectName + ".eval"
    eval.jobOutputFile = eval.out + ".out"
    add(eval)

    val snpEff = new SnpEff
    snpEff.inVcf = combineSNPsIndels.out
    snpEff.config = "/humgen/gsa-pipeline/resources/snpEff/v2_0_5/snpEff.config"
    snpEff.genomeVersion = "GRCh37.64"
    snpEff.onlyCoding = true
    snpEff.outVcf = projectName + ".snpeff.vcf"
    snpEff.jobOutputFile = snpEff.outVcf + ".out"
    snpEff.memoryLimit = 4
    add(snpEff)

    val annotate = new VariantAnnotator with CommandLineGATKArgs with ExpandedIntervals
    annotate.variant = combineSNPsIndels.out
    annotate.snpEffFile = snpEff.outVcf
    annotate.annotation :+= "SnpEff"
    annotate.out = projectName + ".annotated.vcf"
    annotate.jobOutputFile = annotate.out + ".out"
    add(annotate)
  }
}
