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

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport
import org.broadinstitute.sting.queue.extensions.gatk.BaseRecalibrator
import org.broadinstitute.sting.queue.extensions.gatk.PrintReads
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.function._
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils.FilteredRecordMergeType
import org.broadinstitute.variant.variantcontext.VariantContext
import org.broadinstitute.sting.commandline.ClassType

class GeneralCallingPipeline extends QScript {
  qscript =>

  // Required arguments.  All initialized to empty values.

  @Input(shortName="R", doc="The reference file for the bam files.")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(shortName="I", doc="Bam file or .list of BAM files to genotype.")
  var bamFile: File = _

  // The following arguments are all optional.

  @ClassType(classOf[Integer])
  @Argument(shortName="V", doc="The pipeline version", required=false)
  var pipelineVersions: List[Int] = List(2,1)

  @Argument(shortName = "RPath", doc="RPath", required=false)
  var RPath: File = new File("../R")

  @Argument(shortName="sc", doc="set the number of the scattered jobs", required=false)
  var scatterCount: Int = 0

  @Argument(shortName="outputDir", doc="output directory", required=false)
  var outputDir: String = "./tmp"

  @Argument(shortName="L", doc="An optional file with a list of intervals to proccess.",  required=false)
  var intervals: String = null

  @Argument(shortName="filter", doc="A optional list of filter names.", required=false)
  var filterNames: List[String] = Nil // Nil is an empty List, versus null which means a non-existent List.

  @Argument(shortName="noBAQ", doc="turns off BAQ calculation", required=false)
  var noBAQ: Boolean = false

  @Argument(shortName="filterExpression", doc="An optional list of filter expressions.", required=false)
  var filterExpressions: List[String] = Nil

  @Argument(shortName="samples", doc="Samples to include in Variant Eval", required=false)
  var samples: List[String] = Nil

  @Argument(shortName="deletions", doc="Maximum deletion fraction allowed at a site to call a genotype.", required=false)
  var deletions: Double = -1

  @Input(shortName = "XL", doc="Exclude intervals list", required=false)
  var excludeIntervals: List[File] = Nil

  @Input(shortName="ped", fullName="pedigree", doc="PED file of the family", required=false)
  var ped: File = "/broad/hptmp/ami/tmp/CEUTrio.ped"

  @Argument(shortName="bamDir", doc="Subdirectory to store the reduced bams. By default set to 'reduced'.", required=false)
  var mergeBamDir = "reducedBAMs/"

  @Argument(shortName="rrMem", doc="Reduce reads memory limit.",  required=false)
  var reduceReadsMemoryLimit = 4

  @Argument(shortName="callingMem", doc="calling (UG/HC) memory limit.",  required=false)
  var callingMemoryLimit = 6

  /************* include/exclude steps of the pipeline ***********************/

  @Argument(shortName="useBQSR.2.0", doc="turn on a first step of using 2.0 BQSR on the input bam file", required=false)
  var useBQSR2: Boolean = false

  @Argument(shortName="PBT", doc="apply phaseByTrasmission on the output vcf", required=false)
  var usePhaseByTransmission: Boolean = false

  @Argument(shortName="RBP", doc="apply readBackPhasing on the output vcf", required=false)
  var useReadBackPhasing: Boolean = false

  @Argument(shortName="summaryReport", doc="create a pdf file with eval summery tables of the output vcf", required=false)
  var createEvalSummaryReport: Boolean = false

  @Argument(shortName="fullReport", doc="create a pdf file with full QC eval report of the output vcf", required=false)
  var createFullPostQCReport: Boolean = false

  @Argument(shortName="RR", doc="produces reduce reads bam files", required=false)
  var CallReduceReads: Boolean = false

  @Argument(shortName="excludeBadRegions", doc="If provided, we'll skip the really slow bits of the genome", required=false)
  var excludeBadRegions: Boolean = false


  val latestdbSNP = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_137.b37.vcf"  // Best Practices v4
  val hapmapSites = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/hapmap_3.3.b37.vcf"                       // Best Practices v4
  val hapmapGenotypes = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf"                       // Best Practices v4
  val omni_b37_sites = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_2141_samples.b37.vcf"    // Best Practices v4
  val omni_b37_genotypes = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_genotypes_2141_samples.b37.vcf"    // Best Practices v4
  val training_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"  // from the MethodDevelopmentCallingPipeline scala script
  val indelGoldStandardCallset  = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/Mills_and_1000G_gold_standard.indels.b37.vcf" // Best Practices v4
  val dbSNP_129 = "/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.vcf"
  val omni_mono: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_monomorphic_2141_samples.b37.vcf"
  val MVL_FP: String =  "/humgen/gsa-hpprojects/GATK/badSitesDB/Autism.HighMVLR.Exome.indels.vcf"
  val CEUTrio_ped = "/broad/hptmp/ami/tmp/CEUTrio.ped" //todo: move to a proper location
  // TODO -- this really needs to be fixed!
  val BQSRFILE: String = "/humgen/gsa-firehose2/carneiro/agbt/chr20/NA12878-2x250.bwasw.chr20.grp"

  val queueLogDir = ".qlog/"
  val noET_key = "/humgen/gsa-hpprojects/GATK/data/gatk_user_keys/gsamembers_broadinstitute.org.key"

  trait BaseCommandArguments extends CommandLineGATK with RetryMemoryLimit {
    this.logging_level = "INFO"
    phone_home = GATKRunReport.PhoneHomeOption.NO_ET
    gatk_key = noET_key
    this.reference_sequence = qscript.referenceFile
    if ( qscript.intervals != null )
      this.intervalsString :+= qscript.intervals
    else
      this.intervalsString = List("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
        "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")
    this.excludeIntervals = qscript.excludeIntervals
    this.memoryLimit = 2
    if ( excludeBadRegions )
      this.excludeIntervalsString = List("20:25,730,140-26,319,926", "20:29,419,342-29,655,147")
  }

  // 0) BQSR v 2.0
  class BQSR( bam: File, recal: File ) extends BaseRecalibrator with BaseCommandArguments {

    this.input_file :+= bam
    this.out = recal
    this.knownSites ++= List(new File(indelGoldStandardCallset))
    this.knownSites ++= List(new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf"))
    this.knownSites ++= List(new File("/humgen/gsa-hpprojects/dev/carneiro/bqsr/data/projectConsensus.snps.vcf"))
    this.knownSites ++= List(new File(latestdbSNP))
    this.memoryLimit = 8
    this.qq = 0
    this.mcs = 2
    this.ics = 3
    //this.np = true
    this.javaGCThreads = 4
    this.scatterCount = 200 //if(qscript.scatterCount == 0) 200 else qscript.scatterCount
  }

  class RecalBAM( bam: File, recal: File ) extends PrintReads with BaseCommandArguments {

    this.input_file :+= bam
    this.out = qscript.outputDir + "/" + recal +".out"
    this.BQSR = recal
    this.memoryLimit = 6
    this.qq = 0
    this.javaGCThreads = 4
    this.scatterCount = 200 //if(qscript.scatterCount == 0) 200 else qscript.scatterCount
  }

  // 1a.) Call SNPs with UG
  class MyUnifiedGenotyper(name: String, inputBamFile: File) extends UnifiedGenotyper with BaseCommandArguments {
    this.input_file :+= inputBamFile
    if (qscript.deletions >= 0)
      this.max_deletion_fraction = qscript.deletions
    this.out = qscript.outputDir + "/" + name + ".UnifiedGenotyper.unfiltered.vcf"
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    this.analysisName = "CEUTrio_UGs"
    this.jobName =  queueLogDir + "CEU_Trio.ug"
    this.scatterCount = if(qscript.scatterCount == 0) 80 else qscript.scatterCount
    this.stand_call_conf = 30.0
    this.stand_emit_conf = 30.0
    this.D = new File(latestdbSNP)
    this.memoryLimit = qscript.callingMemoryLimit
  }

  // 1.) HaplotypeCaller
  class MyHaplotypeCaller(name:String, inputHCFile: File) extends HaplotypeCaller with BaseCommandArguments {
    this.scatterCount = if(qscript.scatterCount == 0) 200 else qscript.scatterCount
    this.input_file :+= inputHCFile
    this.out = qscript.outputDir + "/" + name + ".HaplotypeCaller.unfiltered.vcf"
    this.analysisName = "HaplotypeCaller"
    this.javaGCThreads = 4
    this.memoryLimit = 4
    this.jobName =  queueLogDir + "CEU_Trio.hc"
    this.dbsnp = latestdbSNP
  }


  /****************************************************************************************
    *                3.)   VariantRecalibrator                                              *
    *****************************************************************************************/

  class VQSRBase(vcf:File) extends VariantRecalibrator with BaseCommandArguments {
    this.input :+= vcf
    this.nt = 4
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.8", "99.7", "99.5", "99.0", "98.5", "98.0", "97.0", "95.0", "90.0")
    this.memoryLimit = 8
    this.tranches_file = swapExt(outputDir, vcf, ".vcf", ".tranches")
    this.recal_file = swapExt(outputDir, vcf, ".vcf", ".recal")
    this.rscript_file = swapExt(outputDir, vcf, ".vcf", ".vqsr.R")

    // these arguments are necessary for the very high quality PCR free data
    this.percentBadVariants = 0.01
    this.minNumBad = 1000
  }

  // 3a)
  class snpRecal(snpVCF: File, useUGAnnotations: Boolean) extends VQSRBase(snpVCF) with BaseCommandArguments{
    this.resource :+= new TaggedFile( hapmapSites, "known=false,training=true,truth=true,prior=15.0" )
    this.resource :+= new TaggedFile( omni_b37_sites, "known=false,training=true,truth=true,prior=12.0" ) // truth=false on the bast practices v4
    this.resource :+= new TaggedFile( training_1000G, "known=false,training=true,prior=10.0" )	// not part of the bast practices v4
    this.resource :+= new TaggedFile( dbSNP_129, "known=true,training=false,truth=false,prior=2.0" )    // prior=6.0 on the bast practices v4
    this.use_annotation ++= List("QD", "FS", "DP", "ReadPosRankSum", "MQRankSum") // , "BaseQRankSum", "LikelihoodRankSum")
    if ( useUGAnnotations )
      this.use_annotation ++= List("HaplotypeScore")
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
    this.analysisName = "CEUTrio_VQSRs"
    this.jobName = queueLogDir + "CEUTrio.snprecal"
  }

  // 3b)
  class indelRecal(indelVCF: String) extends VQSRBase(indelVCF) with BaseCommandArguments {
    this.resource :+= new TaggedFile( indelGoldStandardCallset, "known=false,training=true,truth=true,prior=12.0" ) // known=true on the bast practices v4
    this.resource :+= new TaggedFile( latestdbSNP, "known=true,prior=2.0" )  						// not part of the bast practices v4
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
    this.analysisName = "CEUTrio_VQSRi"
    this.jobName = queueLogDir + "CEUTrio.indelrecal"
    this.use_annotation ++= List("FS", "DP", "ReadPosRankSum", "MQRankSum") // , "BaseQRankSum", "LikelihoodRankSum")
    this.maxGaussians = 4
  }

  /****************************************************************************************
    *               4.) Apply the recalibration table to the appropriate tranches         *
    ***************************************************************************************/

  class applyVQSRBase(vqsr: VariantRecalibrator) extends ApplyRecalibration with BaseCommandArguments  {
    val in = vqsr.input(0)
    this.input :+= in
    this.tranches_file = vqsr.tranches_file
    this.recal_file = vqsr.recal_file
    this.memoryLimit = 6
    this.out = swapExt(outputDir, in, ".unfiltered.vcf", ".recalibrated.vcf")
  }

  class applySnpVQSR(vqsr: VariantRecalibrator) extends applyVQSRBase(vqsr) with BaseCommandArguments {
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
    this.analysisName = "CEUTrio_AVQSRs"
    this.jobName = queueLogDir + "CEUTrio.snpcut"
    this.ts_filter_level = 99.9
  }

  class applyIndelVQSR(vqsr: VariantRecalibrator) extends applyVQSRBase(vqsr) with BaseCommandArguments {
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
    this.analysisName = "CEUTrio_AVQSRi"
    this.jobName = queueLogDir + "CEUTrio.indelcut"
    this.ts_filter_level = 99.9
  }

  def recalibrateSNPsAndIndels(snpsAndIndelsVCF: File, useUGAnnotations: Boolean): File = {
    val selectSNPs = new SelectVariants with BaseCommandArguments
    selectSNPs.V = snpsAndIndelsVCF
    selectSNPs.selectType = List(VariantContext.Type.SNP)
    selectSNPs.out = swapExt(outputDir, snpsAndIndelsVCF, ".vcf", ".snps.vcf")

    val selectIndels = new SelectVariants with BaseCommandArguments
    selectIndels.V = snpsAndIndelsVCF
    selectIndels.selectType = List(VariantContext.Type.INDEL, VariantContext.Type.MIXED, VariantContext.Type.MNP, VariantContext.Type.SYMBOLIC)
    selectIndels.out = swapExt(outputDir, snpsAndIndelsVCF, ".vcf", ".indels.vcf")

    val snpRecalibrator = new snpRecal(selectSNPs.out, useUGAnnotations)
    val snpApplyVQSR = new applySnpVQSR(snpRecalibrator)
    val indelRecalibrator = new indelRecal(selectIndels.out)
    val indelApplyVQSR = new applyIndelVQSR(indelRecalibrator)

    val recal = swapExt(outputDir, snpsAndIndelsVCF, ".unfiltered.vcf", ".recalibrated.vcf")
    val combineSNPsIndels = new CombineSNPsIndels(snpApplyVQSR.out, indelApplyVQSR.out)
    combineSNPsIndels.out = recal

    add(selectSNPs, selectIndels, snpRecalibrator, snpApplyVQSR, indelRecalibrator, indelApplyVQSR, combineSNPsIndels)

    return recal
  }

  /****************************************************************************************
    *               5) Combine Snps and Indels for UG if mode == BOTH                       *
    *****************************************************************************************/

  class CombineSNPsIndels(snpVCF:File, indelVCF:File) extends CombineVariants with BaseCommandArguments {
    this.variant :+= TaggedFile(new File(snpVCF), "snps")
    this.variant :+= TaggedFile(new File(indelVCF), "indels")
    this.filteredrecordsmergetype = FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    this.assumeIdenticalSamples = true
  }

  /****************************************************************************************
    *                6.) Variant Evaluation Base(OPTIONAL)                                  *
    *****************************************************************************************/

  // todo  -- should accept separate indel and snp vcf's, right now script will assume they're combined in one
  class Eval(evalVCF: File, prefix: String, extraStrats: Seq[String]) extends VariantEval with BaseCommandArguments{
    this.eval :+= evalVCF
    this.dbsnp = new File(dbSNP_129)
    this.doNotUseAllStandardModules = true
    this.gold = qscript.indelGoldStandardCallset
    this.comp :+= new TaggedFile( omni_b37_genotypes, "omni" )
    this.comp :+= new TaggedFile( hapmapGenotypes, "hapmap" )
    this.comp :+= new TaggedFile( omni_mono, "omni_mono" )
    this.comp :+= new TaggedFile( qscript.indelGoldStandardCallset, "GSindels" )
    this.comp :+= new TaggedFile( MVL_FP, "fp_MVL" )
    this.sample = List("NA12878", "NA12891", "NA12892")
    this.evalModule = List("TiTvVariantEvaluator", "CountVariants", "CompOverlap", "IndelSummary", "MultiallelicSummary","ValidationReport")
    this.doNotUseAllStandardStratifications = true
    this.stratificationModule = Seq("EvalRod", "CompRod", "Novelty", "FunctionalClass") ++ extraStrats
    this.num_threads = 4
    this.memoryLimit = 8
    this.out = swapExt(outputDir, evalVCF, ".vcf", prefix + ".eval")
    this.sample = samples
    this.analysisName = "CEUTrio_VE"
    this.jobName = queueLogDir + "CEUTrio.eval"
  }

  class QCSummaryRScript(@Input var vcf: File, @Input var bySite: File) extends CommandLineFunction {
    @Output var pdf: File = swapExt(outputDir, vcf, ".vcf", ".pdf")
    private val project = vcf.getName.stripSuffix(".vcf")
    def commandLine = "Rscript %s/variantCallQC_summaryTablesOnly.R %s %s %s".format(RPath, project, bySite, pdf)
  }

  class QCRScript(@Input var vcf: File, @Input var byAC: File, @Input var bySite: File, @Input var indelQC: File) extends CommandLineFunction {
    @Output var pdf: File = swapExt(outputDir, vcf, ".vcf", ".pdf")
    private val project = vcf.getName.stripSuffix(".vcf")
    def commandLine = "Rscript %s/variantCallQC.R %s %s %s %s %s".format(RPath, project, bySite, byAC, indelQC, pdf)
  }

  //------------------------------------------------------------------------------------ //
  //                      7) PhaseByTrasmission                                          //
  //------------------------------------------------------------------------------------ //
  class PBT(inputVCF: File) extends PhaseByTransmission with BaseCommandArguments {
    require(qscript.ped != null && qscript.ped.getName.endsWith(".ped"), "-ped/--pedigree must be specified as <pedigreeFileName>.ped")
    this.pedigree = Seq(qscript.ped)
    this.variant = inputVCF
    this.MendelianViolationsFile = swapExt(outputDir, variant,"vcf","mendelianViolationsFile.txt")
    this.out = swapExt(outputDir, variant,"vcf","phaseByTransmission.vcf")
    this.jobName = queueLogDir + "CEUTrio.phaseBT"
  }

  //------------------------------------------------------------------------------------ //
  //                      8) ReadBackPhasing                                             //
  //------------------------------------------------------------------------------------ //
  class RBP(inputVCF: File, inputBamFile: File) extends ReadBackedPhasing with BaseCommandArguments {
    this.input_file :+= inputBamFile
    this.variant = inputVCF
    this.respectPhaseInInput = true
    this.out = swapExt(outputDir, variant,"vcf","RBphased.vcf")
    this.jobName = queueLogDir + "CEUTrio.RBphasing"
  }

  //------------------------------------------------------------------------------------ //
  //                      9) Assessing against NA12878 KB                                //
  //------------------------------------------------------------------------------------ //
  class MyAssessAgainstKB(overrideIntervals: String = null) extends AssessNA12878 with BaseCommandArguments {
    if ( overrideIntervals != null )
      this.intervalsString = List(overrideIntervals)
    this.excludeCallset = List("CEUTrio_best_practices")
    def ext() = if ( overrideIntervals != null ) ".reviewed_intervals.txt" else ".txt"
  }


  /****************************************************************************************
    *                script                                                                 *
    *****************************************************************************************/

  def script() {
    val callSets: scala.collection.mutable.Map[String, File] = scala.collection.mutable.Map[String,File]()

    for ( pipelineVersion <- pipelineVersions ) {
      var inputBamFile = bamFile
      var name = "noName"
      if(inputBamFile.getName.endsWith(".bam.list")){
        name = inputBamFile.getName.stripSuffix(".bam.list")
      }
      else if (inputBamFile.getName.endsWith(".bam")) {
        name = inputBamFile.getName.stripSuffix(".bam")
      }

      if (useBQSR2){

        val bams = QScriptUtils.createSeqFromFile(inputBamFile)
        var recalBams = Seq.empty[File]
        for (bam <- bams) {
          val recalFile = swapExt(outputDir, bam, ".bam", ".recal.grp")
          add(new BQSR(bam, recalFile))
          val recalBam = new RecalBAM(bam, recalFile)
          add(recalBam)
          recalBams :+= recalBam.out
        }

        val mergeBamList = new ListWriterFunction
        mergeBamList.inputFiles = recalBams
        mergeBamList.listFile = "CEUTrio.recal.bam.list"
        add(mergeBamList)
        inputBamFile = mergeBamList.listFile
      }
      if (CallReduceReads){
        var bams: Seq[(File, File)] = Nil
        var reducedBams = Seq.empty[File]
        if (inputBamFile != null) {
          for (originalBam: File <- io.Source.fromFile(inputBamFile).getLines().toSeq.map(new File(_))) {
            val reducedBam: File = new File(new File(mergeBamDir, "external"), swapExt(originalBam, ".bam", ".reduced.bam").getName)
            bams :+= Tuple2(originalBam, reducedBam)
          }
        }

        for ((originalBam, reducedBam) <- bams) {
          val reduce = new ReduceReads with BadMate with RetryMemoryLimit
          reduce.memoryLimit = qscript.reduceReadsMemoryLimit
          reduce.reference_sequence = qscript.referenceFile
          reduce.input_file = Seq(originalBam)
          //reduce.intervals = Seq(qscript.intervals)
          // reduce.interval_padding = 50
          reduce.out = reducedBam
          reducedBams :+= reduce.out
          add(reduce)
        }

        val mergeReducedBamList = new ListWriterFunction
        mergeReducedBamList.inputFiles = reducedBams
        mergeReducedBamList.listFile = "CEUTrio.reduced.bam.list"
        add(mergeReducedBamList)
        inputBamFile = mergeReducedBamList.listFile

      }

      var finalCallset: File = null

      // Create the functions that we may run depending on options.
      if (pipelineVersion == 1){
        val UGgenotyper = new MyUnifiedGenotyper(name,inputBamFile)
        add(UGgenotyper)
        finalCallset = recalibrateSNPsAndIndels(UGgenotyper.out, true)
        callSets += "UG." + bamFile.getName -> finalCallset
      } else {
        val HCgenotyper = new MyHaplotypeCaller(name,inputBamFile)
        add(HCgenotyper)
        finalCallset = recalibrateSNPsAndIndels(HCgenotyper.out, false)
        callSets += "HC." + bamFile.getName -> finalCallset
      }

      if (createEvalSummaryReport || createFullPostQCReport){
        val bySample = new Eval( finalCallset, ".bySample", Seq("Sample"))
        add(bySample)

        if (createFullPostQCReport){
          val byAC = new Eval(finalCallset, ".byAC", Seq("AlleleCount"))
          add(byAC)
          val indelQC = new Eval(finalCallset, ".indelQC", Seq("Sample"))
          indelQC.stratificationModule = Seq("EvalRod", "CompRod", "Sample", "TandemRepeat", "OneBPIndel")
          indelQC.evalModule = List("IndelSummary", "IndelLengthHistogram")
          add(indelQC)
          val qc = new QCRScript(finalCallset, byAC.out, bySample.out, indelQC.out)
          add(qc)

        }
        if (createEvalSummaryReport){
          val qc = new QCSummaryRScript(finalCallset, bySample.out)  //val qc = new QCRScript(evalVCF, byAC.out, bySample.out, indelQC.out)
          add(qc)
        }
      }

      if (usePhaseByTransmission){
        val phaseBT = new PBT(finalCallset)
        add(phaseBT)
        if (useReadBackPhasing) {
          add(new RBP(phaseBT.out, inputBamFile))}
      }

      for ( assessKB <- List(new MyAssessAgainstKB(), new MyAssessAgainstKB("20:10,000,000-11,000,000"))) {
        assessKB.V :+= finalCallset
        assessKB.out = swapExt(outputDir, finalCallset, ".vcf", ".kb.gatkreport") + assessKB.ext
        assessKB.badSites = swapExt(outputDir, finalCallset, ".vcf", ".kb.badsites") + assessKB.ext + ".vcf"
        add(assessKB)
      }

      toTable(finalCallset)
    }

    // write out a combined KB assessment
    for ( assessKB <- List(new MyAssessAgainstKB(), new MyAssessAgainstKB("20:10,000,000-11,000,000"))) {
      for ( (name, calls) <- callSets )
        assessKB.V :+= new TaggedFile(calls, name)
      assessKB.out = qscript.outputDir + "/NA12878.kb.assessment.gatkreport" + assessKB.ext
      add(assessKB)
    }

    // write out a combined VCF for easy comparisons
    val cv = new CombineVariants with BaseCommandArguments
    for ( (name, calls) <- callSets )
      cv.V :+= new TaggedFile(calls, name)
    cv.out = qscript.outputDir + "/NA12878.ug_hc.combined.vcf"
    add(cv)
    toTable(cv.out)
  }

  def toTable(vcf: File) {
    val vt = new VariantsToTable with BaseCommandArguments
    vt.V :+= vcf
    vt.allowMissingData = true
    vt.raw = true
    vt.F = List("POS", "FILTER", "DP", "FS", "QD", "TYPE", "MQ", "MQRankSum", "ReadPosRankSum", "ClippingRankSum", "set")
    vt.out = swapExt(outputDir, vcf, ".vcf", ".table")
    add(vt)
  }
}

