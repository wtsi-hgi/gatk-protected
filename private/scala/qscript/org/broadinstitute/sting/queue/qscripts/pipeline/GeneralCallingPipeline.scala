package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport
import org.broadinstitute.sting.queue.extensions.gatk.BaseRecalibrator
import org.broadinstitute.sting.queue.extensions.gatk.PrintReads
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.function.ListWriterFunction

class GeneralCallingPipeline extends QScript {
  qscript =>

  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="Bam file to genotype.", shortName="I")
  var bamFile: File = _

  @Argument(shortName="V", doc="The pipeline version", required=true)
  var pipelineVersion: Int = 1

  @Argument(shortName="mode", doc="the class of variation to model: SNP, INDEL, BOTH", required=true)
  var mode: String = _

  // The following arguments are all optional.

  @Argument(shortName = "RPath", doc="RPath", required=false)
  var RPath: File = new File("../R")

  @Argument(shortName="scatterCount", doc="set the number of the scattered jobs", required=false)
  var scatterCount: Int = 0

  @Argument(shortName="useBQSR.2.0", doc="turn on a first step of using 2.0 BQSR on the input bam file", required=false)
  var useBQSR2: Boolean = false

  @Argument(shortName="usePhaseBT", doc="apply phaseByTrasmission on the output vcf", required=false)
  var usePhaseBT: Boolean = false

  @Argument(shortName="useEvalSummary", doc="create a pdf file with eval summery tables of the output vcf", required=false)
  var useEvalSummary: Boolean = false

  @Argument(shortName="outputDir", doc="output directory", required=false)
  var outputDir: String = "./tmp"

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = _

  @Argument(doc="A optional list of filter names.", shortName="filter", required=false)
  var filterNames: List[String] = Nil // Nil is an empty List, versus null which means a non-existent List.

  @Argument(shortName="noBAQ", doc="turns off BAQ calculation", required=false)
  var noBAQ: Boolean = false

  @Argument(shortName="mbq", doc="The minimum Phred-Scaled quality score threshold to be considered a good base.", required=false)
  var minimumBaseQuality: Int = -1

  @Argument(doc="An optional list of filter expressions.", shortName="filterExpression", required=false)
  var filterExpressions: List[String] = Nil

  @Argument(shortName="sample", doc="Samples to include in Variant Eval", required=false)
  var samples: List[String] = Nil

  @Argument(shortName="deletions", doc="Maximum deletion fraction allowed at a site to call a genotype.", required=false)
  var deletions: Double = -1
   
  @Input(doc="Exclude intervals list", shortName = "XL", required=false)
  var excludeIntervals: List[File] = Nil

  @Input(doc="PED file of the family", shortName="ped")
  var ped: File = _

  val dbSNP_135 = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_135.b37.vcf"  // Best Practices v4
  val hapmap = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/hapmap_3.3.b37.sites.vcf"                       // Best Practices v4
  val omni_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_2141_samples.b37.vcf"    // Best Practices v4
  val training_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"  // from the MethodDevelopmentCallingPipeline scala script
  val projectConsensus_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/ALL.wgs.projectConsensus_v2b.20101123.snps.sites.vcf"  // from the MethodDevelopmentCallingPipeline scala script
  val indelGoldStandardCallset  = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/Mills_and_1000G_gold_standard.indels.b37.vcf" // Best Practices v4
  val dbSNP_129 = "/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.vcf"
  val omni_mono: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_monomorphic_2141_samples.b37.vcf"
  val MVL_FP: String =  "/humgen/gsa-hpprojects/GATK/badSitesDB/Autism.HighMVLR.Exome.indels.vcf"



  val queueLogDir = ".qlog/"
  val name = "CEUTrio"
  val noET_key = "/humgen/gsa-hpprojects/GATK/data/gatk_user_keys/gsamembers_broadinstitute.org.key"                // TODO: remove before we make it public!!!!
	

trait BaseCommandArguments extends CommandLineGATK {
    this.logging_level = "INFO"
    phone_home = GATKRunReport.PhoneHomeOption.NO_ET   // TODO: remove before we make it public!!!!
    gatk_key = noET_key                                // TODO: remove before we make it public!!!!
    this.jobQueue = "gsa"

    this.reference_sequence = qscript.referenceFile
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    this.excludeIntervals = qscript.excludeIntervals
    this.memoryLimit = 2
    var inputBamFile = qscript.bamFile
}



trait HaplotypeCallerArguments extends BaseCommandArguments {
	
}

trait BaseBQSR extends BaseCommandArguments{

}

trait CombineIndelSnps extends BaseCommandArguments{
  var snpsOutput = qscript.outputDir + swapExt(inputBamFile, "bam","snp.recalibrated.filtered.vcf")
  var indelsOutput = qscript.outputDir + swapExt(inputBamFile, "bam","indel.recalibrated.filtered.vcf")
  var combineCallsOutput = qscript.outputDir + swapExt(inputBamFile, "bam", "both.recalibrated.filtered.vcf")
}
  
trait snpUnifiedGenotyperArguments extends BaseCommandArguments{
  
    var ugOutput = qscript.outputDir + swapExt(inputBamFile, "bam", "snp.unfiltered.vcf")
    var recalInput = ugOutput
    var applyRecalInput = ugOutput
    var tranches_File = qscript.outputDir + swapExt(inputBamFile, "bam","snp.tranches")
    var recal_File = qscript.outputDir + swapExt(inputBamFile, "bam","snp.recal")
    var recalVCF = qscript.outputDir + swapExt(inputBamFile, "bam","snp.recalibrated.filtered.vcf")
    var evalInput = recalVCF
    var evalOutput = qscript.outputDir + swapExt(recalVCF, "vcf", "eval")
}

trait indelUnifiedGenotyperArguments extends BaseCommandArguments{
    
    var ugOutput = qscript.outputDir + swapExt(inputBamFile, "bam", "indel.unfiltered.vcf")
    var recalInput = ugOutput
    var applyRecalInput = ugOutput
    var tranches_File = qscript.outputDir + swapExt(inputBamFile, "bam","indel.tranches")
    var recal_File = qscript.outputDir + swapExt(inputBamFile, "bam","indel.recal")
    var recalVCF = qscript.outputDir + swapExt(inputBamFile, "bam","indel.recalibrated.filtered.vcf")
    var	evalInput = recalVCF
    var evalOutput = qscript.outputDir + swapExt(recalVCF, "vcf", "eval")
}

// 0) BQSR v 2.0
class BQSR( bam: File, recal: File ) extends BaseRecalibrator with BaseBQSR{

  this.input_file :+= bam
	this.out = recal
	this.knownSites ++= List(new File(indelGoldStandardCallset))
 	this.knownSites ++= List(new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf"))
	this.knownSites ++= List(new File("/humgen/gsa-hpprojects/dev/carneiro/bqsr/data/projectConsensus.snps.vcf"))
	this.knownSites ++= List(new File(dbSNP_135))
	this.memoryLimit = 8
	this.qq = 0
	this.mcs = 2
	this.ics = 3
	//this.np = true
	this.javaGCThreads = 4
	this.scatterCount = 200 //if(qscript.scatterCount == 0) 200 else qscript.scatterCount
}

class RecalBAM( bam: File, recal: File ) extends PrintReads with BaseBQSR{

  this.input_file :+= bam
	this.out = qscript.outputDir + swapExt(bam, ".bam", ".subset.recal.bam")
	this.BQSR = recal
	this.memoryLimit = 6
	this.qq = 0
	this.javaGCThreads = 4
	this.scatterCount = 200 //if(qscript.scatterCount == 0) 200 else qscript.scatterCount
}

// 1.) Unified Genotyper Base
class UGBase extends UnifiedGenotyper with BaseCommandArguments {
    this.input_file :+= inputBamFile
    this.scatterCount = if(qscript.scatterCount == 0) 80 else qscript.scatterCount
    this.nt = 2
    this.dcov = 250 
    this.stand_call_conf = 30.0
    this.stand_emit_conf = 30.0
    this.D = new File(dbSNP_135)
}

  // 1a.) Call SNPs with UG
  class snpCaller extends UGBase with snpUnifiedGenotyperArguments {
    if (qscript.minimumBaseQuality >= 0)
      this.min_base_quality_score = qscript.minimumBaseQuality
    if (qscript.deletions >= 0)
      this.max_deletion_fraction = qscript.deletions
    this.out = ugOutput
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
    this.baq = if (noBAQ) {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF} else {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY}
    this.analysisName = name + "_UGs"
    this.jobName =  queueLogDir + name + ".snpcall"
  }

  // 1b.) Call Indels with UG
  class indelCaller extends UGBase with indelUnifiedGenotyperArguments {
    this.memoryLimit = 6
    this.out = ugOutput
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
    this.analysisName = name + "_UGi"
    this.jobName =  queueLogDir + name + ".indelcall"
  }

// 1.) HaplotypeCaller 
class HCBase(inputVCF: File) extends HaplotypeCaller with HaplotypeCallerArguments {
    this.excludeIntervals = excludeIntervals
    this.scatterCount = if(qscript.scatterCount == 0) 200 else qscript.scatterCount
    this.input_file :+= inputVCF
    this.out = qscript.outputDir + swapExt(inputBamFile, "bam", "HaplotypeCaller.vcf")
    this.analysisName = "HaplotypeCaller"
    this.stand_call_conf = 10.0
    this.stand_emit_conf = 10.0
    this.minPruning = 2
    this.javaGCThreads = 4
    this.memoryLimit = 8 //4	
}
 

// 3)
class VQSRBase extends VariantRecalibrator with BaseCommandArguments {
   
    this.nt = 2
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
}

// 3a)
class snpRecal extends VQSRBase with snpUnifiedGenotyperArguments{

	this.input :+= recalInput
	this.resource :+= new TaggedFile( hapmap, "known=false,training=true,truth=true,prior=15.0" ) 
        this.resource :+= new TaggedFile( omni_b37, "known=false,training=true,truth=true,prior=12.0" ) // truth=false on the bast practices v4
	this.resource :+= new TaggedFile( training_1000G, "known=false,training=true,prior=10.0" )	// not part of the bast practices v4
	this.resource :+= new TaggedFile( dbSNP_135, "known=true,training=false,truth=false,prior=2.0" )    // prior=6.0 on the bast practices v4
	this.resource :+= new TaggedFile( projectConsensus_1000G, "prior=8.0" )				// not part of the bast practices v4
  this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "DP")   //remove DP when generalizing for exomes
 	// this.use_annotation ++= List("InbreedingCoeff")   //need more then 10 samples
	this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
 	this.memoryLimit = 8
	this.tranches_file = tranches_File
	this.recal_file = recal_File
	this.rscript_file = qscript.outputDir + name + ".snp.vqsr.R"
	this.analysisName = name + "_VQSRs"
	this.jobName = queueLogDir + name + ".snprecal"

}

// 3b)
class indelRecal extends VQSRBase with indelUnifiedGenotyperArguments {

	this.input :+= recalInput
	this.resource :+= new TaggedFile(indelGoldStandardCallset, "known=false,training=true,truth=true,prior=12.0" ) // known=true on the bast practices v4
	this.resource :+= new TaggedFile( dbSNP_135, "known=true,prior=2.0" )  						// not part of the bast practices v4
	this.use_annotation ++= List("QD", "HaplotypeScore", "ReadPosRankSum", "FS")
	// this.use_annotation ++= List("InbreedingCoeff")   //need more then 10 samples   	

 	this.tranches_file = tranches_File
 	this.recal_file = recal_File
 	this.rscript_file = qscript.outputDir + name + ".indel.vqsr.R"
 	this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
 	this.analysisName = name + "_VQSRi"
 	this.jobName = queueLogDir + name + ".indelrecal"
}

//3c)
class RecalBoth(inputVCF: File) extends VariantRecalibrator with HaplotypeCallerArguments {

  this.input :+= inputVCF
  this.resource :+= new TaggedFile( hapmap, "known=false,training=true,truth=true,prior=15.0" )
  this.resource :+= new TaggedFile( omni_b37, "known=false,training=true,truth=false,prior=12.0" ) // used to be truth=true
  this.resource :+= new TaggedFile( dbSNP_135, "known=true,training=false,truth=false,prior=6.0" )
  this.resource :+= new TaggedFile( indelGoldStandardCallset, "known=true,training=true,truth=true,prior=12.0" )  //used to be known=false, prior=14.0
  this.use_annotation ++= List("QD", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "DP", "ClippingRankSum")  //remove DP for exomes
  // this.use_annotation ++= List("InbreedingCoeff")   //need more then 10 samples
  this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.BOTH
  this.allPoly = true
  this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.0", "97.9", "97.5", "97.0", "96.0","95.0", "94.0", "93.0", "92.0","90.0")
  this.recalFile = qscript.outputDir + swapExt(inputVCF, "vcf", "both.recalFile")
  this.tranchesFile = qscript.outputDir + swapExt(inputVCF, "vcf", "both.tranches")
  this.rscriptFile = qscript.outputDir + swapExt(inputVCF, "vcf", "both.R")
  this.percentBad = 0.012
  this.mG = 7
  this.ts_filter_level = 97.0
  this.memoryLimit = 8

}


// 4.) Apply the recalibration table to the appropriate tranches
class applyVQSRBase extends ApplyRecalibration with BaseCommandArguments  {
  this.memoryLimit = 6
}

class applySnpVQSR extends applyVQSRBase with snpUnifiedGenotyperArguments {

  this.input :+= applyRecalInput
  this.tranches_file = tranches_File
  this.recal_file = recal_File
  this.ts_filter_level = 99.0
  this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
  this.out = recalVCF
  this.analysisName = name + "_AVQSRs"
  this.jobName = queueLogDir + name + ".snpcut"
}

class applyIndelVQSR extends applyVQSRBase with indelUnifiedGenotyperArguments {

    this.input :+= applyRecalInput
    this.tranches_file = tranches_File
    this.recal_file = recal_File
    this.ts_filter_level = 95.0  							//best practices v4
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
    this.out = recalVCF
    this.analysisName = name + "_AVQSRi"
    this.jobName = queueLogDir + name + ".indelcut"
}

class CutBoth(inputVCF: File) extends applyVQSRBase with HaplotypeCallerArguments {

    this.input :+= inputVCF
    this.memoryLimit = 8
    this.num_threads = 1
    this.recalFile = qscript.outputDir + swapExt(inputVCF, "vcf", "both.recalFile")
    this.tranchesFile = qscript.outputDir + swapExt(inputVCF, "vcf", "both.tranches")
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.BOTH
    this.ts_filter_level = 97.0
    this.out = qscript.outputDir + swapExt(inputVCF, "vcf", "both.recalibrated.filtered.vcf")	
    this.analysisName = name + "_AVQSRb"
    this.jobName = queueLogDir + name + ".bothcut"
}

// 5) Combine Snps and Indels for UG if mode == BOTH
class CombineSNPsIndels extends CombineVariants with CombineIndelSnps {
    this.variant :+= TaggedFile(indelsOutput, "indels")
    this.variant :+= TaggedFile(snpsOutput, "snps")
    this.filteredrecordsmergetype = org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    this.assumeIdenticalSamples = true
    this.out = combineCallsOutput
}




// 6.) Variant Evaluation Base(OPTIONAL)
class EvalBase extends VariantEval with BaseCommandArguments {
      this.memoryLimit = 3
      this.comp :+= new TaggedFile(hapmap, "hapmap" )
      this.D = new File(dbSNP_129)
      this.sample = samples
}

// 6a.) SNP Evaluation (OPTIONAL) based on the recalbrated vcf
class snpEvaluation extends EvalBase with snpUnifiedGenotyperArguments {
     this.comp :+= new TaggedFile( omni_b37, "omni" )
     this.eval :+= evalInput
     this.out = evalOutput
     this.analysisName = name + "_VEs"
     this.jobName = queueLogDir + name + ".snpeval"
}

  // 6b.) Indel Evaluation (OPTIONAL)
class indelEvaluation extends EvalBase with indelUnifiedGenotyperArguments {
    this.eval :+= evalInput
    this.comp :+= new TaggedFile(indelGoldStandardCallset, "indelGS" )
    this.noEV = true
    this.evalModule = List("CompOverlap", "CountVariants", "TiTvVariantEvaluator", "ValidationReport", "IndelSummary")
    this.out = evalOutput
    this.analysisName = name + "_VEi"
    this.jobName = queueLogDir + name + ".indeleval"
}

// todo  -- should accept separate indel and snp vcf's, right now script will assume they're combined in one
class Eval(evalVCF: File, prefix: String, extraStrats: Seq[String]) extends VariantEval with BaseCommandArguments{

    this.eval :+= evalVCF
    this.dbsnp = new File(dbSNP_129)
    this.doNotUseAllStandardModules = true
    this.gold = qscript.indelGoldStandardCallset
    this.comp :+= new TaggedFile( omni_b37, "omni" )
    this.comp :+= new TaggedFile( omni_mono, "omni_mono" )
    this.comp :+= new TaggedFile( qscript.indelGoldStandardCallset, "GSindels" )
    this.comp :+= new TaggedFile( MVL_FP, "fp_MVL" )
    this.evalModule = List("TiTvVariantEvaluator", "CountVariants", "CompOverlap", "IndelSummary", "MultiallelicSummary","ValidationReport")
    this.doNotUseAllStandardStratifications = true
    this.stratificationModule = Seq("EvalRod", "CompRod", "Novelty", "FunctionalClass") ++ extraStrats
    this.num_threads = 1
    this.memoryLimit = 8
    this.out = outputDir + swapExt(evalVCF, ".vcf", prefix + ".eval")
}

class QCRScript(@Input var vcf: File, @Input var bySite: File) extends CommandLineFunction {
    @Output var pdf: File = outputDir + swapExt(vcf, ".vcf", ".pdf")
    private val project = vcf.getName.stripSuffix(".vcf")
    def commandLine = "Rscript %s/variantCallQC_summaryTablesOnly.R %s %s %s".format(RPath, project, bySite, pdf)
}

// 7) PhaseByTrasmission
class PBT(inputVCF: File) extends PhaseByTransmission with BaseCommandArguments {
  this.pedigree = Seq(qscript.ped)
  this.variant = inputVCF
  this.MendelianViolationsFile = qscript.outputDir + swapExt(variant,"vcf","mendelianViolationsFile.txt")
  this.out = qscript.outputDir + swapExt(variant,"vcf","phaseByTransmission.vcf")
  this.jobName = queueLogDir + name + ".phaseBT"
}

def script() {
  // Create the functions that we may run depending on options.
  if (pipelineVersion == 1){
	if (mode == "SNP"){
		val UGgenotyper = new snpCaller
		//val evalUnfiltered = new evalSnp  //VariantEval with UnifiedGenotyperArguments
		val snpRecalibrator = new snpRecal
		val snpApplyVQSR = new applySnpVQSR
		val evaluateSnp = new snpEvaluation
 
		add (UGgenotyper)
		add (snpRecalibrator)
		add (snpApplyVQSR)
		add (evaluateSnp)

    if (useEvalSummary){
      val evalVCF = snpApplyVQSR.out
      val bySample = new Eval( evalVCF, ".bySample", Seq("Sample"))
      add(bySample)
      val qc = new QCRScript(evalVCF, bySample.out)  //val qc = new QCRScript(evalVCF, byAC.out, bySample.out, indelQC.out)
      add(qc)
    }
    if(usePhaseBT){add(new PBT(snpApplyVQSR.out)) }
	}

	if (mode == "INDEL"){
		val UGgenotyper = new indelCaller
		val indelRecalibrator = new indelRecal
		val indelApplyVQSR = new applyIndelVQSR
		val evaluateIndel = new indelEvaluation

		add (UGgenotyper)
		add (indelRecalibrator)
		add (indelApplyVQSR)
		add (evaluateIndel)

    if (useEvalSummary){
      val evalVCF = indelApplyVQSR.out
      val bySample = new Eval( evalVCF, ".bySample", Seq("Sample"))
      add(bySample)
      val qc = new QCRScript(evalVCF, bySample.out)  //val qc = new QCRScript(evalVCF, byAC.out, bySample.out, indelQC.out)
      add(qc)
    }

    if(usePhaseBT){add(new PBT(indelApplyVQSR.out)) }
	}

  // currently running separat runs for SNPs and INDELS
	if (mode == "BOTH"){
		val snpUGgenotyper = new snpCaller
                val snpRecalibrator = new snpRecal
                val snpApplyVQSR = new applySnpVQSR
		val evaluateSnp = new snpEvaluation
		val indelUGgenotyper = new indelCaller
		val indelRecalibrator = new indelRecal
		val indelApplyVQSR = new applyIndelVQSR
	  val evaluateIndel = new indelEvaluation
    val combineSNPsIndels = new CombineSNPsIndels

		add (snpUGgenotyper)
		add (snpRecalibrator)
		add (snpApplyVQSR)
		add (evaluateSnp)
		add (indelUGgenotyper)
		add (indelRecalibrator)
		add (indelApplyVQSR)
		add (evaluateIndel)
    add (combineSNPsIndels)

    if (useEvalSummary){
      val evalVCF = combineSNPsIndels.out
      val bySample = new Eval( evalVCF, ".bySample", Seq("Sample"))
      add(bySample)
      val qc = new QCRScript(evalVCF, bySample.out)  //val qc = new QCRScript(evalVCF, byAC.out, bySample.out, indelQC.out)
      add(qc)
    }

    if(usePhaseBT){add(new PBT(combineSNPsIndels.out)) }
	}


		
 }

 else {
	if (useBQSR2){

    val bams = QScriptUtils.createSeqFromFile(qscript.bamFile)
    var recalBams = Seq.empty[File]
	  for (bam <- bams) {
		      val recalFile = qscript.outputDir + swapExt(bam, ".bam", ".recal.grp")
		      add(new BQSR(bam, recalFile))
          val recalBam = new RecalBAM(bam, recalFile)
		      add(recalBam)
          recalBams :+= recalBam.out
		}

    val mergeBamList = new ListWriterFunction
    mergeBamList.inputFiles = recalBams
    mergeBamList.listFile = "%s.recal.bam.list".format(name)
    add(mergeBamList)

    val HC_input = mergeBamList.listFile
    val HCgenotyper = new HCBase(HC_input)
    val HC_vcfFile = HCgenotyper.out
    val HCRecalibrator = new RecalBoth(HC_vcfFile)
    val HCApplyVQSR = new CutBoth(HC_vcfFile)

    add (HCgenotyper)
		add (HCRecalibrator)
		add (HCApplyVQSR)

    if (useEvalSummary){
      val evalVCF = HCApplyVQSR.out
      val bySample = new Eval( evalVCF, ".bySample", Seq("Sample"))
      add(bySample)
      val qc = new QCRScript(evalVCF, bySample.out)  //val qc = new QCRScript(evalVCF, byAC.out, bySample.out, indelQC.out)
      add(qc)
    }

    if(usePhaseBT){add(new PBT(HCApplyVQSR.out)) }

	}
	else{
		val HCgenotyper = new HCBase(qscript.bamFile)
		val HC_vcfFile = HCgenotyper.out
		val HCRecalibrator = new RecalBoth(HC_vcfFile)
		val HCApplyVQSR = new CutBoth(HC_vcfFile)	

		add (HCgenotyper)
		add (HCRecalibrator)
		add (HCApplyVQSR)

    if (useEvalSummary){
      val evalVCF = HCApplyVQSR.out
      val bySample = new Eval( evalVCF, ".bySample", Seq("Sample"))
      add(bySample)
      val qc = new QCRScript(evalVCF, bySample.out)  //val qc = new QCRScript(evalVCF, byAC.out, bySample.out, indelQC.out)
      add(qc)
    }
    if(usePhaseBT){add(new PBT(HCApplyVQSR.out)) }
	}
 }
}



// ingor this section for now!!
// This trait allows us set the variables below in one place,
// and then reuse this trait on each CommandLineGATK function below.

// currently it is hard coded in the script to run on the CEU trio as the target. In the next version of the script I will add the other built-in targets.
class Target(
              val baseName: String,
              val reference: File,
              val dbsnpFile: String,
              val hapmapFile: String,
              // val maskFile: String,
              val bamList: File,
              val intervals: String,
              // val indelTranchTarget: Double,
              val snpTrancheTarget: Double,
              val isLowpass: Boolean,
              val isExome: Boolean,
              val nSamples: Int) {

}


//val targetDataSets: Map[String, Target] = Map(
//"NA13878_test" -> new Target(	"NA13878.test",
//				qscript.referenceFile,
//				dbSNP,
//				hapmap,
//
//				)
//
//"CEUTrio_wgs_decoy" -> new Target("CEUTrio.HiSeq.WGS.b37_decoy", b37_decoy, dbSNP_b37, hapmap_b37, indelMask_b37,
//              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.list"),
//       //       new File("/humgen/gsa-hpprojects/dev/carneiro/trio/analysis/snps/CEUTrio.WEx.filtered.vcf"),                  // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
//              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg19.intervals", 90.0, 99.0, !lowPass, !exome, 3),
//)


}