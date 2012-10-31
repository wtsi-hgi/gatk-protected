package org.broadinstitute.sting.queue.qscripts.calling

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.walkers.genotyper.{UnifiedGenotyperEngine, GenotypeLikelihoodsCalculationModel}
import org.broadinstitute.sting.gatk.downsampling.DownsampleType

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 6/27/12
 * Time: 3:46 PM
 * To change this template use File | Settings | File Templates.
 */

class LargeScaleValidationCallingSingle extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="intervals to process", shortName="intervals", required=false)
  var intervals: File = _

  @Argument(doc="output path", shortName="outputDir", required=true)
  var outputDir: String = _

  @Argument(doc="base output filename", shortName="baseName", required=false)
  var baseName: String = ""

  @Argument(doc="scatterCount", shortName="sc", required=false)
  var variantCallerScatterCount: Int = 1

  @Argument(doc="chromosomes in pool", shortName="ploidy", required=false)
  var ploidy: Int = 24

  @Argument(doc="doRef", shortName="doRef", required=false)
  var doRefSample: Boolean = false

  //  @Argument(fullName="standard_min_confidence_threshold_for_emitting_and_calling", shortName="stand_conf", doc="The minimum phred-scaled confidence threshold at which variants should be emitted and called", required=false)
  //  var stand_conf: Option[Double] = 10.0

  @Argument(doc="doSingle", shortName="doSingle", required=false)
  var doSingle: Boolean = false

  @Argument(doc="validation set", shortName="vs", required=false)
  var vs: String = _


  val bamListRef: File = new File("/humgen/gsa-hpprojects/dev/largeScaleValidation/dataAnalysis/Pools.bam.list")
  val bamListNoRef: File = new File("/humgen/gsa-hpprojects/dev/largeScaleValidation/dataAnalysis/Pools.noRef.bam.list")

  val originalSites: File = new File("/humgen/gsa-hpprojects/dev/largeScaleValidation/ALL.wgs.allCombinedValidationSites.ACannotated.corr.sites.vcf")

  private val tmpDir: File = new File("/broad/hptmp/delangel/tmp/")
  private val reference: File = new File("/humgen/gsa-hpprojects/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.vcf")
  private val oneKGRelease: File = new File("/humgen/1kg/DCC/ftp/release/20110521/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz")
  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = 2
    this.jobTempDir = qscript.tmpDir
    this.jobQueue = "gsa"
    //    if (!qscript.intervals.isEmpty) this.intervals :+= qscript.intervals


  }

  class PPC(val callName: String, val intervalFile: File) extends UnifiedGenotyper with CommandLineGATKArgs {
    this.scatterCount = qscript.variantCallerScatterCount
    this.sample_ploidy = Some(qscript.ploidy)
    this.dbsnp = qscript.dbSNP
    var refStr:String = ""

    if (doRefSample) {
      this.reference_sample_name = "NA12878"
      this.input_file :+= qscript.bamListRef
      refStr = "withRef"
    }
    else {
      this.input_file :+= qscript.bamListNoRef
      refStr = "noRef"

    }
    this.out = qscript.outputDir + "/"+qscript.baseName + "."+callName+"."+refStr+".vcf"

    this.gt_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY
    this.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_CONFIDENT_SITES
//    this.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES

    this.intervals = Seq(intervalFile)

    this.ignoreLane = true
    this.maxAltAlleles = Some(3)  // memory usage will overflow without this
    this.dt = DownsampleType.NONE

    this.standard_min_confidence_threshold_for_emitting = Some(5.0)
    this.standard_min_confidence_threshold_for_calling= Some(30.0)


  }
  class SNPPC(callName: String, intervalFile: File) extends PPC(callName, intervalFile) {
    this.glm = GenotypeLikelihoodsCalculationModel.Model.SNP
    this.referenceCalls = new File("/humgen/gsa-hpprojects/NA12878Collection/callsets/snps/NA12878.HiSeq.WGS.b37.recalibrated.99_5_cut_for_heng.vcf")
    this.max_deletion_fraction=.1
  }

  class IndelPC(callName: String, intervalFile: File) extends PPC(callName, intervalFile) {
    this.glm = GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.minIndelFrac = Some(0.01)
    this.referenceCalls = new File("/humgen/gsa-scr1/delangel/IndelGoldSet/CEUTrio.HiSeq.WGS.b37_decoy.recal.ts_95.vcf")
    this.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES

  }


  class Eval(evalVCF: File, compVCF: File) extends VariantEval with CommandLineGATKArgs {
    this.eval :+= evalVCF
    this.dbsnp = qscript.dbSNP
    this.doNotUseAllStandardModules = true
    this.evalModule = List("CountVariants", "ValidationReport")
    this.doNotUseAllStandardStratifications = true
    this.stratificationModule = Seq("EvalRod", "Filter") //++ extraStrats
    this.out = swapExt(evalVCF, ".vcf", ".eval")
    this.ploidy = Some(qscript.ploidy)
    this.numSamples = Some(92)

    if (compVCF != null) {
      this.comp :+= compVCF
      this.intervals = Seq(compVCF)
    }
  }

  class Filt(inputVCF: File) extends VariantFiltration with CommandLineGATKArgs {
    this.V = inputVCF
    this.filterExpression = Seq("DP<5000")
    this.filterName = Seq("LowDepth")
    if (qscript.doRefSample) {
      this.filterExpression :+= "REFDEPTH<500"
      this.filterName :+= "LowReferenceSampleDepth"
    }
    this.out = swapExt(inputVCF, ".vcf",".filtered.vcf")
  }
  class SampleEval(evalVCF: File, compVCF: File) extends Eval(evalVCF, compVCF) {
    this.stratificationModule :+= "Sample"
    // this.num_threads = qscript.num_threads
    // this.memoryLimit = 8
    this.out = swapExt(evalVCF, ".vcf", ".bySample.eval")
  }

  class ACEval(evalVCF: File, compVCF: File) extends Eval(evalVCF, compVCF) {
    this.stratificationModule :+= "AlleleCount"
    this.out = swapExt(compVCF, ".vcf", ".byAC.eval")
  }

  class Annot(inputVCF: File) extends VariantAnnotator with CommandLineGATKArgs {
    this.V = inputVCF
    this.out = swapExt(inputVCF, ".vcf",".annotated.vcf")
    this.resource :+= qscript.originalSites
    this.E = Seq("resource.oneKGAC","resource.set" )

  }

  class VToT(inputVCF: File) extends VariantsToTable with CommandLineGATKArgs {
    this.V :+= inputVCF
    this.out = swapExt(inputVCF, ".vcf",".table")
    this.F = Seq("CHROM","POS","REF","ALT","FILTER","AC","QUAL","resource.oneKGAC","resource.set","DP","TYPE","NCALLED")
    this.allowMissingData = true
    this.showFiltered = true

  }

  def script = {

    var jobList:Seq[PPC] = Seq()


    //    if ()
    jobList:+= (new SNPPC("omniMono","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_OmniMono.SNP.sites.interval_list"))
    jobList:+= (new SNPPC("omniPoly","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_OmniPoly.SNP.sites.interval_list"))
    jobList:+= (new SNPPC("exomeChip","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_ExomeChip.SNP.sites.interval_list"))
    jobList:+= (new IndelPC("millsPoly","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_MillsGenotypeInPhase1.INDEL.sites.interval_list"))
    jobList:+= (new SNPPC("lostToImputation","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.2000_lost_to_Imputation.sites.interval_list"))
    //add(new SNPPC(("multiAllelicSNPs","/humgen/gsa-hpprojects/dev/largeScaleValidation/inputSets/triallelics.EricBanks.vcf")
    jobList:+= (new SNPPC("LOFSNP","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/LOF.DanielMacArthur_20120910.sorted.fixed.snps.interval_list"))
    jobList:+= (new IndelPC("LOFINDEL","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/LOF.DanielMacArthur_20120910.sorted.fixed.indels.interval_list"))
    jobList:+= (new IndelPC("afIndels","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.5000_validation_sites_AF_distributed.indels.sites.interval_list"))
    jobList:+= (new IndelPC("unifIndels","/humgen/gsa-hpprojects/dev/largeScaleValidation//outputVCFs/ALL.wgs.5000_validation_sites_Uniformly_distributed.indels.sites.interval_list"))
    jobList:+= (new SNPPC("afSNPs","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.8000_validation_sites_AF_distributed.snp.sites.interval_list"))
    jobList:+= (new SNPPC("unifSNPs","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.8000_validation_sites_Uniformly_distributed.snp.sites.interval_list"))

    for (job <- jobList) {
      add(job)
      val filt = new Filt(job.out)
      add(filt)
      val annot = new Annot(filt.out)
      annot.intervals :+= filt.out
      add(annot)
      val vtot = new VToT(annot.out)
      add(vtot)
    }

    /*
// just to get raw counts of ref sites/variant sites on the whole baits
val ve = new Eval(filter.out, null)
add(ve)

// per-AC validation just in called sites
val veAC = new ACEval(qscript.oneKGRelease, filter.out)
add(veAC)    */
  }

}

