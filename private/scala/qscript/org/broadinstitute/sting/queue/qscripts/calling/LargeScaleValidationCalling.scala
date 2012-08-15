package queue.qscripts.calling

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.walkers.genotyper.{UnifiedGenotyperEngine, AlleleFrequencyCalculationModel, GenotypeLikelihoodsCalculationModel}

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 6/27/12
 * Time: 3:46 PM
 * To change this template use File | Settings | File Templates.
 */

class LargeScaleValidationCalling extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="intervals to process", shortName="intervals", required=false)
  var intervals: String = ""

  @Input(doc="output path", shortName="outputDir", required=true)
  var outputDir: String =  _

  @Input(doc="input bAM list", shortName="bamList", required=true)
  var bamList: File = _

  @Input(doc="base output filename", shortName="baseName", required=false)
  var baseName: String = ""

  @Input(doc="scatterCount", shortName="scatterCount", required=false)
  var variantCallerScatterCount: Int = 1

  @Input(doc="chromosomes in pool", shortName="ploidy", required=false)
  var ploidy: Int = 24

  private val tmpDir: File = new File("/broad/hptmp/delangel/tmp/")
  private val reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.vcf")

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = 2
    this.jobTempDir = qscript.tmpDir
    this.jobQueue = "gsa"
    this.intervalsString :+= qscript.intervals
  }
  class PPC(val callName: String, val allelesFile: String) extends UnifiedGenotyper with CommandLineGATKArgs {
    this.scatterCount = qscript.variantCallerScatterCount
    this.input_file :+= qscript.bamList
    this.sample_ploidy = Some(qscript.ploidy)
    this.dbsnp = qscript.dbSNP
    this.out = qscript.outputDir + "/"+qscript.baseName + "."+callName+".vcf"
    this.reference_sample_name = "NA12878"
    //      this.referenceCalls = new File("/humgen/gsa-hpprojects/NA12878Collection/callsets/snps/NA12878.HiSeq.WGS.b37.recalibrated.99_5_cut_for_heng.vcf")
    this.gt_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
    this.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
    this.alleles = new File(allelesFile)
    this.ignoreLane = true
    this.maxAltAlleles = Some(1)  // memory usage will overflow without this

  }
  class SNPPC(callName: String, allelesFile: String) extends PPC(callName, allelesFile) {
    this.glm = GenotypeLikelihoodsCalculationModel.Model.SNP
    this.referenceCalls = new File("/humgen/gsa-hpprojects/NA12878Collection/callsets/snps/NA12878.HiSeq.WGS.b37.recalibrated.99_5_cut_for_heng.vcf")
  }

  class IndelPC( callName: String,  allelesFile: String) extends PPC(callName, allelesFile) {
    this.glm = GenotypeLikelihoodsCalculationModel.Model.INDEL
    this.referenceCalls = new File("/humgen/gsa-scr1/delangel/IndelGoldSet/CEUTrio.HiSeq.WGS.b37_decoy.recal.ts_95.vcf")
  }

  class BothPC( callName: String,  allelesFile: String) extends PPC(callName, allelesFile) {
    this.glm = GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.referenceCalls = new File("/humgen/gsa-scr1/delangel/IndelGoldSet/CEUTrio.HiSeq.WGS.b37_decoy.recal.ts_95.snp_indel_combined.vcf")
  }

  class Eval(evalVCF: File) extends VariantEval with CommandLineGATKArgs {
    this.eval :+= evalVCF
    this.dbsnp = qscript.dbSNP
    this.doNotUseAllStandardModules = true
    this.evalModule = List("TiTvVariantEvaluator", "CountVariants", "CompOverlap", "IndelSummary")
    this.doNotUseAllStandardStratifications = true
    this.stratificationModule = Seq("EvalRod", "CompRod", "Novelty", "TandemRepeat") //++ extraStrats
    // this.num_threads = qscript.num_threads
    // this.memoryLimit = 8
    this.out = swapExt(evalVCF, ".vcf", ".eval")
  }

  def script = {



    val omnimono = new SNPPC("omniMono","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_OmniMono.SNP.sites.vcf")
    val omnipoly = new SNPPC("omniPoly","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_OmniPoly.SNP.sites.vcf")
    val exomechippoly = new SNPPC("exomeChip","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_ExomeChip.SNP.sites.vcf")
    val millspoly = new IndelPC("millsPoly","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_MillsGenotypeInPhase1.INDEL.sites.vcf")
    val multiAllelicIndels = new IndelPC("multiallelicIndels","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.2000_validation_sites_multiAllelicIndels.sites.vcf")
    multiAllelicIndels.maxAltAlleles = Some(3)
    val lostToImputation = new SNPPC("lostToImputation","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.2000_lost_to_Imputation.sites.vcf")
    val multiAllelicSNPs = new SNPPC("multiAllelicSNPs","/humgen/gsa-hpprojects/dev/largeScaleValidation/inputSets/triallelics.EricBanks.vcf")
    multiAllelicSNPs.maxAltAlleles = Some(3)
    val lof = new BothPC("LOF","/humgen/gsa-hpprojects/dev/largeScaleValidation/inputSets/LOF.DanielMacArthur.vcf")
    val afIndels = new IndelPC("afIndels","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.5000_validation_sites_AF_distributed.indels.sites.vcf")
    val unifIndels = new IndelPC("unifIndels","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.5000_validation_sites_Uniformly_distributed.indels.sites.vcf")
    val afSNPs = new SNPPC("afSNPs","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.8000_validation_sites_AF_distributed.snp.sites.vcf")
    val unifSNPs = new SNPPC("unifSNPs","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.8000_validation_sites_Uniformly_distributed.snp.sites.vcf")
    
    add(omnimono)
    add(omnipoly)
    add(exomechippoly)
    add(millspoly)
    add(lof)
    add(multiAllelicIndels)
    add(lostToImputation)
    add(multiAllelicSNPs)
    add(afIndels)
    add(unifIndels)
    add(afSNPs)
    add(unifSNPs)
    
    // run VE to assess results
    val veOmniMono = new Eval(omnimono.out)
    add(veOmniMono)

    // run VE to assess results
    val veOmniPoly = new Eval(omnipoly.out)
    add(veOmniPoly)

    val veMillsPoly = new Eval(millspoly.out)
    add(veMillsPoly)

    // run VE to assess results
    val veExomeChip = new Eval(exomechippoly.out)
    add(veExomeChip)
    // run VE to assess results
    val veLOF = new Eval(lof.out)
    add(veLOF)
    // run VE to assess results
    val veAFIndels = new Eval(afIndels.out)
    add(veAFIndels)
    // run VE to assess results
    val veUnifIndels = new Eval(unifIndels.out)
    add(veUnifIndels)

    val veAFSNPs = new Eval(afSNPs.out)
    add(veAFSNPs)
    // run VE to assess results
    val veUnifSNPs = new Eval(unifSNPs.out)
    add(veUnifSNPs)
    val veLostToImputation = new Eval(lostToImputation.out)
    add(veLostToImputation)

  }

}

