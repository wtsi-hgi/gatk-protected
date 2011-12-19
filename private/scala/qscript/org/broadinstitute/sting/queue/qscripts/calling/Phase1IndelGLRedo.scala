package org.broadinstitute.sting.queue.qscripts.calling

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.gatk.TaggedFile._
import org.broadinstitute.sting.queue.{QException, QScript}
import scala.io.Source

class Phase1IndelGLRedo extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="output path", shortName="outputDir", required=true)
  var outputDir: String = _

  @Input(doc="queue", shortName="queue", required=false)
  var jobQueue: String = "hour"

  @Input(doc="the chromosome to process", shortName="onlyOneChr", required=false)
  var onlyOneChr: Boolean = false

  @Input(doc="the chromosome to process", shortName="chrToProcess", required=false)
  var chrToProcess: Int = 20

  @Input(doc="sample bam list", shortName="sampleBamList", required=true)
  var sampleBamList: String = "/humgen/gsa-scr1/delangel/Phase1GLFix/allBamsPerSample.bam_list"

  @Input(doc="indel alleles", shortName="indelAlleles", required=false)
  var indelAlleles: String = "/humgen/1kg/processing/production_wgs_phase1/consensus/ALL.indels.combined.chr20.vcf"

  @Input(doc="scatterCount", shortName="scatterCount", required=false)
  var variantCallerScatterCount: Int = 1

  private val reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")

  val populations = List("ASW","CEU","CHB","CHS","CLM","FIN","GBR","IBS","JPT","LWK","MXL","PUR","TSI","YRI")



  val tabixCMD = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/tabix/tabix_0.2.2/bin/tabix"
  val bgzipCMD = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/tabix/tabix_0.2.2/bin/bgzip"

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = Some(2)
    this.jobQueue = qscript.jobQueue

  }
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 2
//    this.isIntermediate = true
  }

  def script = {


    var chrList:List[Int] = List(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
    if (qscript.onlyOneChr) {
      chrList = List(qscript.chrToProcess)
    }

    for(chr <- chrList) {
      var chrStr:String = "%d".format(chr)
      if (chr== 23)
        chrStr = "X"

      // parse input bam file, each line will be one BAM
      val indelCombine = new CombineVariants with CommandLineGATKArgs
      //indelCombine.intervalsString  :+=chrStr
      indelCombine.out = new File(qscript.outputDir + "/calls/ALL.chr"+chr+".VQSR_V3_biallelic_GLs_corrected.phase1.20101123.indels.genotypes.vcf")
     // indelCombine.jobName = qscript.outputDir + "/tmp/ALL.chr"+chr + ".VQSR_V3_biallelic_GLs_corrected.phase1.20101123.indels.genotypes.vcf"
      indelCombine.jobQueue = "gsa"
      indelCombine.memoryLimit = Some(6)
      indelCombine.genotypeMergeOptions = org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.GenotypeMergeType.UNSORTED

      for(line <- Source.fromFile(qscript.sampleBamList).getLines()) {
        val lArray =  line.split("\t")
        val sample = lArray(0)
        val bamName = lArray(1)
        val rawVCFindels = new File(qscript.outputDir + "/calls/chr" + chrStr + "/"+ sample +  ".VQSR_V3_biallelic_GLs_corrected.phase1.20101123.indels.genotypes.vcf")

        val callIndels = new UnifiedGenotyper with CommandLineGATKArgs
        callIndels.out = rawVCFindels
        callIndels.input_file :+= bamName
        callIndels.dcov = 50
        callIndels.stand_call_conf = 0.0
        callIndels.stand_emit_conf = 0.0
        callIndels.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
        callIndels.jobName = qscript.outputDir + "/tmp/chr"+chrStr + "/" +sample  + ".VQSR_V3_biallelic_GLs_corrected.phase1.20101123.indels.genotypes"
        callIndels.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
        callIndels.genotyping_mode = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
        callIndels.out_mode = org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
        callIndels.alleles = qscript.indelAlleles
        callIndels.dbsnp = qscript.dbSNP
        callIndels.ignoreSNPAlleles = true
        callIndels.intervalsString  :+= chrStr
        //        callIndels.scatterCount = qscript.variantCallerScatterCount

        add(callIndels)

        indelCombine.variant :+= TaggedFile(callIndels.out, sample)
      }
      //indelCombine.scatterCount = 50
      // bgzip and tabix index chromosome file
      //val zip = new bgZip(indelCombine.out)
      //val idx = new tabixIndex(indelCombine.out)
      //add(indelCombine)
      //add(zip)
      //add(idx)
    }

  }
  case class bgZip (inVCF: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="vcf file to compress") var inv = inVCF
    def commandLine = qscript.bgzipCMD +  " " + inVCF
    this.jobName =  qscript.outputDir + "/tmp/"+swapExt(inVCF, ".vcf", ".bgzip")
  }

  case class tabixIndex (inVCF: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="vcf file to index") var inv = inVCF
    def commandLine = qscript.tabixCMD + " -p vcf -s 1 -b 2 -e 2 -c \\#  " + inVCF + ".gz"
    this.jobName =  qscript.outputDir + "/tmp/"+swapExt(inVCF, ".vcf", ".tabix")
  }

}
