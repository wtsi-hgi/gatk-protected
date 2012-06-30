package org.broadinstitute.sting.queue.qscripts.calling

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.gatk.CombineVariants
import org.broadinstitute.sting.queue.extensions.gatk.CommandLineGATK

import org.broadinstitute.sting.queue.{QException, QScript}
import scala.io.Source
import org.broadinstitute.sting.queue.function.scattergather.{GatherFunction, CloneFunction, ScatterFunction}
import org.broadinstitute.sting.commandline.ArgumentSource

class Phase1IndelRedoCombining extends QScript {
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
      indelCombine.intervalsString  :+=chrStr
      indelCombine.out = new File(qscript.outputDir + "/calls/ALL.chr"+chrStr+".VQSR_V3_biallelic_GLs_corrected.phase1.20101123.indels.genotypes.vcf")
      //indelCombine.jobOutputFile =  indelCombine.out + ".out"
      indelCombine.jobQueue = "gsa"
      indelCombine.memoryLimit = Some(6)
      indelCombine.genotypeMergeOptions = org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.GenotypeMergeType.UNSORTED

      for(line <- Source.fromFile(qscript.sampleBamList).getLines()) {
        val lArray =  line.split("\t")
        val sample = lArray(0)
        val rawVCFindels = new File(qscript.outputDir + "/calls/chr" + chrStr + "/"+ sample +  ".VQSR_V3_biallelic_GLs_corrected.phase1.20101123.indels.genotypes.vcf")

        indelCombine.variant :+= TaggedFile(rawVCFindels, sample)
      }
      indelCombine.scatterCount = 50

      indelCombine.setupScatterFunction = {
        case scatter: ScatterFunction =>
          scatter.commandDirectory = new File(qscript.outputDir+"/calls/ScatterGather/chr"+chrStr)
          scatter.jobOutputFile = new File(qscript.outputDir+"/calls/ScatterGather/chr"+chrStr +"/Scatter.out")
      }
      indelCombine.setupCloneFunction = {
        case (clone: CloneFunction, index: Int) =>
          clone.commandDirectory = new File(qscript.outputDir + "/calls/ScatterGather/chr"+chrStr  + "/Scatter_%s".format(index))
          clone.jobOutputFile =new File(qscript.outputDir + "/calls/ScatterGather/chr"+chrStr  + "/Scatter_%s.out".format(index))
      }
      indelCombine.setupGatherFunction = {
        case (gather: GatherFunction, source: ArgumentSource) =>
          gather.commandDirectory = new File(qscript.outputDir + "/calls/ScatterGather/chr"+chrStr  +"/Gather_%s".format(source.field.getName))
          gather.jobOutputFile = new File(qscript.outputDir +  "/calls/ScatterGather/chr"+chrStr  +"/Gather_%s.out".format(source.field.getName))

      }


      indelCombine.scatterGatherDirectory = new File(qscript.outputDir + "/tmp/ScatterGather/chr"+chrStr)
      add(indelCombine)
    }

  }
}
