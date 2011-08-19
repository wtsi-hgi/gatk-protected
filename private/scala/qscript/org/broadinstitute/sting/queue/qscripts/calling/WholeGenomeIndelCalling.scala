package org.broadinstitute.sting.queue.qscripts.calling

import org.broadinstitute.sting.commandline.Argument._
import org.broadinstitute.sting.commandline.Input._
import org.broadinstitute.sting.pipeline.Pipeline
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.commandline.Input._
import org.broadinstitute.sting.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk._
import java.io.File
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk.RodBind._
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection
;
/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 8/4/11
 * Time: 11:04 AM
 * To change this template use File | Settings | File Templates.
 */

class WholeGenomeIndelCalling extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=false)
  var gatkJar: File = new File("/humgen/gsa-scr1/delangel/GATK/Sting_unstable/dist/GenomeAnalysisTK.jar")

  @Input(doc="output path", shortName="outputDir", required=true)
  var outputDir: String = _

  @Input(doc="run name", shortName="runName", required=true)
  var runName: String = _

  @Input(doc="queue", shortName="queue", required=false)
  var jobQueue: String = "hour"

  @Input(doc="Do only one chromosome", shortName="onlyOneChr", required=false)
  var onlyOneChr: Boolean = false

  @Input(doc="the chromosome to process", shortName="chrToProcess", required=false)
  var chrToProcess: Int = 20

  @Input(doc="scatter count", shortName="scatterCount", required=false)
  var scatterCount: Int = 50


  @Argument(shortName = "R", doc="B37 reference sequence: defaults to broad standard location", required=false)
  var reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Input(doc="path to tmp space for storing intermediate bam files", shortName="outputTmpDir", required=false)
  var outputTmpDir: String = "/broad/shptmp/delangel"

  @Input(doc="BAM file or list of bam files to call", shortName="bam", required=true)
  var bamList: String = _

  @Argument(shortName = "truth", doc="VQSR truth file", required=false)
  var truthFile: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Mills_Devine_Indels_2011/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.vcf"  )

  @Argument(shortName = "training", doc="VQSR training file", required=false)
  var trainingFile: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Mills_Devine_Indels_2011/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.vcf"  )

  @Argument(shortName = "intervals", doc="intervals", required=false)
  val myIntervals: String = null;

  @Argument(shortName = "doOnlyVQSR", doc="doOnlyVQSR", required=false)
  val doOnlyVQSR: Boolean = false;

  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560)
  //  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,3000000,48129895,51304566,155270560)
  val COMP_MULLIKIN =  "/humgen/gsa-hpprojects/GATK/data/Comparisons/StandardForEvaluation/b37/comps/NA12878.DIPline.NQScm.expanded.chr20.b37.minReads_2_or_gt2bp.vcf"
  val COMP_MILLS =  "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Mills_Devine_Indels_2011/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.vcf"
  val COMP_MILLS_CHIP =  "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Mills_Devine_Indels_2011/indel_genotype_data_hg19_annotated_polymorphic.vcf"

  private var pipeline: Pipeline = _
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = Some(2)
    this.jobQueue = qscript.jobQueue
    if (qscript.myIntervals != null) {
      this.intervalsString = List(myIntervals);
    }
  }

  def script = {
    var projectBase:String = qscript.outputDir + qscript.runName

    var rawVCFIndels = new File(projectBase + ".raw.vcf")
    val callIndels = new UnifiedGenotyper with CommandLineGATKArgs
     callIndels.out = rawVCFIndels
    // callIndels.dcov = 50
     callIndels.stand_call_conf = 4.0
     callIndels.stand_emit_conf = 4.0
     callIndels.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
  //   callIndels.jobName = qscript.outputTmpDir + "/calls/" + qscript.runName
     callIndels.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
     callIndels.dbsnp =  qscript.dbSNP
     callIndels.sites_only = false
     callIndels.scatterCount = qscript.scatterCount
     callIndels.input_file :+= qscript.bamList
//     callIndels.nt=Some(qscript.nt)

    if (!qscript.doOnlyVQSR) {
      add(callIndels)
    }

    val vr = new VariantRecalibrator with CommandLineGATKArgs
    vr.input :+= callIndels.out
    vr.truth :+= new TaggedFile( qscript.truthFile,"prior=15.0")
    vr.training :+=  new TaggedFile( qscript.trainingFile,"prior=12.0")
    vr.known :+= new TaggedFile(qscript.dbSNP,"prior=3.0")
    //vr.trustAllPolymorphic = true
    vr.mode = VariantRecalibratorArgumentCollection.Mode.INDEL

//    vr.use_annotation = List("QD", "HaplotypeScore",  "ReadPosRankSum","FS","InbreedingCoeff")
   // todo - InbreedingCoeff not appropriate for single sample calling, should be extended to command line argument
    vr.use_annotation = List("QD", "HaplotypeScore",  "ReadPosRankSum","FS")
    vr.TStranche = List("99.0", "97.0","95.0")
    vr.tranches_file = projectBase + ".tranches"
    vr.recal_file = projectBase + ".recal"
    vr.rscriptFile = projectBase + ".plots.R"
    vr.jobOutputFile = vr.recal_file + ".out"
    vr.memoryLimit = 32
    vr.nt = 16
    vr.mG = 8
    vr.percentBad = 0.03
    //vr.std = 14
    add(vr)

    for (tranche <- vr.TStranche) {
      val ar = new ApplyRecalibration with CommandLineGATKArgs
      ar.input :+= (callIndels.out)
      ar.tranches_file = vr.tranches_file
      ar.recal_file = vr.recal_file
      ar.ts_filter_level = tranche.toDouble
      ar.out = projectBase + ".recalibrated." + tranche + ".vcf"
      ar.jobOutputFile = ar.out + ".out"
      ar.memoryLimit = 32
      ar.mode = VariantRecalibratorArgumentCollection.Mode.INDEL
      add(ar)

      val eval = new VariantEval with CommandLineGATKArgs
      eval.tranchesFile = vr.tranches_file
      eval.eval :+= ( ar.out)
      eval.dbsnp = qscript.dbSNP
      eval.doNotUseAllStandardStratifications = true
      eval.doNotUseAllStandardModules = true
      eval.evalModule = List("SimpleMetricsByAC", "CountVariants","CompOverlap")
      eval.stratificationModule = List("EvalRod", "CompRod", "Novelty","Sample")
      eval.out = swapExt(ar.out, ".vcf", ".eval")
      eval.jobOutputFile = eval.out + ".out"
      eval.comp :+= new TaggedFile(COMP_MILLS,"Mills")
      eval.comp :+= new TaggedFile(COMP_MILLS_CHIP,"MillsChip")
      eval.comp :+= new TaggedFile(COMP_MULLIKIN,"Mullikin")

      eval.memoryLimit = 32
      add(eval)
    }

  }

}