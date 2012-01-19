package org.broadinstitute.sting.queue.qscripts.misc

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.text.XReadLines
import scala.collection.JavaConversions._
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.library.ipf.vcf.VCFExtractIntervals
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantQualityScore
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection
import org.broadinstitute.sting.gatk.walkers.genotyper.{GenotypeLikelihoodsCalculationModel, UnifiedGenotyperEngine}
import org.broadinstitute.sting.utils.baq.BAQ

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/18/12
 * Time: 9:36 PM
 * To change this template use File | Settings | File Templates.
 */

class T2D_VQSR extends QScript {

  @Input(shortName="BI",fullName="BroadCalls",required=true,doc="Broad unfiltered calls")
  var broadCalls : File = _

  @Input(shortName="U", fullName="Union",required=true,doc="U-Michigan Filtered SVM Union calls")
  var unionCalls : File = _

  @Input(shortName="I",fullName="Bams",required=true,doc="File listing all bam files for squaring off")
  var bamList : File = _

  @Input(shortName="NC",fullName="NumBamChunks",required=false,doc="Number of samples in each bam chunk")
  var nBamChunk : Int = 150

  val callDir : File = new File("/broad/shptmp/chartl/t2d/calling/")
  val bamDir : File = new File("/broad/shptmp/chartl/t2d/calling/bams/")
  val vqsrDir : File = new File("/broad/shptmp/chartl/t2d/vqsr/")
  val vqsrBase : File = new File(vqsrDir,"vqsr.base")

  val ref : File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  val hapMap : File = new TaggedFile("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf",
    "hapmap,vcf,known=false,training=true,truth=true,prior=15.0")
  val dbSNP : File = new TaggedFile("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.vcf",
    "dbsnp,vcf,known=true,training=false,truth=false,prior=8.0")
  // todo -- perhaps restrict ourselves only to EUR-polymorphic omni sites?
  val omni : File = new TaggedFile("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1856_samples.b37.vcf",
  "omni,vcf,known=false,training=true,truth=false,prior=12.0")

  val MAX_GAUSS = List(8,10,15,20,30)
  val DIRICH = List(0.0001,0.001,0.01)

  trait STAND_ARGS extends CommandLineGATK {
    this.reference_sequence = ref
  }

  def VQSR(inVCF: File, inList : File,  maxGauss: Int, dir: Double) : VariantRecalibrator = {
    val vqsr = new VariantRecalibrator with STAND_ARGS
    vqsr.an ++= List("QD","HaplotypeScore","MQRankSum","ReadPosRankSum","FS","MQ","InbreedingCoeff","DP")
    vqsr.mode = VariantRecalibratorArgumentCollection.Mode.SNP
    vqsr.input :+= inVCF
    vqsr.resource ++= List(hapMap,dbSNP,omni)
    vqsr.recalFile = swapExt(vqsrDir,vqsrBase, ".base", ".mg%d.dir%f.recal.txt".format(maxGauss,dir) )
    vqsr.tranchesFile = swapExt(vqsrDir,vqsrBase,".base",".mg%d.dir%f.tranche.txt".format(maxGauss,dir) )
    vqsr.intervals :+= inList
    vqsr.maxGaussians = Some(maxGauss)
    vqsr.dirichlet = Some(dir)
    vqsr.tranche ++= (900 to 1000).map(u => (u.toDouble/1000).toString)
    vqsr
  }

  def script = {

    // step 0: generate BI interval list
    val biCaInt = new VCFExtractIntervals(broadCalls,new File(callDir,"BI_calls.intervals.list"),true)
    add(biCaInt)

    // step 1: determine MI-unique calls by removing broad calls
    var remBroad = new SelectVariants with STAND_ARGS
    remBroad.variant = unionCalls
    remBroad.excludeIntervals :+= biCaInt.listOut
    remBroad.out = new File(callDir,"MI_Unique_calls.vcf")
    add(remBroad)

    // make sure to unfilter them
    var filt = new VariantFiltration with STAND_ARGS
    filt.variant = remBroad.out
    filt.out = new File(callDir,"MI_Unique_calls.noFilt.vcf")
    filt.invalidatePreviousFilters = true
    add(filt)

    val miUqInt = new VCFExtractIntervals(filt.out,new File(callDir,"MI_Unique_calls.intervals.list"),true)
    add(miUqInt)

    // step 2: merge the bam files into sub-chunks at MI-unique sites
    val printReads : List[PrintReads] = asScalaIterator(new XReadLines(bamList)).grouped(nBamChunk).zipWithIndex.map( (u : (Seq[String],Int)) => {
      val pr = new PrintReads with STAND_ARGS
      pr.input_file ++= u._1
      pr.out = new File(bamDir,"bam_chunk_%d.bam".format(u._2))
      pr.intervals :+= miUqInt.listOut
      pr
    }).toList
    addAll(printReads)

    // step 3: recall using alleles only at MI unique sites
    var callMI = new UnifiedGenotyper with STAND_ARGS
    callMI.alleles = filt.out
    callMI.input_file ++= printReads.map( u => u.out )
    callMI.intervals :+= miUqInt.listOut
    callMI.genotyping_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
    callMI.out = new File(callDir,"BI_Calls_MI_Unique_Sites.vcf")
    callMI.stand_call_conf = Some(4.0)
    callMI.stand_emit_conf = Some(4.0)
    callMI.baq = BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    callMI.scatterCount = 50

    add(callMI)

    // step 4: combine initial broad calls with the unique ones
    var combineBroad = new CombineVariants with STAND_ARGS
    combineBroad.variant :+= broadCalls
    combineBroad.variant :+= callMI.out
    combineBroad.intervals :+= broadCalls
    combineBroad.intervals :+= miUqInt.listOut
    combineBroad.out = new File(callDir,"BI_Union_calls.vcf")

    add(combineBroad)

    // step 5: make an interval list
    var combInt = new VCFExtractIntervals(combineBroad.out,new File(vqsrDir,"Combined_Calls.intervals.list"),true)
    add(combInt)

    // step 6: run the VQSR for various selections of parameters
    for ( mg <- MAX_GAUSS ) {
      for ( dir <- DIRICH ) {
        add(VQSR(combineBroad.out,combInt.listOut,mg,dir))
      }
    }
  }
}