package org.broadinstitute.sting.queue.qscripts.analysis

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction
import org.broadinstitute.sting.utils.baq.BAQ

class ReducedBAMEvaluation extends QScript {



  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "bam", doc="BAM", required=true)
  val bam: File = null;

  @Argument(shortName = "reduceIntervals", doc="Interval to reduce at", required=false)
  val REDUCE_INTERVAL: String = null;

  @Argument(shortName = "dcov", doc="dcov", required=false)
  val DCOV: Int = 250;

  @Argument(shortName = "minimalVCF", doc="", required=false)
  val minimalVCF: Boolean = false;

  @Argument(shortName = "CS", doc = "ContextSize", required = false)
  val CS: String = "7-7";

  @Argument(shortName = "ADAV", doc = "AverageDepthAtVariableSites", required = false)
  val ADAV: String = "30-30";

  @Argument(shortName = "MBRC", doc = "MinimumBasepairsForRunningConsensus", required = false)
  val MBRC: String = "50-50";

  val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_132.b37.vcf")

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = referenceFile;
    this.memoryLimit = 4
  }

  //Makes the file name sorting less confusing "CS-9" -> "CS-09"
  def toStringLength(list: List[Int], num: Int): String = {
    val length = list.max.toString().length()
    var Num = num.toString
    if (Num.length() >= length)
      Num
    else {
      while(Num.length() < length) { Num = "0" + Num }
      Num
    }
  }

  trait CoFoJa extends JavaCommandLineFunction {
    override def javaOpts = super.javaOpts //  + " -javaagent:lib/cofoja.jar"
  }

  //Convert CS, ADAV, MBRC to lists
  object Int {
  def unapply(s : String) : Option[Int] = try {
    Some(s.toInt)
  } catch {
    case _ : java.lang.NumberFormatException => None
  }
}

  def getValue(s: String): Int = s match {
    case Int(x) => x
    case _ => error("String conversion FAILED! : Not a number!")
  }

  /* todo: allow user to type 7-8,10,11-12. SiIngle value w/o range is not working */
  // Take out commas to make individual elements and convert dashes to there respective ranges
  def breakUpString(string: String): List[Int] = {
    var result: List[Int] = List();
    //if (string.matches("[0-9,-]"))
    val splitString = string.split(',')
    println(splitString)
    for{
      num <- splitString
    } yield {
      println(num)
      if ( ! num.contains('-')) {
        println (getValue(num))
        getValue(num)
      }
      if ( num.contains('-') ) {
        val range = num.split('-')
        val low: Int = getValue(range(0))
        val hi: Int = getValue(range(1))
        val newRange = (low to hi).toList
        result = List(result, newRange).flatMap ( x => x )
      }
    }
    result
  }

  val cs = breakUpString(CS)
  val adav = breakUpString(ADAV)
  val mbrc = breakUpString(MBRC)

  def script = {

    val sliceBAM =  swapExt(bam,".bam",".printreads.bam")
    val sliceVCF = swapExt(sliceBAM,".bam",".filtered.vcf")
    add(SliceBAM(bam, sliceBAM))
    callAndEvaluateBAM(sliceBAM, sliceVCF)

    for {
      a <- cs
      b <- adav
      c <- mbrc
    } yield {
    val header = ".CS-" + toStringLength(cs, a) + ".mrav-" + toStringLength(adav, b) + ".MBRC-" + toStringLength(mbrc, c)
    val reduceBAM = swapExt(bam, ".bam", header + ".reduced.bam")
    val reduceVCF = swapExt(reduceBAM,".bam",".filtered.vcf")
    val combineVCF = swapExt(reduceVCF, ".bam", ".filtered.combined.vcf")

    // Generate the new BAMs
    add(ReduceBAM(bam, reduceBAM, a, b, c ))

    // Call SNPs, filter and get variant eval report
    callAndEvaluateBAM(reduceBAM, reduceVCF)

    // Combine the callsets
    add(combine(reduceVCF, sliceVCF, combineVCF))

    // Get a report of the combined callsets
    val eval = new Eval(combineVCF) // evaluate the combined VCF
    eval.select = List("'set==\"Intersection\"'", "'set==\"fullBAM\"'", "'set==\"reducedBAM\"'", "'set==\"filterInreducedBAM-fullBAM\"'", "'set==\"reducedBAM-filterInfullBAM\"'")
    eval.selectName = List("Intersection", "fullBAM", "reducedBAM", "filterInreducedBAM-fullBAM", "reducedBAM-filterInfullBAM")
    add(eval)
    }

  }

  case class ReduceBAM(bam: File, outVCF: File, a: Int, b: Int, c: Int) extends ReduceReads with UNIVERSAL_GATK_ARGS with CoFoJa {
    this.memoryLimit = 3
    this.input_file = List(bam)
    this.o = outVCF
    this.CS = a   // Best safe value 7
    this.ADAV = b     //Best 30
    this.MBRC = c
    this.baq = BAQ.CalculationMode.OFF
    this.intervalsString = List(REDUCE_INTERVAL);
  }

  case class SliceBAM(bam: File, outVCF: File) extends PrintReads with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.input_file = List(bam)
    this.baq = BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    this.o = outVCF
    this.intervalsString = List(REDUCE_INTERVAL);
  }

  case class combine (reducedBAMVCF: File, fullBAMVCF: File, outVCF: File) extends CombineVariants with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("fullBAM", "VCF", fullBAMVCF)
    this.rodBind :+= RodBind("reducedBAM", "VCF", reducedBAMVCF)
    this.rod_priority_list = "reducedBAM,fullBAM"
    this.filteredrecordsmergetype = org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    this.out = outVCF
  }

  def callAndEvaluateBAM(inBAM: File, outVCF: File) {
    val rawVCF   = swapExt(inBAM, ".bam", ".vcf")
    val recalVCF = swapExt(inBAM, ".bam", ".filtered.vcf")
    val recal    = swapExt(inBAM, ".bam", ".recal")
    val tranches = swapExt(inBAM, ".bam", ".tranches")

    add(Call(inBAM, rawVCF),
        HardFilter(rawVCF, recalVCF),
//        VQSR(rawVCF, recal, tranches),
//        applyVQSR(rawVCF, recal, tranches, recalVCF),
        Eval(recalVCF)
//        DiffableTable(rawVCF),   // for convenient diffing
//        DiffableTable(recalVCF)  // for convenient diffing
    )
  }

  case class Eval(@Input vcf: File) extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("eval", "VCF", vcf)
    if ( dbSNP.exists() )
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.doNotUseAllStandardStratifications = true
    this.doNotUseAllStandardModules = true
    this.evalModule = List("TiTvVariantEvaluator", "CountVariants")
    this.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "JexlExpression")
    this.out = swapExt(vcf,".vcf",".eval")
    this.intervalsString = List(REDUCE_INTERVAL);
  }

  case class Call(inBAM: File, outVCF: File) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.input_file = List(inBAM)
    this.stand_call_conf = 50.0
    this.stand_emit_conf = 50.0
    this.dcov = DCOV;
    this.baq = BAQ.CalculationMode.OFF
    this.o = outVCF
    this.intervalsString = List(REDUCE_INTERVAL);

    if ( dbSNP.exists() )
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)

    if ( minimalVCF )
      this.group = List("none")

    if ( REDUCE_INTERVAL != null ) {
      this.intervalsString = List(REDUCE_INTERVAL)
    }
  }

  case class HardFilter (inVCF: File, outVCF: File) extends VariantFiltration with UNIVERSAL_GATK_ARGS {
    this.variantVCF = inVCF
    this.filterName = List("SNP_SB", "SNP_QD", "SNP_HRun")
    this.filterExpression = List("\"SB>=0.10\"", "\"QD<5.0\"", "\"HRun>=4\"")
    this.clusterWindowSize = 10
    this.clusterSize = 3
    this.out = outVCF
  }

  case class VQSR(inVCF: File, outRecal: File, outTranches: File) extends VariantRecalibrator with UNIVERSAL_GATK_ARGS {
    val hapmap_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
    val dbSNP_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_132_b37.leftAligned.vcf"
    val omni_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf"

    this.intervalsString = List(REDUCE_INTERVAL)
    this.rodBind :+= RodBind("input", "VCF", inVCF)
    this.rodBind :+= RodBind("hapmap", "VCF", hapmap_b37, "known=false,training=true,truth=true,prior=15.0")
    this.rodBind :+= RodBind("omni", "VCF", omni_b37, "known=false,training=true,truth=false,prior=12.0")
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP_b37, "known=true,training=false,truth=false,prior=4.0")
    this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "DP")
    this.tranches_file = outTranches
    this.recal_file = outRecal
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
    this.analysisName = outTranches + "_VQSR"
    this.jobName =  outTranches + ".VQSR"
    this.mG = 5
    this.std = 14
    this.percentBad = 0.04
    this.nt = 5
}

  // 4.) Apply the recalibration table to the appropriate tranches
  case class applyVQSR (inVCF: File, inRecal: File, inTranches: File, outVCF: File) extends ApplyRecalibration with UNIVERSAL_GATK_ARGS {
    this.intervalsString = List(REDUCE_INTERVAL)
    this.rodBind :+= RodBind("input", "VCF", inVCF )
    this.tranches_file = inTranches
    this.recal_file = inRecal
    this.ts_filter_level = 99.0
    this.out = outVCF
    this.analysisName = outVCF + "_AVQSR"
    this.jobName =  outVCF + ".applyVQSR"
  }


  case class DiffableTable(@Input vcf: File) extends CommandLineFunction {
    @Output var out: File = swapExt(vcf,".vcf",".table")
    def commandLine = "cut -f 1,2,4,5,7 %s > %s".format(vcf, out)
  }
}

