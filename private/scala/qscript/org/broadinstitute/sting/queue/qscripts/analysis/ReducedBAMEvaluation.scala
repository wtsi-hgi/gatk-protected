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

  //@Argument(shortName = "callingIntervals", doc="intervals", required=false)
  val CALLING_INTERVAL: String = REDUCE_INTERVAL;

  @Argument(shortName = "dcov", doc="dcov", required=false)
  val DCOV: Int = 250;

  @Argument(shortName = "minimalVCF", doc="", required=false)
  val minimalVCF: Boolean = false;

  @Argument(shortName = "CS", doc = "ContextSize", required = false)
  val CS: String = "6-6";

  @Argument(shortName = "MRAV", doc = "MaximumReadsAtVariableSite", required = false)
  val MRAV: String = "30-30";

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

  //Convert CS, MRAV, MBRC to lists
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
  val mrav = breakUpString(MRAV)
  val mbrc = breakUpString(MBRC)

  def script = {

    val sliceBAM =  swapExt(bam,".bam",".printreads.bam")
    val sliceVCF = swapExt(sliceBAM,".bam",".filtered.vcf")
    add(SliceBAM(bam, sliceBAM))
    callAndEvaluateBAM(sliceBAM, sliceVCF)

    for {
      a <- cs
      b <- mrav
      c <- mbrc
    } yield {
    val header = ".CS-" + toStringLength(cs, a) + ".mrav-" + toStringLength(mrav, b) + ".mbrc-" + toStringLength(mbrc, c)
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
    this.CS = a   // Best value 5
    this.mravs = b     //Best 30
    this.mbrc = c
    this.baq = BAQ.CalculationMode.OFF

    if ( REDUCE_INTERVAL != null )
      this.intervalsString = List(REDUCE_INTERVAL);
  }

  case class SliceBAM(bam: File, outVCF: File) extends PrintReads with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.input_file = List(bam)
    this.baq = BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    this.o = outVCF
    if ( REDUCE_INTERVAL != null )
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
    val vcf = swapExt(inBAM, ".bam", ".vcf")

    val filterSNPs = new VariantFiltration with UNIVERSAL_GATK_ARGS
    filterSNPs.variantVCF = vcf
    filterSNPs.filterName = List("SNP_SB", "SNP_QD", "SNP_HRun")
    filterSNPs.filterExpression = List("\"SB>=0.10\"", "\"QD<5.0\"", "\"HRun>=4\"")
    filterSNPs.clusterWindowSize = 10
    filterSNPs.clusterSize = 3
    filterSNPs.out = outVCF

    add(Call(inBAM, vcf),
        filterSNPs,
        Eval(filterSNPs.out),          // create a variant eval for us
        DiffableTable(vcf),            // for convenient diffing
        DiffableTable(filterSNPs.out)) // for convenient diffing
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
    if ( CALLING_INTERVAL != null )
      this.intervalsString = List(CALLING_INTERVAL);
  }

  case class Call(inBAM: File, outVCF: File) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.input_file = List(inBAM)
    this.stand_call_conf = 50.0
    this.stand_emit_conf = 50.0
    this.dcov = DCOV;
    this.baq = BAQ.CalculationMode.OFF
    this.o = outVCF

    if ( dbSNP.exists() )
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)

    if ( minimalVCF )
      this.group = List("none")

    if ( CALLING_INTERVAL != null ) {
      this.intervalsString = List(CALLING_INTERVAL)
    }
  }

  case class DiffableTable(@Input vcf: File) extends CommandLineFunction {
    @Output var out: File = swapExt(vcf,".vcf",".table")
    def commandLine = "cut -f 1,2,4,5,7 %s > %s".format(vcf, out)
  }
}

