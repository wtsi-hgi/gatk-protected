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

  @Argument(shortName = "callingIntervals", doc="intervals", required=false)
  val CALLING_INTERVAL: String = null;

  @Argument(shortName = "dcov", doc="dcov", required=false)
  val DCOV: Int = 250;

  @Argument(shortName = "minimalVCF", doc="", required=false)
  val minimalVCF: Boolean = false;

  val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_132.b37.vcf")

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = referenceFile;
    this.memoryLimit = 4
  }

  trait CoFoJa extends JavaCommandLineFunction {
    override def javaOpts = super.javaOpts //  + " -javaagent:lib/cofoja.jar"
  }

  def script = {

    val reduceBAM = swapExt(bam, ".bam", ".reduced.bam")
    val reduceVCF = swapExt(reduceBAM,".bam",".filtered.vcf")
    val sliceBAM =  swapExt(bam,".bam",".printreads.bam")
    val sliceVCF = swapExt(sliceBAM,".bam",".filtered.vcf")
    val combineVCF = swapExt(reduceVCF, ".bam", ".filtered.combined.vcf")

    // Generate the new BAMs
    add(ReduceBAM(bam, reduceBAM),
        SliceBAM(bam, sliceBAM))

    // Call SNPs, filter and get variant eval report
    callAndEvaluateBAM(reduceBAM, reduceVCF)
    callAndEvaluateBAM(sliceBAM, sliceVCF)

    // Combine the callsets
    add(combine(reduceVCF, sliceVCF, combineVCF))

    // Get a report of the combined callsets
    val eval = new Eval(combineVCF) // evaluate the combined VCF
    eval.select = List("'set==\"Intersection\"'", "'set==\"fullBAM\"'", "'set==\"reducedBAM\"'", "'set==\"filterInreducedBAM-fullBAM\"'", "'set==\"reducedBAM-filterInfullBAM\"'")
    eval.selectName = List("Intersection", "fullBAM", "reducedBAM", "filterInreducedBAM-fullBAM", "reducedBAM-filterInfullBAM")
    add(eval)
  }

  case class ReduceBAM(bam: File, outVCF: File) extends ReduceReads with UNIVERSAL_GATK_ARGS with CoFoJa {
    this.memoryLimit = 3
    this.input_file = List(bam)
    this.o = outVCF
    this.CS = 20
    this.mravs = 50
    this.mbrc = 10000
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

