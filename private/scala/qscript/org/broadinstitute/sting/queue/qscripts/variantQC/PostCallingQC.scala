import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

class PostCallingQC extends QScript {
  qscript =>

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  // todo  -- should accept separate indel and snp vcf's, right now script will assume they're combined in one
  @Argument(shortName = "eval", doc="VCFs to evaluate")
  var evalVCFs: List[File] = Nil

  @Argument(shortName = "L", doc="intervals", required=false)
  val myIntervals: List[String] = Nil

  @Argument(shortName = "XL", doc="exclude intervals", required=false)
  var myExcludeIntervals: List[String] = Nil

  @Argument(shortName = "RPath", doc="RPath", required=false)
  var RPath: File = new File("../R")

  @Argument(shortName = "dbSNP", doc="dbSNP", required=false)
  val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_132.b37.vcf")

  @Argument(shortName = "nt", doc="nt", required=false)
  val nt: Int = 1

  trait UniversalGATKArgs extends CommandLineGATK {
    this.reference_sequence = referenceFile
    this.memoryLimit = 2
  }

  // TODO -- should include "standard" eval for plotting expectations

  def script() {
    for ( evalVCF <- evalVCFs ) {
      // The basic summary eval
      // The basic summary eval, by AF
      val byAC = createEval(evalVCF, ".byAC",
        List("TiTvVariantEvaluator", "CountVariants", "CompOverlap"),
        List("AlleleCount"))

      // The basic summary eval, broken down by sample as well as functional class
      val bySampleByFunctionalClass = createEval(evalVCF, ".bySampleByFunctionalClass",
        List("TiTvVariantEvaluator", "CountVariants", "CompOverlap"),
        List("Sample","FunctionalClass"))

      val evals = List(byAC.out, bySampleByFunctionalClass.out)

      val qc = new QCRScript(evalVCF, evals)
      qc.jobOutputFile = qc.pdf + ".out"
      add(qc)

      add(qc)
    }
  }

  def createEval(evalVCF: File, prefix: String, evalModules: List[String], extraStrats: List[String]) = {
    val eval = new Eval(evalVCF)
    eval.out = swapExt(evalVCF,".vcf", prefix + ".eval")
    eval.evalModule = evalModules
    eval.stratificationModule = List("EvalRod", "CompRod", "Novelty") ::: extraStrats
    eval.nt = qscript.nt
    add(eval)
    eval
  }

  class Eval(@Input var vcf: File) extends VariantEval with UniversalGATKArgs {
    this.eval :+= vcf
    this.dbsnp = dbSNP
    this.doNotUseAllStandardStratifications = true
    this.doNotUseAllStandardModules = true
    this.intervalsString = myIntervals
    this.excludeIntervalsString = myExcludeIntervals
  }

  class QCRScript(@Input var vcf: File, @Input var evals: List[File]) extends CommandLineFunction {
    @Output var pdf: File = swapExt(vcf, ".vcf", ".pdf")
    val root = swapExt(vcf,".vcf", "") // remove the prefix
    def commandLine = "Rscript %s/variantCallQC.R %s %s %s".format(RPath, root, root, pdf)
  }
}
