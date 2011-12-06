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
  val num_cpu_threads: Int = 1

  def script() {
    for ( evalVCF <- evalVCFs ) {
      // The basic summary eval
      // The basic summary eval, by AF
      val byAC = new Eval(evalVCF, ".byAC", List("AlleleCount"))
      add(byAC)

      // The basic summary eval, broken down by sample as well as functional class
      val bySampleByFunctionalClass = new Eval(evalVCF, ".bySampleByFunctionalClass", List("Sample","FunctionalClass"))
      add(bySampleByFunctionalClass)

      val qc = new QCRScript(evalVCF, byAC.out, bySampleByFunctionalClass.out)
      add(qc)
    }
  }

  class Eval(evalVCF: File, prefix: String, extraStrats: List[String]) extends VariantEval {
    this.reference_sequence = referenceFile
    this.intervalsString = myIntervals
    this.excludeIntervalsString = myExcludeIntervals
    this.eval :+= evalVCF
    this.dbsnp = dbSNP
    this.doNotUseAllStandardModules = true
    this.evalModule = List("TiTvVariantEvaluator", "CountVariants", "CompOverlap")
    this.doNotUseAllStandardStratifications = true
    this.stratificationModule = List("EvalRod", "CompRod", "Novelty") ::: extraStrats
    this.num_cpu_threads = qscript.num_cpu_threads
    this.memoryLimit = 2
    this.out = swapExt(evalVCF, ".vcf", prefix + ".eval")
    this.jobOutputFile = out + ".out"
  }

  class QCRScript(@Input var vcf: File, @Input var byAC: File, @Input var bySite: File) extends CommandLineFunction {
    @Output var pdf: File = swapExt(vcf, ".vcf", ".pdf")
    private val project = vcf.getName.stripSuffix(".vcf")
    this.jobOutputFile = pdf + ".out"
    def commandLine = "Rscript %s/variantCallQC.R %s %s %s".format(RPath, project, bySite, byAC, pdf)
  }
}
