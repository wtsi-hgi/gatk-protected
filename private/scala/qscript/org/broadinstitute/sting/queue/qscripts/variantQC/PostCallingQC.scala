package oneoffs.depristo

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction
import org.broadinstitute.sting.queue.extensions.gatk._

class PostCallingQC extends QScript {
  qscript =>

  @Argument(doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("/home/radon01/depristo/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar")

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  // todo  -- should accept separate indel and snp vcf's, right now script will assume they're combined in one
  @Argument(shortName = "eval", doc="VCFs to evaluate", required=true)
  var evalVCFs: List[File]= _

  @Argument(shortName = "intervals", doc="intervals", required=true)
  val myIntervals: String = null;

  @Argument(shortName = "RPath", doc="RPath", required=false)
  var RPath: File = new File("../R")

  @Argument(shortName = "dbSNP", doc="dbSNP", required=false)
  val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_132.b37.vcf")

  @Argument(shortName = "nt", doc="nt", required=false)
  val nt: Int = 1;

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    this.reference_sequence = referenceFile;
    this.memoryLimit = 4
    this.jobQueue = "gsa"
  }

  // TODO -- should include "standard" eval for plotting expectations

  def script = {
    for ( evalVCF <- evalVCFs ) {
      // The basic summary eval
      createEval(evalVCF, ".summary",
        List("TiTvVariantEvaluator", "CountVariants", "CompOverlap"),
        List("FunctionalClass"))

      // The basic summary eval, by AF
      createEval(evalVCF, ".byAC",
        List("TiTvVariantEvaluator", "CountVariants", "CompOverlap"),
        List("AlleleCount"))

      // By sample
      createEval(evalVCF, ".bySample",
        List("TiTvVariantEvaluator", "CountVariants", "CompOverlap"),
        List("Sample"))
      add(new QCRScript(evalVCF))
    }
  }

  def createEval(evalVCF: File, prefix: String, evalModules: List[String], extraStrats: List[String]) {
    val eval = new Eval(evalVCF)
    eval.out = swapExt(evalVCF,".vcf", prefix + ".eval")
    eval.evalModule = evalModules
    eval.stratificationModule = List("EvalRod", "CompRod", "Novelty") ::: extraStrats
    eval.nt = qscript.nt
    add(eval)
  }

  class Eval(@Input vcf: File) extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.eval :+= TaggedFile(vcf,":VCF")
    if ( dbSNP.exists() )
      this.dbsnp = dbSNP
    this.doNotUseAllStandardStratifications = true
    this.doNotUseAllStandardModules = true
    this.intervalsString = List(myIntervals);
  }

  class QCRScript(vcf: File) extends CommandLineFunction {
    @Output var pdf: File = swapExt(vcf,".vcf", ".pdf")
    val root = swapExt(vcf,".vcf", "") // remove the prefix
    def commandLine = "Rscript %s/variantCallQC.R %s %s %s".format(RPath, root, root, pdf)
  }
}

