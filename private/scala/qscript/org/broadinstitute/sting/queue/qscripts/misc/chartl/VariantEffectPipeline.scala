import org.broadinstitute.sting.queue.extensions.gatk.ByTranscriptEvaluator
import org.broadinstitute.sting.queue.library.ipf.vcf.VCFExtractIntervals
import org.broadinstitute.sting.queue.QScript

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 6/28/12
 * Time: 2:24 PM
 * To change this template use File | Settings | File Templates.
 */

class VariantEffectPipeline extends QScript {

  @Input(doc="The raw (un-annotated) VCF file",required=true,shortName="vcf",fullName="vcfToAnalyse")
  var rawInputVCF : File = _

  @Argument(doc="Project name",required=true,shortName="pname",fullName="projectName")
  var projectName : String = _

  @Argument(doc="Predict temporary dir",required=false,fullName="veptmpdir")
  var veptmpdir : String = "/broad/hptmp/chartl/vep/"

  override def script() {
    var ext : String = null
    if ( rawInputVCF.getAbsolutePath.endsWith(".vcf")) {
      ext = ".vcf"
    } else if ( rawInputVCF.getAbsolutePath.endsWith(".vcf.gz") ){
      ext = ".vcf.gz"
    } else {
      ext = ".bcf"
    }

    // compute the interval list
    var extract : VCFExtractIntervals = new VCFExtractIntervals(rawInputVCF)
    add(extract)
    // generate the predictions
    var predictRefseq : PVE = new PVE()
    predictRefseq.memoryLimit = Some(4)
    predictRefseq.variants = rawInputVCF
    predictRefseq.intervals = extract.listOut
    predictRefseq.tempdir = veptmpdir
    predictRefseq.annotVCF = swapExt(rawInputVCF,ext,".predicted.otherDB.vcf")
    predictRefseq.vCommand = "source /humgen/gsa-hphome1/chartl/projects/variantEffect/resources/runPredictorRefseqOnline.src | grep -v \"0+0k\" | sort -nk2,2 | perl /humgen/gsa-hphome1/chartl/sting/public/perl/sortByRef.pl - /humgen/1kg/reference/human_g1k_v37.fasta.fai | python /humgen/gsa-hphome1/chartl/projects/variantEffect/resources/headerToTop.py"
    //add(predictRefseq)

    var predictStandard :PVE = new PVE()
    predictStandard.memoryLimit=Some(4)
    predictStandard.variants = rawInputVCF
    predictStandard.intervals = extract.listOut
    predictStandard.tempdir = veptmpdir
    predictStandard.annotVCF = swapExt(rawInputVCF,ext,".predicted.standard.vcf")
    predictStandard.vCommand = "source /humgen/gsa-hphome1/chartl/projects/variantEffect/resources/runPredictor.src | grep -v \"0+0k\" | sort -nk2,2 | perl /humgen/gsa-hphome1/chartl/sting/public/perl/sortByRef.pl - /humgen/1kg/reference/human_g1k_v37.fasta.fai | python /humgen/gsa-hphome1/chartl/projects/variantEffect/resources/headerToTop.py"
    add(predictStandard)

    val stdEvals : List[ByTranscriptEvaluator] = generateEvals(predictStandard.annotVCF)
    val refEvals : List[ByTranscriptEvaluator] = generateEvals(predictRefseq.annotVCF)

    val stdParsed : List[ParseInfo] = stdEvals.map( u => {
      var pi = new ParseInfo
      pi.bteval = u.out
      pi.perGeneAggroFile = new File(u.out.getAbsolutePath+".genes.txt")
      pi.geneComposite = new File(u.out.getAbsolutePath+".worst.txt")
      pi
    })

    //addAll(stdParsed)

    val refParse : List[ParseInfo] = refEvals.map(u => {
      var pi = new ParseInfo
      pi.bteval = u.out
      pi.perGeneAggroFile = new File(u.out.getAbsolutePath+".genes.txt")
      pi.geneComposite = new File(u.out.getAbsolutePath+".worst.txt")
      pi
    })

    //addAll(refParse)

  }

  def generateEvals(inVCF : File ) : List[ByTranscriptEvaluator] = {
  // evaluate the VCF
    trait StdArgs extends ByTranscriptEvaluator {
      this.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
      this.eval = inVCF
      this.memoryLimit = Some(4)
    }

    var all_all : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    all_all.out = swapExt(inVCF,".vcf",".allTranscripts.allFrequency.eval")
    var all_01 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    all_01.out = swapExt(inVCF,".vcf",".allTranscripts.minAAF_0_01.eval")
    var all_05 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    all_05.out = swapExt(inVCF,".vcf",".allTranscripts.minAAF_0_05.eval")
    //add(all_all,all_01,all_05)
    var ccds_all : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    ccds_all.out = swapExt(inVCF,".vcf",".ccds.allFrequency.eval")
    var ccds_01 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    ccds_01.out = swapExt(inVCF,".vcf",".ccds.minAAF_0_01.eval")
    var ccds_05 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    ccds_05.out = swapExt(inVCF,".vcf",".ccds.minAAF_0_05.eval")
    //add(ccds_all,ccds_01,ccds_05)
    var refseq_all : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    refseq_all.out = swapExt(inVCF,".vcf",".refseq.allFrequency.eval")
    var refseq_01 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    refseq_01.out = swapExt(inVCF,".vcf",".refseq.minAAF_0_01.eval")
    var refseq_05 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    refseq_05.out = swapExt(inVCF,".vcf",".refseq.minAAF_0_05.eval")
    //add(refseq_all,refseq_01,refseq_05)

    return List(all_all,all_01,all_05,ccds_all,ccds_01,ccds_05,refseq_all,refseq_01,refseq_05)
  }

  class PVE extends CommandLineFunction {
    @Input(doc="Queue jar file")
    var queueJar : File = new File("/humgen/gsa-hphome1/chartl/sting/dist/Queue.jar")

    @Input(doc="The input VCF file. Preferably sites-only.")
    var variants : File = _

    @Input(doc="Intervals file for the VCF")
    var intervals : File = _

    @Argument(doc="The number of variants to scatter into")
    var numVariantsPerScatter : Int = 10000

    @Argument(doc="temp dir")
    var tempdir : String = "./"

    @Argument(doc="javaTempDir")
    var javatmp = "/broad/hptmp/chartl/tmp/"

    @Output(doc="annotated vcf")
    var annotVCF : File = _

    @Argument(doc="The script for running Variant Effect Predictor. Must read from stdin and write to stdout.")
    var vCommand : String = _

    def commandLine : String = {
      val scriptPath : String = "/humgen/gsa-hphome1/chartl/sting/private/scala/qscript/org/broadinstitute/sting/queue/qscripts/misc/chartl/PredictVariantEffects.scala"
      "java -Xmx4g -Djava.io.tmpdir=%s -jar %s -S %s -V %s --numVariantsPerScatter %d -L %s --tempdir %s --vScript \"%s\" -out %s -run -bsub -jobQueue week".format(javatmp.getAbsolutePath,
        queueJar.getAbsolutePath,scriptPath,variants.getAbsolutePath,numVariantsPerScatter,
        intervals.getAbsolutePath,tempdir, vCommand, annotVCF.getAbsolutePath)
    }
  }

  class ParseInfo extends CommandLineFunction {
    @Input(doc="The eval file")
    var bteval : File = _

    @Output(doc="The per-gene aggregation")
    var perGeneAggroFile : File = _

    @Output(doc="The worst composite file")
    var geneComposite : File = _

    def commandLine : String = {
      "python /humgen/gsa-hphome1/chartl/projects/variantEffect/resources/runRefseq/merged %s %s %s".format(bteval.getAbsolutePath,
      perGeneAggroFile.getAbsolutePath,geneComposite.getAbsolutePath)
    }
  }


}

