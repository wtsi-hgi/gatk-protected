package org.broadinstitute.sting.queue.qscripts.misc.chartl

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.commandline.{Argument, Output, Input}
import org.broadinstitute.sting.queue.function.{InProcessFunction, CommandLineFunction}
import org.broadinstitute.sting.utils.text.XReadLines
import scala.collection.JavaConversions._
import java.io.{PrintStream, File}
import collection.immutable.HashMap
import org.broadinstitute.sting.utils.exceptions.UserException

/**
 * A script to verify that the null distribution of GCTA is as it should be, by permuting case/control labels.
 *
 * In addition, explore the effect of stratification by permuting case/control labels only within cohorts.
 */

class VerifyNull extends QScript {
  val randgen : scala.util.Random = new scala.util.Random

  @Input(fullName = "phenotype",shortName="phe",doc="The input phenotype file",required=true)
  var phenotypeFile : File = _

  @Input(fullName = "binaryPed",shortName="bped",doc="The input bed file",required=true)
  var genotypeFile : File = _

  @Input(fullName = "covariate",shortName="cov",doc="The input covariate file",required=false)
  var covariateFile : File = _

  @Input(fullName = "quantCovariate", shortName="qcov",doc="The input quantitative covariate file",required=false)
  var quantCovariateFile : File = _

  @Argument(fullName = "cohortOffset",shortName="co",doc="The offest of the cohort column in the covariate file",required=false)
  var cohortOffset : Int = 2

  @Argument(fullName = "nPermutations",doc="The number of permutations of the phenotype file to use",required=false)
  var nPermutations : Int = 20

  override def script = {
    // switches here for multiple analysis
    //scriptNullTest
    scriptStratifiedNull
  }

  def scriptNullTest = {
    // want to calculate the GRM, then permute the phenotype file a bunch of times, and then calculate associations
    var grmCalc : CalculateRelationshipMatrix = new CalculateRelationshipMatrix(genotypeFile,genotypeFile.getName.stripSuffix(".bed"))
    add(grmCalc)
    (new Range(0,nPermutations,1)).foreach( p => {
      // permute
      var perm : PermutePhenotypeLabels = new PermutePhenotypeLabels(phenotypeFile,phenotypeFile.getName+".perm%d.txt".format(p))
      add(perm)
      var varEst : EstimateVariance = new EstimateVariance(grmCalc.grm_gz,perm.phenoOut,genotypeFile.getName.stripSuffix(".bed")+".null%d".format(p))
      add(varEst)
    })
  }

  def scriptStratifiedNull = {
    if ( null != covariateFile || null != quantCovariateFile ) {
      // calculate if necessary
      var grmCalc : CalculateRelationshipMatrix = new CalculateRelationshipMatrix(genotypeFile,genotypeFile.getName.stripSuffix(".bed"))
      add(grmCalc)
      (new Range(0,nPermutations,1)).foreach( p => {
        // permute
        var perm : PermutePhenotypeLabels = new PermutePhenotypeLabels(phenotypeFile,phenotypeFile.getName+".strat.perm%d.txt".format(p),covariateFile,cohortOffset)
        add(perm)
        var varEst : EstimateVariance = new EstimateVariance(grmCalc.grm_gz,perm.phenoOut,quantCovariateFile,covariateFile,genotypeFile.getName.stripSuffix(".bed")+".strat.null%d".format(p))
        add(varEst)
      })
    }
  }


  // this class is for permuting phenotype labels: e.g. the third column of the phenotype file
  // private for access to RNG
  class PermutePhenotypeLabels(inPheno:File,outPheno:File,covar:File,covarOffset:Int) extends InProcessFunction {
    def this(ip:File,op:File) = this(ip,op,null,2)
    def this(ip:File,op:String) = this(ip,new File(op))

    @Input(doc="Unmodified phenotype file")
    var pheno : File = inPheno

    @Input(doc="Covariate file (for permuting within cohorts)",required=false)
    var cohort : File = covar

    @Argument(doc="Offset for particular covariate in the cohort file")
    var cohortOffset : Int = covarOffset

    @Output(doc="The permuted output phenotype file")
    var phenoOut : File = outPheno

    @Argument(doc="Permute within cohorts?")
    var permuteInCohort = false;

    override def run() {
      // simply read the file twice
      if ( null == cohort || ! permuteInCohort ) {
        var pheList : List[String] = asScalaIterator(new XReadLines(pheno)).map(u => u.split("\t")(2)).toList
        val permIter = randgen.shuffle(pheList).iterator
        val out : PrintStream = new PrintStream(phenoOut)
        asScalaIterator(new XReadLines(pheno)).foreach( u => {
          val t = u.split("\t")
          out.printf("%s\t%s\t%s%n",t(0),t(1),permIter.next())
        })
        out.close()
      } else {
        // first: read in the cohort file
        var offset : Int = 2
        if ( null != cohortOffset )
          offset = cohortOffset
        val samCohortMap = asScalaIterator(new XReadLines(cohort)).foldLeft(new HashMap[String,String])( (u: HashMap[String,String],v:String) => {
          val sp = v.split("\t")
          val s : String = sp(1)
          val c : String = sp(offset)
          val d: HashMap[String,String] = u + (s -> c)
          d
        })
        // now read in the pheno file
        val cohortPhenoMap = asScalaIterator(new XReadLines(pheno)).foldLeft(new HashMap[String,List[String]])( (u: HashMap[String,List[String]],v:String) => {
          val sp = v.split("\t")
          val s : String = sp(1)
          val c : String = samCohortMap(s)
          val p : String = sp(2)
          var d : HashMap[String,List[String]] = null
          if ( u.contains(c) ) {
            d = u.updated(c, (u.get(c) ++ List(p)).map(_.toString).toList )
          } else {
            d = u + (c -> List(p))
          }
          d
        })
        // generate iterators for each cohort
        val cohortIterMap = cohortPhenoMap.mapValues( randgen.shuffle(_).iterator )
        // dump to output
        val out : PrintStream = new PrintStream(phenoOut)
        asScalaIterator(new XReadLines(pheno)).foreach( u => {
          val t = u.split("\t")
          val phe = cohortIterMap(samCohortMap(t(1))).next()
          out.printf("%s\t%s\t%s%n",t(0),t(1),phe)
        })
        out.close()
      }
    }
  }

}

class GCTA(bin : File) extends CommandLineFunction {

  val GCTA_STANDARD_FILE : File = new File("/humgen/gsa-hphome1/chartl/projects/t2d/gcta/resources/bin/gcta64")

  @Input(doc="The GCTA binary")
  var binary : File = bin

  def this() = this(new File("/humgen/gsa-hphome1/chartl/projects/t2d/gcta/resources/bin/gcta64"))
  def this(binPath: String) = this(new File(binPath))

  def commandLine : String = {
    // just run "gcta64" which generates help text
    binary.getAbsolutePath
  }
}

class CalculateRelationshipMatrix(binary:File,inputBed:File,outputBase:String) extends GCTA(binary) {
  def this(inBed:File,outBase:String) = this(new File("/humgen/gsa-hphome1/chartl/projects/t2d/gcta/resources/bin/gcta64"),inBed,outBase)

  @Input(doc="The bed file. Will strip \".bed\" from filename")
  var bedFile : File = inputBed

  @Output(doc="The output GRM file. Note that these values are filled in from the output base argument. "+
              "Will dump to the working directory.")
  var grm_gz : File = new File(outputBase+".grm.gz")

  @Output(doc="The output GRM id file. Note that these values are filled in from the output base argument. "+
              "Will dump to the working directory.")
  var grm_id : File = new File(outputBase+".grm .id")

  @Argument(doc="The minimum minor allele frequency.")
  var min_maf : Double = 0.0

  @Argument(doc="The maximum minor allele frequency.")
  var max_maf : Double = 0.5

  @Argument(doc="Chromosomes to inlcude.",required=false)
  var chromosomes : List[String] = Nil

  override def commandLine : String = {
    // build up the chromosome string
    val chrArg : String = chromosomes.foldLeft("")( (u,z) => {
      u + " --chr %s".format(z)
    })
    val bedFileNameBase : String = bedFile.getAbsolutePath.stripSuffix(".bed")
    super.commandLine + " --make-grm --bfile %s --out %s --maf %f --max-maf %f%s".format(bedFileNameBase,
       outputBase,min_maf,max_maf,chrArg)
  }

}

class EstimateVariance(binary:File,inGRM:File,phenoFile:File,qcovarFile:File,covarFile:File,outBase:String) extends GCTA(binary) {
 def this(inGRM:File,phenoFile:File,outBase:String) = this(new File("/humgen/gsa-hphome1/chartl/projects/t2d/gcta/resources/bin/gcta64"),inGRM,phenoFile,null,null,outBase)
  def this(inGRM:File,phenoFile:File,qcovar:File,covar:File,ob:String) = this(new File("/humgen/gsa-hphome1/chartl/projects/t2d/gcta/resources/bin/gcta64"),inGRM,phenoFile,qcovar,covar,ob)

  @Input(doc="A genetic relationship matrix. Can be a file specifying multiple GRMs.")
  var grm : File = inGRM

  @Input(doc="A file listing family id, individual id, and phenotype",required=true)
  var pheno : File = phenoFile

  @Input(doc="A file listing family id, individual id, and quantitative covariates",required=false)
  var qCovar : File = qcovarFile

  @Input(doc="A file listing family id, individual id, and categorical covariates",required=false)
  var covar : File = covarFile

  @Output(doc="The output reml result")
  var reml : File = new File(outBase+".hsq")

  override def commandLine : String = {
    var grmStr = grm.getAbsolutePath.stripSuffix(".grm.gz")
    var phst = pheno.getAbsolutePath
    var cmd : String = super.commandLine + " --reml --grm %s --pheno %s --out %s".format(grmStr,phst,outBase)
    if ( null != qCovar ) {
      cmd += " --qcovar %s".format(qCovar.getAbsolutePath)
    }
    if ( null != covar ) {
      cmd += " --covar %s".format(covar.getAbsolutePath)
    }

    cmd
  }
}


