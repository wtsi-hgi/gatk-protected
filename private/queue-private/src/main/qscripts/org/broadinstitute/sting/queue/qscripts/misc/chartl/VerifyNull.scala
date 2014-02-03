/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

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


