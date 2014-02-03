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

package org.broadinstitute.sting.queue.qscripts.calling

import org.broadinstitute.sting.commandline.Argument._
import org.broadinstitute.sting.commandline.Input._
import org.broadinstitute.sting.pipeline.Pipeline
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.commandline.Input._
import org.broadinstitute.sting.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk._
import java.io.File
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk.RodBind._
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection
;
/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 8/4/11
 * Time: 11:04 AM
 * To change this template use File | Settings | File Templates.
 */

class WholeGenomeIndelCalling extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=false)
  var gatkJar: File = new File("/humgen/gsa-scr1/delangel/GATK/Sting_unstable/dist/GenomeAnalysisTK.jar")

  @Input(doc="output path", shortName="outputDir", required=true)
  var outputDir: String = _

  @Input(doc="run name", shortName="runName", required=true)
  var runName: String = _

  @Input(doc="queue", shortName="queue", required=false)
  var jobQueue: String = "hour"

  @Input(doc="Do only one chromosome", shortName="onlyOneChr", required=false)
  var onlyOneChr: Boolean = false

  @Input(doc="the chromosome to process", shortName="chrToProcess", required=false)
  var chrToProcess: Int = 20

  @Input(doc="scatter count", shortName="scatterCount", required=false)
  var scatterCount: Int = 50


  @Argument(shortName = "R", doc="B37 reference sequence: defaults to broad standard location", required=false)
  var reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Input(doc="path to tmp space for storing intermediate bam files", shortName="outputTmpDir", required=false)
  var outputTmpDir: String = "/broad/shptmp/delangel"

  @Input(doc="BAM file or list of bam files to call", shortName="bam", required=true)
  var bamList: String = _

  @Argument(shortName = "truth", doc="VQSR truth file", required=false)
  var truthFile: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/GoldStandardIndel/gold.standard.indel.MillsAnd1000G.b37.vcf"  )

  @Argument(shortName = "training", doc="VQSR training file", required=false)
  var trainingFile: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/GoldStandardIndel/gold.standard.indel.MillsAnd1000G.b37.vcf"  )

  @Argument(shortName = "intervals", doc="intervals", required=false)
  val myIntervals: String = null;

  @Argument(shortName = "doOnlyVQSR", doc="doOnlyVQSR", required=false)
  val doOnlyVQSR: Boolean = false;

  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560)
  //  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,3000000,48129895,51304566,155270560)
  val COMP_MULLIKIN =  "/humgen/gsa-hpprojects/GATK/data/Comparisons/StandardForEvaluation/b37/comps/NA12878.DIPline.NQScm.expanded.chr20.b37.minReads_2_or_gt2bp.vcf"
  val COMP_MILLS =  "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Mills_Devine_Indels_2011/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.vcf"
  val COMP_MILLS_CHIP =  "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Mills_Devine_Indels_2011/indel_genotype_data_hg19_annotated_polymorphic.vcf"

  private var pipeline: Pipeline = _
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = Some(2)
    this.jobQueue = qscript.jobQueue
    if (qscript.myIntervals != null) {
      this.intervalsString = List(myIntervals);
    }
  }

  def script = {
    var projectBase:String = qscript.outputDir + qscript.runName

    var rawVCFIndels = new File(projectBase + ".raw.vcf")
    val callIndels = new UnifiedGenotyper with CommandLineGATKArgs
     callIndels.out = rawVCFIndels
    // callIndels.dcov = 50
     callIndels.stand_call_conf = 4.0
     callIndels.stand_emit_conf = 4.0
     callIndels.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
  //   callIndels.jobName = qscript.outputTmpDir + "/calls/" + qscript.runName
     callIndels.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
     callIndels.dbsnp =  qscript.dbSNP
     callIndels.sites_only = false
     callIndels.scatterCount = qscript.scatterCount
     callIndels.input_file :+= qscript.bamList
//     callIndels.nt=Some(qscript.nt)

    if (!qscript.doOnlyVQSR) {
      add(callIndels)
    }

    val vr = new VariantRecalibrator with CommandLineGATKArgs
    vr.input :+= callIndels.out
    vr.truth :+= new TaggedFile( qscript.truthFile,"prior=15.0")
    vr.training :+=  new TaggedFile( qscript.trainingFile,"prior=12.0")
    vr.known :+= new TaggedFile(qscript.dbSNP,"prior=3.0")
    //vr.trustAllPolymorphic = true
    vr.mode = VariantRecalibratorArgumentCollection.Mode.INDEL

//    vr.use_annotation = List("QD", "HaplotypeScore",  "ReadPosRankSum","FS","InbreedingCoeff")
   // todo - InbreedingCoeff not appropriate for single sample calling, should be extended to command line argument
    vr.use_annotation = List("QD", "HaplotypeScore",  "ReadPosRankSum","FS")
    vr.TStranche = List("99.0", "97.0","95.0")
    vr.tranches_file = projectBase + ".tranches"
    vr.recal_file = projectBase + ".recal"
    vr.rscriptFile = projectBase + ".plots.R"
    vr.jobOutputFile = vr.recal_file + ".out"
    vr.memoryLimit = 32
    vr.nt = 16
    vr.mG = 8
    vr.percentBad = 0.03
    //vr.std = 14
    add(vr)

    for (tranche <- vr.TStranche) {
      val ar = new ApplyRecalibration with CommandLineGATKArgs
      ar.input :+= (callIndels.out)
      ar.tranches_file = vr.tranches_file
      ar.recal_file = vr.recal_file
      ar.ts_filter_level = tranche.toDouble
      ar.out = projectBase + ".recalibrated." + tranche + ".vcf"
      ar.jobOutputFile = ar.out + ".out"
      ar.memoryLimit = 32
      ar.mode = VariantRecalibratorArgumentCollection.Mode.INDEL
      add(ar)

      val eval = new VariantEval with CommandLineGATKArgs
      eval.tranchesFile = vr.tranches_file
      eval.eval :+= ( ar.out)
      eval.dbsnp = qscript.dbSNP
      eval.doNotUseAllStandardStratifications = true
      eval.doNotUseAllStandardModules = true
      eval.evalModule = List("CountVariants","CompOverlap")
      eval.stratificationModule = List("EvalRod", "CompRod", "Novelty","Sample")
      eval.out = swapExt(ar.out, ".vcf", ".eval")
      eval.jobOutputFile = eval.out + ".out"
      eval.comp :+= new TaggedFile(COMP_MILLS,"Mills")
      eval.comp :+= new TaggedFile(COMP_MILLS_CHIP,"MillsChip")
      eval.comp :+= new TaggedFile(COMP_MULLIKIN,"Mullikin")

      eval.memoryLimit = 32
      eval.nt = 16
      add(eval)
    }

  }

}