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

import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.IndelStatistics
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.apache.commons.io.FilenameUtils
import scala.Some
;

class Phase1IndelVQSR extends QScript {
  qscript =>
  // todo -- update to released version when things stabilize
  @Argument(shortName = "gatk",doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("/humgen/gsa-scr1/delangel/GATK/Sting_unstable/dist/GenomeAnalysisTK.jar")

  @Argument(shortName = "R", doc="B37 reference sequence: defaults to broad standard location", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")


  @Argument(shortName = "dataDir", doc="Path to the standard evaluation data files", required=false)
  val DATA_DIR = "/humgen/gsa-hpprojects/GATK/data/Comparisons/StandardForEvaluation/b37/"
  @Argument(shortName = "baseDir", doc="Path to the standard evaluation data files", required=false)
  val baseDir = "/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/20110820VQSRConsensus_WG_DP"

  @Argument(shortName = "outDir", doc="Path to the output files", required=false)
  val OUT_DIR = "/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/20110820VQSRConsensus_WG_DP"

  @Argument(shortName = "rawCalls", doc="VQSR raw input file", required=false)
  var rawCalls: File = new File("/humgen/1kg/processing/production_wgs_phase1/consensus_wgs/indel_v1/calls/combined.phase1.chr1.raw.indels.vcf")

  @Argument(shortName = "truth", doc="VQSR truth file", required=false)
  var truthFile: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Mills_Devine_Indels_2011/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.vcf"  )

  @Argument(shortName = "training", doc="VQSR training file", required=false)
  var trainingFile: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Mills_Devine_Indels_2011/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.vcf"  )

  @Argument(shortName = "intervals", doc="intervals", required=false)
  val myIntervals: String = null;


  val populations = List("EUR","AMR","ASN","AFR")

  @Argument(shortName = "evalStandard1000GCalls", doc="If provided, we'll include some standard 1000G data for evaluation", required=false)
  val EVAL_STANDARD_1000G_CALLS: Boolean = true

  @Argument(shortName = "numG", doc="If provided, we'll include some standard 1000G data for evaluation", required=false)
  val numG: Int = 8

  @Argument(shortName = "pctBad", doc="If provided, we'll include some standard 1000G data for evaluation", required=false)
  val pctBad: Double = 0.03

  @Argument(shortName = "runName", doc="Run Name", required=false)
  val runName:String = "mills100"
  val COMPS_DIR = DATA_DIR + "comps/"
  val EVALS_DIR = DATA_DIR + "evals/"

  @Argument(shortName = "createAllPos", doc="If provided, create all POPS file", required=false)
  val CREATE_ALL_POPS_FILE: Boolean = false

  @Argument(shortName = "pops", doc="Populations to do", required=false)
  val moreIndelsToEval: List[String] = List("EUR","ASN","AFR","AMR")


  val VARIANT_TYPES: List[String] = List("indels", "snps")
      /*
  val VARIANT_TYPE_VT: Map[String, List[org.broad.tribble.util.variantcontext.VariantContext.Type]] = Map(
    "indels" -> List(org.broad.tribble.util.variantcontext.VariantContext.Type.INDEL, org.broad.tribble.util.variantcontext.VariantContext.Type.MIXED, org.broad.tribble.util.variantcontext.VariantContext.Type.NO_VARIATION),
    "snps" -> List(org.broad.tribble.util.variantcontext.VariantContext.Type.SNP, org.broad.tribble.util.variantcontext.VariantContext.Type.NO_VARIATION)
  )
         */
  val SITES_DIR: String = "sitesFiles"

  // path to b37 DBSNP
  val MY_DBSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf")
  val dindelCalls: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/AFR+EUR+ASN+1KG.dindel_august_release_merged_pilot1.20110126.sites.vcf"

  val DINDEL:String = "/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/officialCalls/ALL.wgs.dindel.20101123.leftAlignedPassingIndels.sites.vcf"
  val BI:String ="/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/officialCalls/ALL.wg.bi.20101123.leftAlignedPassingIndelsCollapsed.sites.vcf"
  val SI:String ="/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/officialCalls/ALL.SI.20101123.LeftAlignedPassingIndels.sites.vcf"
  val OX:String ="/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/officialCalls/ALL.wg.oxford.20101123.leftAlignedPassingIndelsCollapsed.sites.vcf "
  val BC:String ="/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/officialCalls/ALL.wg.bc.20101123.leftAlignedPassingIndelsCollapsed.sites.vcf"
  val WG2of5:String =  "/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/officialCalls/ALL.wgs.2of5_consensus_v2.20101123.indels.sites.vcf"

  val basePath:String = "/humgen/1kg/processing/production_wgs_phase1/consensus_wgs/indel_v1/calls/"
  val baseNewPath:String = "/humgen/1kg/processing/production_wgs_phase1/consensus_wgs/indel_v3/calls/"
  var COMPS: List[Comp] = Nil
  def addComp(comp: Comp) { COMPS = comp :: COMPS }

  var EVALS: List[Eval] = Nil
  def addEval(eval: Eval) { EVALS = eval :: EVALS }
  def addEvalFromCMD(file: File, t: String) { addEval(new Eval(file.getName, t, file.getName)) }

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJarFile
    this.reference_sequence = qscript.referenceFile
    this.memoryLimit = Some(32)
    // this.rodBind :+= RodBind("dbsnp", "VCF", qscript.dbSNP )
    this.jobQueue = "gsa"

    if (qscript.myIntervals != null)
      this.intervalsString = List(qscript.myIntervals )

  }
  class Comp(val name: String, val evalType: String, val filename: String) {
    val file: File = new File(filename)
  }

  class Eval(val name: String, val evalType: String, val filename: String ) {
    val file: File = new File(filename)
  }

  def initializeStandardDataFiles() = {
    //
    // Standard evaluation files for indels
    //
    //addComp(new Comp("CG.38samples", "indels", COMPS_DIR+"CG.Indels.leftAligned.b37.vcf"))
    addComp(new Comp("g1k.pilot1.validation", "indels", COMPS_DIR+"pilot1_indel_validation_2009.b37.vcf"))
    //addComp(new Comp("NA12878.hand_curated", "indels", "NA12878.validated.curated.polymorphic.indels.vcf"))
    addComp(new Comp("NA12878.Mullikin", "indels", COMPS_DIR+"NA12878.DIPline.NQScm.expanded.chr20.b37.minReads_2_or_gt2bp.vcf"))
    //addComp(new Comp("Mills.25pct", "indels", "/humgen/gsa-scr1/delangel/devine_data/indel_hg19_051711_leftAligned_25percent_chr20.vcf"))
    //addComp(new Comp("Phase1Validation", "indels", "/humgen/gsa-scr1/delangel/VQSRIndels/1KG_Validation_Phase1_SNPs_05032011.HG19.finalized.vcf"))
   addComp(new Comp("Phase1Validation", "indels", "/humgen/gsa-scr1/delangel/VQSRIndels/1000G.20101123.validation_set_v1.QCed.indels.vcf"))
   addComp(new Comp("MillsChipValidationPoly", "indels", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Mills_Devine_Indels_2011/indel_genotype_data_hg19_annotated_polymorphic.vcf"))
   addComp(new Comp("MillsChipValidationMono", "indels", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Mills_Devine_Indels_2011/indel_genotype_data_hg19_annotated_monomorphic.vcf"))
    addComp(new Comp("MillsAll", "indels", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Mills_Devine_Indels_2011/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed.sites.vcf"))

    //
    // INDEL call sets
    //

    if ( EVAL_STANDARD_1000G_CALLS ) {
      addEval(new Eval("dindel", "indels",qscript.DINDEL))
      addEval(new Eval("si", "indels",qscript.SI))
      addEval(new Eval("bi", "indels", qscript.BI))
      addEval(new Eval("bc", "indels", qscript.BC))
      addEval(new Eval("ox", "indels", qscript.OX))
      addEval(new Eval("2of5", "indels",qscript.WG2of5))
      addEval(new Eval("union", "indels", "/humgen/gsa-hpprojects/dev/delangel/Phase1Calls/officialCalls/ALL.wgs.union_v2.20101123.indels.sites.vcf"))
    }

  }

  def script = {

    initializeStandardDataFiles();

    var ts:Double = 0.0
    var tranches =  List("96.0","95.0","94.0","93.0","92.0","91.0")

    var numG:Int = qscript.numG
    var pctBad:Double = qscript.pctBad
    val runName:String = qscript.runName +  "_mG%d_pb%1.2f_QD_FS_HS_RP_IC_DP".format(numG,pctBad)

    val rawCalls = qscript.rawCalls
    var tranchesFile = new File(qscript.baseDir +"/VQSRData/%s.tranches".format(runName))
    var recalFile = new File(qscript.baseDir +"/VQSRData/%s.recal".format(runName))
    var rscriptFile = new File(qscript.baseDir +"/VQSRData/%s.plots.R".format(runName))



    var vr = new VariantRecalibrator with CommandLineGATKArgs

    var chrList =  List(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)

    for(chr <- chrList) {
      var chrStr:String = chr.toString;
      if (chr == 23)
        chrStr = "X"

      vr.rodBind :+= RodBind("input"+chrStr, "VCF", qscript.basePath+"combined.phase1.chr"+chrStr+".raw.indels.vcf")
      vr.rodBind :+= RodBind("inputz"+chrStr, "VCF", qscript.baseNewPath+"combined.phase1.chr"+chrStr+".raw.indels.vcf")

    }

    vr.rodBind :+= RodBind("truth", "VCF",qscript.truthFile,"known=true,training=false,truth=true,prior=15.0" )
    vr.rodBind :+= RodBind("training", "VCF",qscript.trainingFile,"known=true,training=true,truth=false,prior=12.0" )
    vr.rodBind :+= RodBind("dbsnp", "VCF",qscript.MY_DBSNP,"known=true,training=false,truth=false,prior=8.0" )

    vr.rodBind :+= RodBind("BC", "VCF",qscript.BC,"consensus=true" )
    vr.rodBind :+= RodBind("BI", "VCF",qscript.BI,"consensus=true" )
    vr.rodBind :+= RodBind("SI", "VCF",qscript.SI,"consensus=true" )
    vr.rodBind :+= RodBind("DINDEL", "VCF",qscript.DINDEL,"consensus=true" )
    vr.rodBind :+= RodBind("OXFORD", "VCF",qscript.OX,"consensus=true" )

    vr.mode = VariantRecalibratorArgumentCollection.Mode.INDEL
    vr.tranchesFile = tranchesFile
    vr.recalFile = recalFile
    vr.rscriptFile = rscriptFile

    vr.an = List("QD","FS","HaplotypeScore","ReadPosRankSum","InbreedingCoeff","DP")
    vr.maxGaussians = Some(numG)
    vr.tranche = tranches
    vr.nt = Some(12)
    vr.percentBad = Some(pctBad)
    vr.std = Some(14.0)
    vr.jobName = qscript.baseDir +"/tmp/%s.vr".format(runName)
    vr.memoryLimit = Some(32)
    vr.projectConsensus = true
    add(vr)

    val VE = new MyEval()
    val  ve2 = new MyEval

   // VE.VT = VARIANT_TYPE_VT("indels")
    VE.o = new File(OUT_DIR+"/"+ runName + ".eval")
    VE.jobName = qscript.baseDir +"/tmp/"+runName + ".eval"

    for (tas: String <- tranches) {
      ts = tas.toDouble
      val outFile = new File(OUT_DIR+"/calls/phase1.WG.recal_%s_ts_%4.1f.indels.sites.vcf".format(runName,ts))

      var ar = new ApplyRecalibration with CommandLineGATKArgs
      for(chr <- chrList) {
//      for(chr <- List(20)) {
        var chrStr:String = chr.toString;
        if (chr == 23)
          chrStr = "X"

        ar.rodBind :+= RodBind("input"+chrStr, "VCF", qscript.basePath+"combined.phase1.chr"+chrStr+".raw.indels.vcf")
       ar.rodBind :+= RodBind("inputz"+chrStr, "VCF", qscript.baseNewPath+"combined.phase1.chr"+chrStr+".raw.indels.vcf")

      }
      ar.mode = VariantRecalibratorArgumentCollection.Mode.INDEL
      ar.tranchesFile = tranchesFile
      ar.recalFile = recalFile
      ar.ts_filter_level = Some(ts)
      ar.sites_only = true
      ar.o = outFile
      ar.jobName = qscript.baseDir +"/tmp/%s_ts%4.1f.ar".format(runName,ts)

      add(ar)

      VE.rodBind :+= RodBind("eval_ts%4.1f".format(ts), "VCF", ar.o)
      ve2.rodBind :+= RodBind("eval_ts%4.1f".format(ts), "VCF", ar.o)

    }


    // add evals
    for ( calls <- EVALS )
      VE.rodBind :+= RodBind("eval_" + calls.name, "VCF", calls.file)

    // add comps
    for ( comp <- COMPS )
      VE.rodBind :+= RodBind("comp_" + comp.name, "VCF", comp.file)

    VE.jobName = qscript.baseDir +"/tmp/"+runName + ".eval"
    add(VE)



    // comps are now other callsets to measure overlap
    ve2.rodBind :+= RodBind("comp_dindel", "VCF",qscript.DINDEL)
    ve2.rodBind :+= RodBind("comp_bc", "VCF", qscript.BC)
    ve2.rodBind :+= RodBind("comp_bi", "VCF", qscript.BI)
    ve2.rodBind :+= RodBind("comp_ox", "VCF", qscript.OX)
    ve2.rodBind :+=  RodBind("comp_si", "VCF", qscript.SI)
    ve2.rodBind :+= RodBind("comp_2of5", "VCF", qscript.WG2of5)
    //ve2.VT = VARIANT_TYPE_VT("indels")
    ve2.o = new File(OUT_DIR+"/"+ runName + ".comps.eval")
    ve2.jobName = qscript.baseDir +"/tmp/"+runName + ".comps.eval"

    add(ve2)
  }


  /**
   * Base class for VariantEval used here
   */
  class MyEval() extends VariantEval with CommandLineGATKArgs {
    this.noST = true
    this.noEV = true
    this.nt = Some(8)
    this.memoryLimit = Some(32)
    this.evalModule :+= "ValidationReport"
    //this.evalModule :+= "IndelMetricsByAC"
    //this.evalModule :+= "IndelStatistics"
    this.evalModule :+= "CountVariants"
    this.evalModule :+= "CompOverlap"
    //this.evalModule :+= "IndelClasses"
  }



}
