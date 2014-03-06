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

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.walkers.genotyper.{UnifiedGenotyperEngine, GenotypeLikelihoodsCalculationModel}
import org.broadinstitute.sting.gatk.downsampling.DownsampleType

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 6/27/12
 * Time: 3:46 PM
 * To change this template use File | Settings | File Templates.
 */

class LargeScaleValidationCallingSingle extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Argument(doc="intervals to process", shortName="intervals", required=false)
  var intervals: String = ""

  @Argument(doc="output path", shortName="outputDir", required=true)
  var outputDir: String =  _

  @Input(doc="input BAM list", shortName="bamList", required=true)
  var bamList: File = _

  @Argument(doc="base output filename", shortName="baseName", required=false)
  var baseName: String = ""

  @Argument(doc="scatterCount", shortName="sc", required=false)
  var variantCallerScatterCount: Int = 1

  @Argument(doc="chromosomes in pool", shortName="ploidy", required=false)
  var ploidy: Int = 24

  @Argument(doc="doRef", shortName="doRef", required=false)
  var doRefSample: Boolean = false

  @Argument(fullName="standard_min_confidence_threshold_for_emitting_and_calling", shortName="stand_conf", doc="The minimum phred-scaled confidence threshold at which variants should be emitted and called", required=false)
  var stand_conf: Double = 30.0

  @Argument(doc="validation set", shortName="vs", required=true)
  var vs: String = _


  private val tmpDir: File = new File("/broad/hptmp/farjoun/tmp/")
  private val reference: File = new File("/humgen/gsa-hpprojects/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.vcf")

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = 2
    this.jobTempDir = qscript.tmpDir
    this.jobQueue = "gsa"
    if (!qscript.intervals.isEmpty) this.intervalsString :+= qscript.intervals


  }

  class PPC(val callName: String, val allelesFile: String) extends UnifiedGenotyper with CommandLineGATKArgs {
    this.scatterCount = qscript.variantCallerScatterCount
    this.input_file :+= qscript.bamList
    this.sample_ploidy = Some(qscript.ploidy)
    this.dbsnp = qscript.dbSNP
    this.out = qscript.outputDir + "/"+qscript.baseName + "."+callName+".vcf"
    if (qscript.doRefSample) {
      this.reference_sample_name = "NA12878"
    }
    //      this.referenceCalls = new File("/humgen/gsa-hpprojects/NA12878Collection/callsets/snps/NA12878.HiSeq.WGS.b37.recalibrated.99_5_cut_for_heng.vcf")
    this.gt_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY
    this.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_CONFIDENT_SITES
    this.max_deletion_fraction=.1

    this.intervals :+= new File(allelesFile)
    if(this.intervals.length + this.intervalsString.length > 1) {
      this.isr =  org.broadinstitute.sting.utils.interval.IntervalSetRule.INTERSECTION
      this.logger.info("Setting IntervalSetRule to INTERSECTION")
    }

    //this.alleles = new File(allelesFile)
    this.ignoreLane = true
    this.maxAltAlleles = Some(1)  // memory usage will overflow without this
    this.dt = DownsampleType.NONE

    this.standard_min_confidence_threshold_for_emitting=stand_conf
    this.standard_min_confidence_threshold_for_calling=stand_conf

    
  }
  class SNPPC(callName: String, allelesFile: String) extends PPC(callName, allelesFile) {
    this.glm = GenotypeLikelihoodsCalculationModel.Model.SNP
    this.referenceCalls = new File("/humgen/gsa-hpprojects/NA12878Collection/callsets/snps/NA12878.HiSeq.WGS.b37.recalibrated.99_5_cut_for_heng.vcf")
  }

  class IndelPC( callName: String,  allelesFile: String) extends PPC(callName, allelesFile) {
    this.glm = GenotypeLikelihoodsCalculationModel.Model.INDEL
    this.referenceCalls = new File("/humgen/gsa-scr1/delangel/IndelGoldSet/CEUTrio.HiSeq.WGS.b37_decoy.recal.ts_95.vcf")
  }

  class BothPC( callName: String,  allelesFile: String) extends PPC(callName, allelesFile) {
    this.glm = GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.referenceCalls = new File("/humgen/gsa-scr1/delangel/IndelGoldSet/CEUTrio.HiSeq.WGS.b37_decoy.recal.ts_95.snp_indel_combined.vcf")
  }


  class Eval(evalVCF: File) extends VariantEval with CommandLineGATKArgs {
    this.eval :+= evalVCF
    this.dbsnp = qscript.dbSNP
    this.doNotUseAllStandardModules = true
    this.evalModule = List("TiTvVariantEvaluator", "CountVariants", "CompOverlap", "IndelSummary")
    this.doNotUseAllStandardStratifications = true
    this.stratificationModule = Seq("EvalRod", "CompRod", "Novelty", "TandemRepeat","Filter") //++ extraStrats
    // this.num_threads = qscript.num_threads
    // this.memoryLimit = 8
    this.out = swapExt(evalVCF, ".vcf", ".eval")
  }

  class Filt(inputVCF: File) extends VariantFiltration with CommandLineGATKArgs {
    this.V = inputVCF
    this.filterExpression = Seq("DP<5000","REFDEPTH<500")
    this.filterName = Seq("LowDepth","LowReferenceSampleDepth")
    this.out = swapExt(inputVCF, ".vcf",".filtered.vcf")
  }
  class SampleEval(evalVCF: File) extends Eval(evalVCF) {
    this.stratificationModule :+= "Sample"
    // this.num_threads = qscript.num_threads
    // this.memoryLimit = 8
    this.out = swapExt(evalVCF, ".vcf", ".bySample.eval")
  }

  class ACEval(evalVCF: File) extends Eval(evalVCF) {
    this.stratificationModule :+= "AlleleCount"
    this.out = swapExt(evalVCF, ".vcf", ".byAC.eval")
  }

  def script = {

    val genotyper:PPC= qscript.vs match{
      case "omniMono" => new SNPPC("omniMono","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_OmniMono.SNP.sites.vcf")
      case "omniPoly" => new SNPPC("omniPoly","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_OmniPoly.SNP.sites.vcf")
      case "exomeChipPoly" => new SNPPC("exomeChip","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_ExomeChip.SNP.sites.vcf")
      case "millsPoly" => new IndelPC("millsPoly","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.1000_control_sites_MillsGenotypeInPhase1.INDEL.sites.vcf")
      //case "multiAllelicIndels" => { val temp=new IndelPC("multiallelicIndels","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.2000_validation_sites_multiAllelicIndels.sites.vcf")
      //  temp.maxAltAlleles = Some(3)
      //  temp
      //}
      case "lostToImputation" => new SNPPC("lostToImputation","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.2000_lost_to_Imputation.sites.vcf")
      //case "multiAllelicSNPs" => {val temp=new SNPPC("multiAllelicSNPs","/humgen/gsa-hpprojects/dev/largeScaleValidation/inputSets/triallelics.EricBanks.vcf")
      //   temp.maxAltAlleles = Some(3)
      //}
      case "LoF" => new BothPC("LOF","/humgen/gsa-hpprojects/dev/largeScaleValidation/inputSets/LOF.DanielMacArthur.vcf")
      case "afIndels" => new IndelPC("afIndels","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.5000_validation_sites_AF_distributed.indels.sites.vcf")
      case "unifIndels" => new IndelPC("unifIndels","/humgen/gsa-hpprojects/dev/largeScaleValidation//outputVCFs/ALL.wgs.5000_validation_sites_Uniformly_distributed.indels.sites.vcf")
      case "afSNPs" => new SNPPC("afSNPs","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.8000_validation_sites_AF_distributed.snp.sites.vcf")
      case "unifSNPs" => new SNPPC("unifSNPs","/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/ALL.wgs.8000_validation_sites_Uniformly_distributed.snp.sites.vcf")
      case _=>{ this.logger.error("vs must be one of the allowed validation sites. no match found for \"%s\"".format( qscript.vs)) ; throw new Error("Unknown Set") }
    }
    add(genotyper)


    val filter = new Filt(genotyper.out)
    add(filter)

  }

}

