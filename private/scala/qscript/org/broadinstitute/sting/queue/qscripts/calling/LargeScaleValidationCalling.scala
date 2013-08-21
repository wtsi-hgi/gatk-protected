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

class LargeScaleValidationCalling extends QScript {
  qscript =>

//  @Input(doc="path to GATK jar", shortName="gatk", required=true)
//  var gatkJar: File = _

  @Input(doc="intervals to process", shortName="intervals", required=false)
  var intervals: File = _

  @Argument(doc="output path", shortName="outputDir", required=true)
  var outputDir: String = _

  @Argument(doc="base output filename", shortName="baseName", required=false)
  var baseName: String = ""

  @Argument(doc="scatterCount", shortName="sc", required=false)
  var variantCallerScatterCount: Int = 1

  @Argument(doc="chromosomes in pool", shortName="ploidy", required=false)
  var ploidy: Int = 24

  @Argument(doc="max alt alleles", shortName="altAlleles", required=false)
  var altAlleles: Int = 3


  @Argument(doc="doRef", shortName="doRef", required=false)
  var doRefSample: Boolean = false

  //  @Argument(fullName="standard_min_confidence_threshold_for_emitting_and_calling", shortName="stand_conf", doc="The minimum phred-scaled confidence threshold at which variants should be emitted and called", required=false)
  //  var stand_conf: Option[Double] = 10.0

  @Argument(doc="Use Genotype Given Alleles mode for calling", shortName="doGGA", required=false)
  var doGGA: Boolean = false

  @Argument(doc="validation set", shortName="vs", required=false)
  var vs: String = _

  val baseDir = "/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/"
  val bamListRef: File = new File(baseDir+"dataAnalysis/Pools.bam.list")
  val bamListNoRef: File = new File(baseDir+"dataAnalysis/Pools.noRef.bam.list")

  val originalSites: File = new File(baseDir+"ALL.wgs.allCombinedValidationSites.ACannotated.corr.sites.vcf")

  private val tmpDir: File = new File("/broad/hptmp/delangel/tmp/")
  private val reference: File = new File("/humgen/gsa-hpprojects/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.vcf")
  private val oneKGRelease: File = new File("/humgen/1kg/DCC/ftp/release/20110521/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz")
  private val oneKGReleaseMinusBadPools: File = new File("/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/ALL.wgs.phase1_release_v3.20101123.snps_indels_svs.excludingSamplesInBadPools.sites.vcf")
  private val lofFile: File = new File(baseDir+"outputVCFs/LOF.DanielMacArthur_20130212.fixed.vcf")
  private val exomeChip: File = new File(baseDir+"finalPaperData/1000G.exomechip.20121009.subsetToGoodPoolSamples.sites.vcf")
  private val axiomChip: File = new File(baseDir+"finalPaperData/ALL.wex.axiom.20120206.snps_and_indels.subsetToGoodPoolSamples.sites.vcf")
  private val omniChip: File = new File(baseDir+"finalPaperData/Omni25_genotypes_2141_samples.b37.samplesInGoodPools.vcf")

  trait CommandLineGATKArgs extends CommandLineGATK {
//    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = 2
    this.jobTempDir = qscript.tmpDir
    this.jobQueue = "hour"
    //    if (!qscript.intervals.isEmpty) this.intervals :+= qscript.intervals


  }

  class PPC(val callName: String, val intervalFile: File) extends UnifiedGenotyper with CommandLineGATKArgs {
    this.scatterCount = qscript.variantCallerScatterCount
    this.sample_ploidy = Some(qscript.ploidy)
    this.dbsnp = qscript.dbSNP
    var refStr:String = ""

    if (doRefSample) {
      this.reference_sample_name = "NA12878"
      this.input_file :+= qscript.bamListRef
      refStr = "withRef"
    }
    else {
      this.input_file :+= qscript.bamListNoRef
      refStr = "noRef"

    }
    this.out = qscript.outputDir + "/"+qscript.baseName + "."+callName+"."+refStr+".vcf"

    if (doGGA) {
      this.gt_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
      this.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
      this.alleles = swapExt(baseDir+"outputVCFs/",intervalFile,"interval_list","vcf")

    } else {
      this.gt_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY
      this.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
    }

    this.intervals = Seq(intervalFile)

    this.ignoreLane = true
    if (callName.contains("multiAllelic"))
      this.maxAltAlleles = Some(altAlleles)  // memory usage will overflow without this
    else
      this.maxAltAlleles = 1;

    this.dt = DownsampleType.NONE

    this.standard_min_confidence_threshold_for_emitting = Some(5.0)
    this.standard_min_confidence_threshold_for_calling= Some(5.0)

  }
  class SNPPC(callName: String, intervalFile: File) extends PPC(callName, intervalFile) {
    this.glm = GenotypeLikelihoodsCalculationModel.Model.SNP
//    this.referenceCalls = new File("/humgen/gsa-hpprojects/NA12878Collection/callsets/snps/NA12878.HiSeq.WGS.b37.recalibrated.99_5_cut_for_heng.vcf")
    this.referenceCalls = new File("/humgen/1kg/DCC/ftp/technical/working/20130610_ceu_hc_trio/broad/CEU.wgs.UnifiedGenotyper_bi.20130520.snps_indels.high_coverage_pcr_free.genotypes.vcf.gz")
    this.max_deletion_fraction=.1
  }

  class IndelPC(callName: String, intervalFile: File) extends PPC(callName, intervalFile) {
    this.glm = GenotypeLikelihoodsCalculationModel.Model.INDEL
    this.minIndelFrac = Some(0.01)
//    this.referenceCalls = new File(baseDir+"inputSets/CEUTrio.HiSeq.WGS.b37_decoy.recal.ts_95.vcf")
    this.referenceCalls = new File("/humgen/1kg/DCC/ftp/technical/working/20130610_ceu_hc_trio/broad/CEU.wgs.UnifiedGenotyper_bi.20130520.snps_indels.high_coverage_pcr_free.genotypes.vcf.gz")
    this.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES

  }


  class Eval(evalVCF: File, compVCF: File) extends VariantEval with CommandLineGATKArgs {
    this.eval :+= evalVCF
    this.dbsnp = qscript.dbSNP
    this.doNotUseAllStandardModules = true
    this.evalModule = List("CountVariants", "ValidationReport")
    this.doNotUseAllStandardStratifications = true
    this.stratificationModule = Seq("EvalRod", "Filter") //++ extraStrats
    this.out = swapExt(evalVCF, ".vcf", ".eval")
    this.ploidy = Some(qscript.ploidy)
    this.numSamples = Some(92)

    if (compVCF != null) {
      this.comp :+= compVCF
      this.intervals = Seq(compVCF)
    }
  }

  class Filt(inputVCF: File) extends VariantFiltration with CommandLineGATKArgs {
    this.V = inputVCF
    this.filterExpression = Seq("DP<5000")
    this.filterName = Seq("LowDepth")
    if (qscript.doRefSample) {
      this.filterExpression :+= "REFDEPTH<500"
      this.filterName :+= "LowReferenceSampleDepth"
    }
    this.filterExpression :+= "QUAL<100"
    this.filterName :+= "LowQual"
    this.filterExpression :+= "FS > 200"
    this.filterName :+= "FisherStrand"

    this.out = swapExt(inputVCF, ".vcf",".filtered.vcf")
  }
  class SampleEval(evalVCF: File, compVCF: File) extends Eval(evalVCF, compVCF) {
    this.stratificationModule :+= "Sample"
    // this.num_threads = qscript.num_threads
    // this.memoryLimit = 8
    this.out = swapExt(evalVCF, ".vcf", ".bySample.eval")
  }

  class ACEval(evalVCF: File, compVCF: File) extends Eval(evalVCF, compVCF) {
    this.stratificationModule :+= "AlleleCount"
    this.out = swapExt(compVCF, ".vcf", ".byAC.eval")
  }

  class Annot(inputVCF: File) extends VariantAnnotator with CommandLineGATKArgs {
    this.V = inputVCF
    this.intervals :+= inputVCF
    this.out = swapExt(inputVCF, ".vcf",".annotated.vcf")
    this.resource :+= qscript.originalSites
    this.E = Seq("resource.set","resource.ALT" )

  }

  class OneKGAnnot(inputVCF: File) extends VariantAnnotator with CommandLineGATKArgs {
    this.V = inputVCF
    this.out = swapExt(inputVCF, ".vcf",".a1000g.vcf")
    this.resource :+= qscript.oneKGReleaseMinusBadPools
    this.intervals :+= inputVCF
    this.E = Seq("resource.AC" )

  }

  class AxiomAnnot(inputVCF: File) extends VariantAnnotator with CommandLineGATKArgs {
    this.V = inputVCF
    this.out = swapExt(inputVCF, ".vcf",".axiomAnnot.vcf")
    this.resource :+= qscript.axiomChip
    this.intervals :+= inputVCF
    this.E = Seq("resource.AxiomAC" )
  }

  class OmniAnnot(inputVCF: File) extends VariantAnnotator with CommandLineGATKArgs {
    this.V = inputVCF
    this.out = swapExt(inputVCF, ".vcf",".omniAnnot.vcf")
    this.resource :+= qscript.omniChip
    this.E = Seq("resource.OmniAC" )
  }
  class ExomeChipAnnot(inputVCF: File) extends VariantAnnotator with CommandLineGATKArgs {
    this.V = inputVCF
    this.out = swapExt(inputVCF, ".vcf",".exomeChipAnnot.vcf")
    this.resource :+= qscript.exomeChip
    this.intervals :+= inputVCF
    this.E = Seq("resource.ExomeChipAC" )

  }

  class LOFAnnot(inputVCF: File) extends VariantAnnotator with CommandLineGATKArgs {
    this.V = inputVCF
    this.out = swapExt(inputVCF, ".vcf",".LOF.vcf")
    this.resource :+= qscript.lofFile
    this.intervals :+= inputVCF
    this.E = Seq("resource.LOF" )

  }

  class VToT(inputVCF: File) extends VariantsToTable with CommandLineGATKArgs {
    this.V :+= inputVCF
    this.out = swapExt(inputVCF, ".vcf",".table")
    this.F = Seq("CHROM","POS","REF","ALT","FILTER","AC","QUAL","resource.set","DP","TYPE","NCALLED")
    this.GF :+= "MLPSAC"
    this.allowMissingData = true
    this.showFiltered = true

  }

  def script = {

    var jobList:Seq[PPC] = Seq()


    jobList:+= (new SNPPC("omniMono",baseDir+"outputVCFs/ALL.wgs.1000_control_sites_OmniMono.SNP.sites.interval_list"))
    jobList:+= (new SNPPC("omniPoly",baseDir+"outputVCFs/ALL.wgs.1000_control_sites_OmniPoly.SNP.sites.interval_list"))
    jobList:+= (new SNPPC("exomeChip",baseDir+"outputVCFs/ALL.wgs.1000_control_sites_ExomeChip.SNP.sites.interval_list"))
    jobList:+= (new IndelPC("millsPoly",baseDir+"outputVCFs/ALL.wgs.1000_control_sites_MillsGenotypeInPhase1.INDEL.sites.interval_list"))
    jobList:+= (new SNPPC("lostToImputation",baseDir+"outputVCFs/ALL.wgs.2000_lost_to_Imputation.sites.interval_list"))
    jobList:+= (new IndelPC("afIndels",baseDir+"outputVCFs/ALL.wgs.5000_validation_sites_AF_distributed.indels.sites.interval_list"))
    jobList:+= (new IndelPC("unifIndels",baseDir+"outputVCFs/ALL.wgs.5000_validation_sites_Uniformly_distributed.indels.sites.interval_list"))
    jobList:+= (new SNPPC("afSNPs",baseDir+"outputVCFs/ALL.wgs.8000_validation_sites_AF_distributed.snp.sites.interval_list"))
    jobList:+= (new SNPPC("unifSNPs",baseDir+"outputVCFs/ALL.wgs.8000_validation_sites_Uniformly_distributed.snp.sites.interval_list"))
    jobList:+= (new SNPPC("multiAllelicSNPs",baseDir+"outputVCFs/triallelics.EricBanks.interval_list"))
    jobList:+= (new SNPPC("LOFSNP",baseDir+"outputVCFs/LOF.DanielMacArthur_20130212.fixed.snps.interval_list"))
    jobList:+= (new IndelPC("LOFINDEL",baseDir+"outputVCFs/LOF.DanielMacArthur_20130212.fixed.indels.interval_list"))
    jobList:+= (new IndelPC("multiAllelicIndels",baseDir+"outputVCFs/ALL.wgs.2000_validation_sites_multiAllelicIndels.sites.interval_list"))

    for (job <- jobList) {
      add(job)
      val filt = new Filt(job.out)
      add(filt)
      val annot = new Annot(filt.out)
      add(annot)

      if (job.callName.contains("LOF")) {
        val lofAnnot = new LOFAnnot(annot.out)
        lofAnnot.intervals :+= filt.out
        add(lofAnnot)
        val vtot = new VToT(lofAnnot.out)
        vtot.F :+= "resource.LOF"
        add(vtot)
      }
      else {
        val vtot = new VToT(annot.out)
        add(vtot)
      }
    }
  }
}

