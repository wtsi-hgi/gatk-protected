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

package org.broadinstitute.sting.queue.qscripts.analysis

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import java.io.FileWriter
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.utils.variantcontext.VariantContext
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel

class MultiSampleCallingManuscript extends QScript {
  qscript =>

  @Argument(shortName="outputDir", doc="output directory", required=false)
  var outputDir: String = "results/"

  @Argument(shortName="extraSamples", doc="BAM list of extra samples", required=true)
  var extraSamples: File = _

  @Argument(shortName="sample", doc="Samples to include in Variant Eval", required=false)
  var samples: List[String] = List("NA12878")

  @Argument(shortName = "L", fullName = "intervals", doc="intervals", required=false)
  val myIntervals: List[String] = null;

  @Argument(shortName = "sc", fullName = "scatterCount", doc = "Scatter/Gather jobs to use", required=false)
  val myScatterCount: Int = 1;

  @Argument(shortName = "full", fullName = "full", doc = "Scatter/Gather jobs to use", required=false)
  val fullRun: Boolean = false;

  @Argument(shortName = "indel", fullName = "indel", doc = "Scatter/Gather jobs to use", required=false)
  val includeIndels: Boolean = false;

//  val b37_decoy = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  val bundle = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/")
  val b37 = new File(bundle.getPath + "/human_g1k_v37.fasta")
  val dbSNP_b37 = new File(bundle.getPath + "/dbsnp_132.b37.vcf")
  val dbSNP_b37_129 = new File(bundle.getPath + "/dbsnp_132.b37.excluding_sites_after_129.vcf")
  val hapmap_b37 = new File(bundle.getPath + "/hapmap_3.3.b37.vcf")
  val omni_b37 = new File(bundle.getPath + "/1000G_omni2.5.b37.vcf")

  val resources = new File("resources")
  val CGI_NA12878 = new File(resources.getPath + "/NA12878.CG.b37.snps.vcf")
  val DENOVO_FPs = new File(resources.getPath + "/ceu_merged_validation_data_240610.annotated.b37.FPs.vcf")
  val G1K_PILOT_VALIDATION = new File(resources.getPath + "/1000G.snp.validation.b37.vcf")
  val G1K_PHASEI_VALIDATION = new File(resources.getPath + "/PacBioVsSequenomStatus.vcf")
  val OMNI_MONO = new File(resources.getPath + "/Omni25_monomorphic_1525_samples.b37.vcf")

  var na12878HapMap: File = _
  var na12878OmniTraining: File = _
  var na12878OmniEval: File = _

//  val training_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"
//  val badSites_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.terrible.vcf"
//  val projectConsensus_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/ALL.wgs.projectConsensus_v2b.20101123.snps.sites.vcf"
  val queueLogDir = outputDir.getPath + "/log/"

  val NA12878_BAM_DIRECTORY = "resources/"

  class Target(
          val baseName: String,
          val trancheTarget: Double,
          val nSamples: Int,
          val na12878Depth: Int,
          val nAdditionalSamples: Int,
          val additionalSamples: File) {
    val name = "%s%s_d%d_w%d".format(qscript.outputDir, baseName, na12878Depth, nAdditionalSamples)
    val rawVCF = new File(name + ".raw.vcf")
    val VQSRName = "%s.tranche_%.1f".format(name, trancheTarget)
    val clusterFile = new File(VQSRName + ".clusters")
    val recalibratedVCF = new File(VQSRName + ".recalibrated.vcf")
    val tranchesFile = new File(VQSRName + ".tranches")
    val vqsrRscript = VQSRName + ".vqsr.r"
    val recalFile = new File(VQSRName + ".tranches.recal")
    val evalFile = new File(VQSRName + ".snp.eval")
    val evalIndelFile = new File(VQSRName + ".indel.eval")
    def useBAQ = true

    def na12878bam = new File(NA12878_BAM_DIRECTORY.getPath + "/na12878.%dx.bam".format(na12878Depth))
  }

  def createTargets: List[Target] = {
    var targets: List[Target] = List()
    def addTargets(nAdditionalSamplesList: List[Int], na12878Depths: List[Int], sensitivityThresholds: List[Double]) = {
      for ( nAdditionalSamples <- nAdditionalSamplesList ) {
        for ( na12878Depth <- na12878Depths ) {
          for ( sensitivityThreshold <- sensitivityThresholds ) {
            val file = if ( nAdditionalSamples > 0 ) createAdditionalSampleFile(nAdditionalSamples) else null
            val t = new Target("NA12878", sensitivityThreshold, 1, na12878Depth, nAdditionalSamples, file)
            targets :+= t
          }
        }
      }
    }

    if ( fullRun ) {
      // TODO -- come up with complete list of parameters
      addTargets(List(0, 5, 10, 81, 324), List(4, 64), List(99.0, 99.5, 99.7, 99.9))
    } else {
      addTargets(List(0, 10), List(4,64), List(99.0))
    }

    targets
  }

  def createAdditionalSampleFile(nAdditionalSamples: Int): File = {
    val lines = scala.io.Source.fromFile(extraSamples).getLines().slice(0, nAdditionalSamples).toList

    if ( nAdditionalSamples > lines.size )
      throw new UserException("not enough extra samples provided")

    val destFile = new File(outputDir + "/additional_samples_%d.bam.list".format(nAdditionalSamples))
    val fw = new FileWriter(destFile)
    for ( line <- lines ) fw.write("%s%n".format(line))
    fw.close()

    destFile
  }


  def script = {
    if ( ! outputDir.exists() ) outputDir.mkdirs()

    // set up training and eval data sets
    val na12878HapMapCmd = new SelectNA12878Variants(hapmap_b37)
    val na12878OmniCmd = new SelectNA12878Variants(omni_b37)
    val omniTrainingCmd = new SplitVCF(na12878OmniCmd.out, 0.5)

    add(na12878HapMapCmd, na12878OmniCmd, omniTrainingCmd)

    na12878HapMap = na12878HapMapCmd.out
    na12878OmniTraining = omniTrainingCmd.out1
    na12878OmniEval = omniTrainingCmd.out2

    val targets = createTargets
    for (target <- targets) {
      add(new callVariants(target))
      //add(new indelFilter(target), new indelEvaluation(target))
      add(new VQSR(target))
      add(new applyVQSR(target))
      val na12878NR = new SelectNA12878Variants(target.recalibratedVCF)
      add(na12878NR, new evalVariants(target, na12878NR.out))
    }
  }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    logging_level = "INFO";
    memoryLimit = 2;
    reference_sequence = b37
    intervalsString = myIntervals
  }

  // 1.) Unified Genotyper Base
  class callVariants (t: Target) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.scatterCount = myScatterCount
    //this.nt = 2
    this.dcov = 250
    this.stand_call_conf = 30.0
    this.stand_emit_conf = 30.0
    this.input_file :+= t.na12878bam
    this.input_file :+= t.additionalSamples
    this.D = dbSNP_b37
    this.out = t.rawVCF
    this.glm = if ( includeIndels ) GenotypeLikelihoodsCalculationModel.Model.BOTH else GenotypeLikelihoodsCalculationModel.Model.SNP
    this.baq = if (t.useBAQ) {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY} else {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF}
    this.analysisName = t.name + "_UGs"
    this.jobName =  queueLogDir + t.name + ".snpcall"
  }

//  // todo -- I'd like this to be one process -> snp + indel -> snp + indel (filtered) -> snp (filtered) + indel (filtered)
//  // 2.) Hard Filtering for indels
//  class indelFilter (t: Target) extends VariantFiltration with UNIVERSAL_GATK_ARGS {
//    this.memoryLimit = 2
//    this.scatterCount = 10
//    this.V = t.rawIndelVCF
//    this.out = t.filteredIndelVCF
//    this.filterName ++= List("IndelQD", "IndelReadPosRankSum", "IndelFS")
//    this.filterExpression ++= List("\"QD < 2.0\"", "\"ReadPosRankSum < -20.0\"", "\"FS > 200.0\"")
//    if (t.nSamples >= 10) {
//        this.filterName ++= List("IndelInbreedingCoeff")
//        this.filterExpression ++= List("\"InbreedingCoeff < -0.8\"")
//    }
//    this.analysisName = t.name + "_VF"
//    this.jobName =  queueLogDir + t.name + ".indelfilter"
//  }

  // 3.) Variant Quality Score Recalibration - Generate Recalibration table
  class VQSR(t: Target) extends VariantRecalibrator with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 6
    this.nt = 2
    this.input :+= t.rawVCF
    this.resource :+= new TaggedFile( hapmap_b37, "training=true,truth=true,prior=15.0" )
    this.resource :+= new TaggedFile( na12878OmniTraining, "training=true,truth=true,prior=12.0" )
    //this.resource :+= new TaggedFile( training_1000G, "training=true,prior=10.0" )
    this.resource :+= new TaggedFile( dbSNP_b37_129, "known=true,prior=2.0" )
    //this.resource :+= new TaggedFile( projectConsensus_1000G, "prior=8.0" )
    this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS")
    if(t.nSamples >= 10) { // InbreedingCoeff is a population-wide statistic that requires at least 10 samples to calculate
        this.use_annotation ++= List("InbreedingCoeff")
    }
    this.use_annotation ++= List("DP")
    //this.resource :+= new TaggedFile( badSites_1000G, "bad=true,prior=2.0" )
    this.tranches_file = t.tranchesFile
    this.recal_file = t.recalFile
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.7", "99.5", "99.3", "99.0", "98.0", "97.0", "95.0", "90.0")
    this.rscript_file = t.vqsrRscript
    this.analysisName = t.name + "_VQSR"
    this.jobName = queueLogDir + t.VQSRName + ".VQSR"
  }

  // 4.) Apply the recalibration table to the appropriate tranches
  class applyVQSR (t: Target) extends ApplyRecalibration with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 6
    this.input :+= t.rawVCF
    this.tranches_file = t.tranchesFile
    this.recal_file = t.recalFile
    this.ts_filter_level = t.trancheTarget
    this.out = t.recalibratedVCF
    this.analysisName = t.name + "_AVQSR"
    this.jobName = queueLogDir + t.VQSRName + ".applyVQSR"
  }

  // 5.) Variant Evaluation Base(OPTIONAL)
  class evalVariants(t: Target, vcf: File) extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.eval :+= vcf
    this.eval :+= new TaggedFile( CGI_NA12878, "CGI" )
    this.D = new File(dbSNP_b37_129)

    this.comp :+= new TaggedFile(na12878HapMap, "hapmap3" )
    this.comp :+= new TaggedFile( dbSNP_b37, "dbsnp132" )
    this.comp :+= new TaggedFile( na12878OmniTraining, "omniTraining" )
    this.comp :+= new TaggedFile( na12878OmniEval, "omniEval" )
    this.comp :+= new TaggedFile( CGI_NA12878, "CGI" )
    this.comp :+= new TaggedFile( DENOVO_FPs, "DeNovoFP" )
    this.comp :+= new TaggedFile( G1K_PILOT_VALIDATION, "1000G_pilot_val" )
    this.comp :+= new TaggedFile( G1K_PHASEI_VALIDATION, "1000G_phase1_val" )
    this.comp :+= new TaggedFile( OMNI_MONO, "OmniMono" )

    //this.sample = List("NA12878")
    this.out =  swapExt(outputDir, vcf, ".vcf", ".eval")
    this.analysisName = t.name + "_VEs"
    //this.evalModule :+= "IndelStatistics"
    this.ST = List("Filter")
    this.jobName = queueLogDir + t.VQSRName + ".eval"
  }

  class SelectNA12878Variants(vcf: File) extends SelectVariants with UNIVERSAL_GATK_ARGS {
    this.variant = vcf
    this.out = swapExt(outputDir, vcf, ".vcf", ".na12878.vcf")
    this.sample_name = List("NA12878")
    this.selectType = List(VariantContext.Type.SNP)
    this.excludeNonVariants = true // we only want sites polymorphic in the sample
  }

  class SplitVCF(vcf: File, fractionForTraining: Double) extends RandomlySplitVariants with UNIVERSAL_GATK_ARGS {
    this.variant = vcf
    this.out1 = swapExt(outputDir, vcf, ".vcf", ".training.vcf")
    this.out2 = swapExt(outputDir, vcf, ".vcf", ".eval.vcf")
    this.fraction = fractionForTraining
  }
}
