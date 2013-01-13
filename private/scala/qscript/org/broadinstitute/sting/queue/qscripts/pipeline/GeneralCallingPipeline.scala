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

package org.broadinstitute.sting.queue.qscripts.pipeline

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport
import org.broadinstitute.sting.queue.extensions.gatk.BaseRecalibrator
import org.broadinstitute.sting.queue.extensions.gatk.PrintReads
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.function.{RetryMemoryLimit, ListWriterFunction}
import org.broadinstitute.sting.queue.function._


class GeneralCallingPipeline extends QScript {
  qscript =>

  // Required arguments.  All initialized to empty values.

  @Input(shortName="R", doc="The reference file for the bam files.")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(shortName="I", doc="Bam file to genotype.")
  var bamFile: File = _

  @Argument(shortName="V", doc="The pipeline version", required=true)
  var pipelineVersion: Int = 1

  @Argument(shortName="mode", doc="the class of variation to model: SNP, INDEL, BOTH", required=true)
  var mode: String = _

  // The following arguments are all optional.

  @Argument(shortName = "RPath", doc="RPath", required=false)
  var RPath: File = new File("../R")

  @Argument(shortName="scatterCount", doc="set the number of the scattered jobs", required=false)
  var scatterCount: Int = 0

  @Argument(shortName="outputDir", doc="output directory", required=false)
  var outputDir: String = "./tmp"

  @Input(shortName="L", doc="An optional file with a list of intervals to proccess.",  required=false)
  var intervals: File = _

  @Argument(shortName="filter", doc="A optional list of filter names.", required=false)
  var filterNames: List[String] = Nil // Nil is an empty List, versus null which means a non-existent List.

  @Argument(shortName="noBAQ", doc="turns off BAQ calculation", required=false)
  var noBAQ: Boolean = false

  @Argument(shortName="mbq", doc="The minimum Phred-Scaled quality score threshold to be considered a good base.", required=false)
  var minimumBaseQuality: Int = -1

  @Argument(shortName="filterExpression", doc="An optional list of filter expressions.", required=false)
  var filterExpressions: List[String] = Nil

  @Argument(shortName="sample", doc="Samples to include in Variant Eval", required=false)
  var samples: List[String] = Nil

  @Argument(shortName="deletions", doc="Maximum deletion fraction allowed at a site to call a genotype.", required=false)
  var deletions: Double = -1

  @Input(shortName = "XL", doc="Exclude intervals list", required=false)
  var excludeIntervals: List[File] = Nil

  @Input(shortName="ped", fullName="pedigree", doc="PED file of the family", required=false)
  var ped: File = _

  @Argument(doc="Subdirectory to store the reduced bams. By default set to 'reduced'.", shortName="bamDir", required=false)
  var bamDir = "reducedBAMs/"

  @Argument(doc="Reduce reads memory limit.", shortName="rrMem", required=false)
  var reduceReadsMemoryLimit = 4

  /************* invlude/exclude steps of the pipeline ***********************/

  @Argument(shortName = "doNotCall", doc="don't call any events", required=false )
  var doNotCall: Boolean = false

  @Argument(shortName = "doNotUseVQSR", doc="don't call preform VQSR on the called vcf", required=false )
  var doNotUseVQSR: Boolean = false


  @Argument(shortName="useBQSR.2.0", doc="turn on a first step of using 2.0 BQSR on the input bam file", required=false)
  var useBQSR2: Boolean = false

  @Argument(shortName="usePhaseBT", doc="apply phaseByTrasmission on the output vcf", required=false)
  var usePhaseBT: Boolean = false

  @Argument(shortName="createEvalSummaryReport", doc="create a pdf file with eval summery tables of the output vcf", required=false)
  var createEvalSummaryReport: Boolean = false

  @Argument(shortName="createFullPostQCReport", doc="create a pdf file with full QC eval report of the output vcf", required=false)
  var createFullPostQCReport: Boolean = false

  @Argument(shortName="CallReduceReads", doc="produces reduce reads bam files", required=false)
  var CallReduceReads: Boolean = false



  val dbSNP_135 = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_135.b37.vcf"  // Best Practices v4
  val hapmap = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/hapmap_3.3.b37.sites.vcf"                       // Best Practices v4
  val omni_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_2141_samples.b37.vcf"    // Best Practices v4
  val training_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"  // from the MethodDevelopmentCallingPipeline scala script
  val projectConsensus_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/ALL.wgs.projectConsensus_v2b.20101123.snps.sites.vcf"  // from the MethodDevelopmentCallingPipeline scala script
  val indelGoldStandardCallset  = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/Mills_and_1000G_gold_standard.indels.b37.vcf" // Best Practices v4
  val dbSNP_129 = "/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.vcf"
  val omni_mono: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_monomorphic_2141_samples.b37.vcf"
  val MVL_FP: String =  "/humgen/gsa-hpprojects/GATK/badSitesDB/Autism.HighMVLR.Exome.indels.vcf"
  val CEUTrio_ped = "/broad/hptmp/ami/tmp/CEUTrio.ped" //todo: move to a proper location


  val queueLogDir = ".qlog/"
  val noET_key = "/humgen/gsa-hpprojects/GATK/data/gatk_user_keys/gsamembers_broadinstitute.org.key"                // TODO: remove before we make it public!!!!


trait BaseCommandArguments extends CommandLineGATK {
    this.logging_level = "INFO"
    phone_home = GATKRunReport.PhoneHomeOption.NO_ET   // TODO: remove before we make it public!!!!
    gatk_key = noET_key                                // TODO: remove before we make it public!!!!
    this.jobQueue = "gsa"

    this.reference_sequence = qscript.referenceFile
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    this.excludeIntervals = qscript.excludeIntervals
    this.memoryLimit = 2

}



trait HaplotypeCallerArguments extends BaseCommandArguments {
	
}

trait BaseBQSR extends BaseCommandArguments{

}

// 0) BQSR v 2.0
class BQSR( bam: File, recal: File ) extends BaseRecalibrator with BaseBQSR{

  this.input_file :+= bam
	this.out = recal
	this.knownSites ++= List(new File(indelGoldStandardCallset))
 	this.knownSites ++= List(new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.indels.b37.vcf"))
	this.knownSites ++= List(new File("/humgen/gsa-hpprojects/dev/carneiro/bqsr/data/projectConsensus.snps.vcf"))
	this.knownSites ++= List(new File(dbSNP_135))
	this.memoryLimit = 8
	this.qq = 0
	this.mcs = 2
	this.ics = 3
	//this.np = true
	this.javaGCThreads = 4
	this.scatterCount = 200 //if(qscript.scatterCount == 0) 200 else qscript.scatterCount
}

class RecalBAM( bam: File, recal: File ) extends PrintReads with BaseBQSR{

  this.input_file :+= bam
	this.out = qscript.outputDir + recal +".out"
	this.BQSR = recal
	this.memoryLimit = 6
	this.qq = 0
	this.javaGCThreads = 4
	this.scatterCount = 200 //if(qscript.scatterCount == 0) 200 else qscript.scatterCount
}

// 1.) Unified Genotyper Base
class UGBase extends UnifiedGenotyper with BaseCommandArguments {

    this.scatterCount = if(qscript.scatterCount == 0) 80 else qscript.scatterCount
    this.nt = 2
    this.dcov = 250 
    this.stand_call_conf = 30.0
    this.stand_emit_conf = 30.0
    this.D = new File(dbSNP_135)
}

  // 1a.) Call SNPs with UG
  class snpCaller(name: String, inputBamFile: File) extends UGBase with BaseCommandArguments {
    this.input_file :+= inputBamFile
    if (qscript.minimumBaseQuality >= 0)
      this.min_base_quality_score = qscript.minimumBaseQuality
    if (qscript.deletions >= 0)
      this.max_deletion_fraction = qscript.deletions
    this.out = qscript.outputDir + name + "snp.unfiltered.vcf"
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
    this.baq = if (noBAQ) {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF} else {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY}
    this.analysisName = "CEUTrio_UGs"
    this.jobName =  queueLogDir + "CEU_Trio.snpcall"
  }

  // 1b.) Call Indels with UG
  class indelCaller(name: String, inputBamFile: File) extends UGBase with BaseCommandArguments {
    this.input_file :+= inputBamFile
    this.memoryLimit = 6
    this.out = qscript.outputDir + name + "indel.unfiltered.vcf"
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
    this.analysisName = "CEUTrio_UGi"
    this.jobName =  queueLogDir + "CEUTrio.indelcall"
  }

// 1.) HaplotypeCaller 
class HCBase(name:String, inputHCFile: File) extends HaplotypeCaller with HaplotypeCallerArguments {
    this.excludeIntervals = excludeIntervals
    this.scatterCount = if(qscript.scatterCount == 0) 200 else qscript.scatterCount
    this.input_file :+= inputHCFile
    this.out = qscript.outputDir + name + "HaplotypeCaller.vcf"
    this.analysisName = "HaplotypeCaller"
    this.stand_call_conf = 10.0
    this.stand_emit_conf = 10.0
    this.minPruning = 2
    this.javaGCThreads = 4
    this.memoryLimit = 8 //4	
}
 

  /****************************************************************************************
  *                3.)   VariantRecalibrator                                              *
  *****************************************************************************************/

  class VQSRBase extends VariantRecalibrator with BaseCommandArguments {
   
    this.nt = 2
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
}

// 3a)
class snpRecal(name: String) extends VQSRBase with BaseCommandArguments{
  this.input :+= qscript.outputDir + name + "snp.unfiltered.vcf"
	this.resource :+= new TaggedFile( hapmap, "known=false,training=true,truth=true,prior=15.0" ) 
  this.resource :+= new TaggedFile( omni_b37, "known=false,training=true,truth=true,prior=12.0" ) // truth=false on the bast practices v4
	this.resource :+= new TaggedFile( training_1000G, "known=false,training=true,prior=10.0" )	// not part of the bast practices v4
	this.resource :+= new TaggedFile( dbSNP_135, "known=true,training=false,truth=false,prior=2.0" )    // prior=6.0 on the bast practices v4
	this.resource :+= new TaggedFile( projectConsensus_1000G, "prior=8.0" )				// not part of the bast practices v4
  this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "DP")   //remove DP when generalizing for exomes
 	// this.use_annotation ++= List("InbreedingCoeff")   //need more then 10 samples
	this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
 	this.memoryLimit = 8
	this.tranches_file = qscript.outputDir + name + "snp.tranches"
	this.recal_file = qscript.outputDir +  name + "snp.recal"
	this.rscript_file = qscript.outputDir + name + ".snp.vqsr.R"
	this.analysisName = "CEUTrio_VQSRs"
	this.jobName = queueLogDir + "CEUTrio.snprecal"

}

// 3b)
class indelRecal(name:String) extends VQSRBase with BaseCommandArguments {

	this.input :+= qscript.outputDir + name + "indel.unfiltered.vcf"
	this.resource :+= new TaggedFile(indelGoldStandardCallset, "known=false,training=true,truth=true,prior=12.0" ) // known=true on the bast practices v4
	this.resource :+= new TaggedFile( dbSNP_135, "known=true,prior=2.0" )  						// not part of the bast practices v4
	this.use_annotation ++= List("QD", "HaplotypeScore", "ReadPosRankSum", "FS")
	// this.use_annotation ++= List("InbreedingCoeff")   //need more then 10 samples   	

 	this.tranches_file = qscript.outputDir + name + "indel.tranches"
 	this.recal_file = qscript.outputDir + name + "indel.recal"
 	this.rscript_file = qscript.outputDir + name + ".indel.vqsr.R"
 	this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
  this.analysisName = "CEUTrio_VQSRi"
 	this.jobName = queueLogDir + "CEUTrio.indelrecal"
}

//3c)
class RecalBoth(inputVCF: File) extends VariantRecalibrator with HaplotypeCallerArguments {

  this.input :+= inputVCF
  this.resource :+= new TaggedFile( hapmap, "known=false,training=true,truth=true,prior=15.0" )
  this.resource :+= new TaggedFile( omni_b37, "known=false,training=true,truth=false,prior=12.0" ) // used to be truth=true
  this.resource :+= new TaggedFile( dbSNP_135, "known=true,training=false,truth=false,prior=6.0" )
  this.resource :+= new TaggedFile( indelGoldStandardCallset, "known=true,training=true,truth=true,prior=12.0" )  //used to be known=false, prior=14.0
  this.use_annotation ++= List("QD", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "DP", "ClippingRankSum")  //remove DP for exomes
  // this.use_annotation ++= List("InbreedingCoeff")   //need more then 10 samples
  this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.BOTH
  this.allPoly = true
  this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.0", "97.9", "97.5", "97.0", "96.0","95.0", "94.0", "93.0", "92.0","90.0")
  this.recalFile = qscript.outputDir + swapExt(inputVCF, "vcf", "both.recalFile")
  this.tranchesFile = qscript.outputDir + swapExt(inputVCF, "vcf", "both.tranches")
  this.rscriptFile = qscript.outputDir + swapExt(inputVCF, "vcf", "both.R")
  this.percentBad = 0.012
  this.mG = 7
  this.ts_filter_level = 97.0
  this.memoryLimit = 8

}

    /****************************************************************************************
    *               4.) Apply the recalibration table to the appropriate tranches           *
    *****************************************************************************************/

class applyVQSRBase extends ApplyRecalibration with BaseCommandArguments  {
  this.memoryLimit = 6
}

class applySnpVQSR(name:String) extends applyVQSRBase with BaseCommandArguments {

  this.input :+= qscript.outputDir + name + "snp.unfiltered.vcf"
  this.tranches_file = qscript.outputDir + name + "snp.tranches"
  this.recal_file = qscript.outputDir +  name + "snp.recal"
  this.ts_filter_level = 99.0
  this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
  this.out = qscript.outputDir +  name + "snp.recalibrated.filtered.vcf"
  this.analysisName = "CEUTrio_AVQSRs"
  this.jobName = queueLogDir + "CEUTrio.snpcut"
}

class applyIndelVQSR(name:String) extends applyVQSRBase with BaseCommandArguments {

    this.input :+= qscript.outputDir + name + "indel.unfiltered.vcf"
    this.tranches_file = qscript.outputDir + name +"indel.tranches"
    this.recal_file = qscript.outputDir + name +"indel.recal"
    this.ts_filter_level = 95.0  							//best practices v4
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
    this.out = qscript.outputDir + name + "indel.recalibrated.filtered.vcf"
    this.analysisName = "CEUTrio_AVQSRi"
    this.jobName = queueLogDir + "CEUTrio.indelcut"
}

class CutBoth(inputVCF: File) extends applyVQSRBase with HaplotypeCallerArguments {

    this.input :+= inputVCF
    this.memoryLimit = 8
    this.num_threads = 1
    this.recalFile = qscript.outputDir + swapExt(inputVCF, "vcf", "both.recalFile")
    this.tranchesFile = qscript.outputDir + swapExt(inputVCF, "vcf", "both.tranches")
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.BOTH
    this.ts_filter_level = 97.0
    this.out = qscript.outputDir + swapExt(inputVCF, "vcf", "both.recalibrated.filtered.vcf")	
    this.analysisName = "CEUTrio_AVQSRb"
    this.jobName = queueLogDir + "CEUTrio.bothcut"
}

    /****************************************************************************************
    *               5) Combine Snps and Indels for UG if mode == BOTH                       *
    *****************************************************************************************/

class CombineSNPsIndels(name:String) extends CombineVariants with BaseCommandArguments {
    this.variant :+= TaggedFile(new File(qscript.outputDir + name + "indel.recalibrated.filtered.vcf"), "indels")
    this.variant :+= TaggedFile(new File(qscript.outputDir + name + "snp.recalibrated.filtered.vcf"), "snps")
    this.filteredrecordsmergetype = org.broadinstitute.variant.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    this.assumeIdenticalSamples = true
    this.out = qscript.outputDir + name +  "both.recalibrated.filtered.vcf"
}



/****************************************************************************************
*                6.) Variant Evaluation Base(OPTIONAL)                                  *
*****************************************************************************************/

class EvalBase extends VariantEval with BaseCommandArguments {
      this.memoryLimit = 3
      this.comp :+= new TaggedFile(hapmap, "hapmap" )
      this.D = new File(dbSNP_129)
      this.sample = samples
}


  // 6a.) SNP Evaluation (OPTIONAL) based on the recalbrated vcf
class snpEvaluation(name:String) extends EvalBase with BaseCommandArguments {
     this.comp :+= new TaggedFile( omni_b37, "omni" )
     this.eval :+= qscript.outputDir +  name + "snp.recalibrated.filtered.vcf"
     this.out = qscript.outputDir + name+ "eval"
     this.analysisName = "CEUTrio_VEs"
     this.jobName = queueLogDir + "CEUTrio.snpeval"
}

  // 6b.) Indel Evaluation (OPTIONAL)
class indelEvaluation(name:String) extends EvalBase with BaseCommandArguments{
    this.eval :+= qscript.outputDir + name +"indel.recalibrated.filtered.vcf"
    this.comp :+= new TaggedFile(indelGoldStandardCallset, "indelGS" )
    this.noEV = true
    this.evalModule = List("CompOverlap", "CountVariants", "TiTvVariantEvaluator", "ValidationReport", "IndelSummary")
    this.out = qscript.outputDir + name + "eval"
    this.analysisName = "CEUTrio_VEi"
    this.jobName = queueLogDir + "CEUTrio.indeleval"
}

// todo  -- should accept separate indel and snp vcf's, right now script will assume they're combined in one
class Eval(evalVCF: File, prefix: String, extraStrats: Seq[String]) extends VariantEval with BaseCommandArguments{

    this.eval :+= evalVCF
    this.dbsnp = new File(dbSNP_129)
    this.doNotUseAllStandardModules = true
    this.gold = qscript.indelGoldStandardCallset
    this.comp :+= new TaggedFile( omni_b37, "omni" )
    this.comp :+= new TaggedFile( omni_mono, "omni_mono" )
    this.comp :+= new TaggedFile( qscript.indelGoldStandardCallset, "GSindels" )
    this.comp :+= new TaggedFile( MVL_FP, "fp_MVL" )
    this.evalModule = List("TiTvVariantEvaluator", "CountVariants", "CompOverlap", "IndelSummary", "MultiallelicSummary","ValidationReport")
    this.doNotUseAllStandardStratifications = true
    this.stratificationModule = Seq("EvalRod", "CompRod", "Novelty", "FunctionalClass") ++ extraStrats
    this.num_threads = 1
    this.memoryLimit = 8
    this.out = outputDir + swapExt(evalVCF, ".vcf", prefix + ".eval")
    this.sample = samples
    this.analysisName = "CEUTrio_VE"
    this.jobName = queueLogDir + "CEUTrio.eval"
}

class QCSummaryRScript(@Input var vcf: File, @Input var bySite: File) extends CommandLineFunction {
    @Output var pdf: File = outputDir + swapExt(vcf, ".vcf", ".pdf")
    private val project = vcf.getName.stripSuffix(".vcf")
    def commandLine = "Rscript %s/variantCallQC_summaryTablesOnly.R %s %s %s".format(RPath, project, bySite, pdf)
}

class QCRScript(@Input var vcf: File, @Input var byAC: File, @Input var bySite: File, @Input var indelQC: File) extends CommandLineFunction {
    @Output var pdf: File = swapExt(vcf, ".vcf", ".pdf")
    private val project = vcf.getName.stripSuffix(".vcf")
    def commandLine = "Rscript %s/variantCallQC.R %s %s %s %s %s".format(RPath, project, bySite, byAC, indelQC, pdf)
}

//------------------------------------------------------------------------------------ //
//                      7) PhaseByTrasmission                                          //
//------------------------------------------------------------------------------------ //
class PBT(inputVCF: File) extends PhaseByTransmission with BaseCommandArguments {
  require(qscript.ped != null && qscript.ped.getName.endsWith(".ped"), "-ped/--pedigree must be specified as <pedigreeFileName>.ped")
  this.pedigree = Seq(qscript.ped)
  this.variant = inputVCF
  this.MendelianViolationsFile = qscript.outputDir + swapExt(variant,"vcf","mendelianViolationsFile.txt")
  this.out = qscript.outputDir + swapExt(variant,"vcf","phaseByTransmission.vcf")
  this.jobName = queueLogDir + "CEUTrio.phaseBT"
}


    /****************************************************************************************
    *                script                                                                 *
    *****************************************************************************************/

 def script() {
    var inputBamFile = qscript.bamFile //todo make sure it work both for bam file and bam list file
    var name = "noName"_
    if(inputBamFile.getName.endsWith(".bam.list")){
      name = inputBamFile.getName.stripSuffix(".bam.list")
    }
    else if (inputBamFile.getName.endsWith(".bam")) {
       name = inputBamFile.getName.stripSuffix(".bam")
    }

  if (useBQSR2){

    val bams = QScriptUtils.createSeqFromFile(inputBamFile)
    var recalBams = Seq.empty[File]
    for (bam <- bams) {
          val recalFile = qscript.outputDir + swapExt(bam, ".bam", ".recal.grp")
          add(new BQSR(bam, recalFile))
          val recalBam = new RecalBAM(bam, recalFile)
          add(recalBam)
          recalBams :+= recalBam.out
    }

    val mergeBamList = new ListWriterFunction
    mergeBamList.inputFiles = recalBams
    mergeBamList.listFile = "CEUTrio.recal.bam.list"
    add(mergeBamList)
    inputBamFile = mergeBamList.listFile
  }

  if (CallReduceReads){
      var bams: Seq[Tuple2[File, File]] = Nil
      var reducedBams = Seq.empty[File]
      if (inputBamFile != null) {
          for (originalBam: File <- io.Source.fromFile(inputBamFile).getLines().toSeq.map(new File(_))) {
            val reducedBam: File = new File(new File(bamDir, "external"), swapExt(originalBam, ".bam", ".reduced.bam").getName)
            bams :+= Tuple2(originalBam, reducedBam)
          }
        }

        for ((originalBam, reducedBam) <- bams) {
          val reduce = new ReduceReads with BadMate with RetryMemoryLimit
          reduce.memoryLimit = qscript.reduceReadsMemoryLimit
          reduce.reference_sequence = qscript.referenceFile
          reduce.input_file = Seq(originalBam)
          reduce.intervals = Seq(qscript.intervals)
          // reduce.interval_padding = 50
          reduce.out = reducedBam
          reducedBams :+= reduce.out
          add(reduce)
        }

      val mergeReducedBamList = new ListWriterFunction
      mergeReducedBamList.inputFiles = reducedBams
      mergeReducedBamList.listFile = "CEUTrio.reduced.bam.list"
      add(mergeReducedBamList)
      inputBamFile = mergeReducedBamList.listFile

    }




  // Create the functions that we may run depending on options.
  if (pipelineVersion == 1){
	if (mode == "SNP"){
		val UGgenotyper = new snpCaller(name,inputBamFile)
		val snpRecalibrator = new snpRecal(name)
		val snpApplyVQSR = new applySnpVQSR(name)
		val evaluateSnp = new snpEvaluation(name)
 
		if (! doNotCall){
      add (UGgenotyper)
    }
    if (! doNotUseVQSR){
      add (snpRecalibrator)
      add (snpApplyVQSR)
		  add (evaluateSnp)
    }
    if (createEvalSummaryReport){
      val evalVCF = snpApplyVQSR.out
      val bySample = new Eval( evalVCF, ".bySample", Seq("Sample"))
      add(bySample)
      val qc = new QCSummaryRScript(evalVCF, bySample.out)
      add(qc)
    }
    if(usePhaseBT){add(new PBT(snpApplyVQSR.out)) }
	}

	if (mode == "INDEL"){
		val UGgenotyper = new indelCaller(name,inputBamFile)
		val indelRecalibrator = new indelRecal(name)
		val indelApplyVQSR = new applyIndelVQSR(name)
		val evaluateIndel = new indelEvaluation(name)

    if (! doNotCall){
      add (UGgenotyper)
    }
    if (! doNotUseVQSR){
      add (indelRecalibrator)
      add (indelApplyVQSR)
		  add (evaluateIndel)
    }

    if (createEvalSummaryReport){
      val evalVCF = indelApplyVQSR.out
      val bySample = new Eval( evalVCF, ".bySample", Seq("Sample"))
      add(bySample)
      val qc = new QCSummaryRScript(evalVCF, bySample.out)
      add(qc)
    }

    if(usePhaseBT){add(new PBT(indelApplyVQSR.out)) }
	}

  // currently running separate runs for SNPs and INDELS
	if (mode == "BOTH"){
		val snpUGgenotyper = new snpCaller(name,inputBamFile)
                val snpRecalibrator = new snpRecal(name)
                val snpApplyVQSR = new applySnpVQSR(name)
		val indelUGgenotyper = new indelCaller(name,inputBamFile)
		val indelRecalibrator = new indelRecal(name)
		val indelApplyVQSR = new applyIndelVQSR(name)
	  val combineSNPsIndels = new CombineSNPsIndels(name)

    if (! doNotCall){
      add (snpUGgenotyper)
      add (indelUGgenotyper)
    }
    if (! doNotUseVQSR){
      add (snpRecalibrator)
		  add (snpApplyVQSR)

		  add (indelRecalibrator)
		  add (indelApplyVQSR)
		  add (combineSNPsIndels)
    }
    if (createEvalSummaryReport || createFullPostQCReport){
      val evalVCF: File = combineSNPsIndels.out

      val bySample = new Eval( evalVCF, ".bySample", Seq("Sample"))
      add(bySample)

      if (createFullPostQCReport){
        val byAC = new Eval(evalVCF, ".byAC", Seq("AlleleCount"))
        add(byAC)
        val indelQC = new Eval(evalVCF, ".indelQC", Seq("Sample"))
        indelQC.stratificationModule = Seq("EvalRod", "CompRod", "Sample", "TandemRepeat", "OneBPIndel")
        indelQC.evalModule = List("IndelSummary", "IndelLengthHistogram")
        add(indelQC)
        val qc = new QCRScript(evalVCF, byAC.out, bySample.out, indelQC.out)
        add(qc)

      }
      if (createEvalSummaryReport){
        val qc = new QCSummaryRScript(evalVCF, bySample.out)
        add(qc)
      }
    }


    if(usePhaseBT){add(new PBT(combineSNPsIndels.out)) }
	}
 }

 else {


		val HCgenotyper = new HCBase(name,inputBamFile)
		val HC_vcfFile = HCgenotyper.out
		val HCRecalibrator = new RecalBoth(HC_vcfFile)
		val HCApplyVQSR = new CutBoth(HC_vcfFile)

    if (! doNotCall){
      add (HCgenotyper)
    }
    if (! doNotUseVQSR){
		  add (HCRecalibrator)
		  add (HCApplyVQSR)
    }
    if (createEvalSummaryReport || createFullPostQCReport){
      val evalVCF: File = HCApplyVQSR.out

      val bySample = new Eval( evalVCF, ".bySample", Seq("Sample"))
      add(bySample)

      if (createFullPostQCReport){
        val byAC = new Eval(evalVCF, ".byAC", Seq("AlleleCount"))
        add(byAC)
        val indelQC = new Eval(evalVCF, ".indelQC", Seq("Sample"))
        indelQC.stratificationModule = Seq("EvalRod", "CompRod", "Sample", "TandemRepeat", "OneBPIndel")
        indelQC.evalModule = List("IndelSummary", "IndelLengthHistogram")
        add(indelQC)
        val qc = new QCRScript(evalVCF, byAC.out, bySample.out, indelQC.out)
        add(qc)

      }
      if (createEvalSummaryReport){
        val qc = new QCSummaryRScript(evalVCF, bySample.out)  //val qc = new QCRScript(evalVCF, byAC.out, bySample.out, indelQC.out)
        add(qc)
      }
    }
    if(usePhaseBT){add(new PBT(HCApplyVQSR.out)) }

 }
}



// ingor this section for now!!
// This trait allows us set the variables below in one place,
// and then reuse this trait on each CommandLineGATK function below.

// currently it is hard coded in the script to run on the CEU trio as the target. In the next version of the script I will add the other built-in targets.
class Target(
              val baseName: String,
              val reference: File,
              val dbsnpFile: String,
              val hapmapFile: String,
              // val maskFile: String,
              val bamList: File,
              val intervals: String,
              // val indelTranchTarget: Double,
              val snpTrancheTarget: Double,
              val isLowpass: Boolean,
              val isExome: Boolean,
              val nSamples: Int) {

}


//val targetDataSets: Map[String, Target] = Map(
//"NA13878_test" -> new Target(	"NA13878.test",
//				qscript.referenceFile,
//				dbSNP,
//				hapmap,
//
//				)
//
//"CEUTrio_wgs_decoy" -> new Target("CEUTrio.HiSeq.WGS.b37_decoy", b37_decoy, dbSNP_b37, hapmap_b37, indelMask_b37,
//              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.list"),
//       //       new File("/humgen/gsa-hpprojects/dev/carneiro/trio/analysis/snps/CEUTrio.WEx.filtered.vcf"),                  // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
//              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg19.intervals", 90.0, 99.0, !lowPass, !exome, 3),
//)


}