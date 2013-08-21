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

package org.broadinstitute.sting.queue.qscripts.incrementalCaller

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function._
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils.FilteredRecordMergeType
import org.broadinstitute.variant.variantcontext.VariantContext

class AssessJointCaller extends QScript {
  qscript =>

  var referenceFile: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37_decoy.fasta")
  var exomeIntervals: File = new File("/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list")

  // Required arguments.  All initialized to empty values.

  @Argument(shortName="outputDir", doc="output directory", required=false)
  var outputDir: String = "./results"

  @Input(shortName="I", doc="Bam file or .list of BAM files to genotype.")
  var jointBams: File = _

  @Input(shortName="FIN", doc="Bam file or .list of BAM files for NA12878.")
  var NA12878FullBam: File = _

  @Input(shortName="IN", doc="Bam file or .list of BAM files for NA12878.")
  var NA12878ReducedBam: File = _

  @Input(shortName="CE", doc="Bam file or .list of BAM files for NA12878.", required=false)
  var ClinicalNA12878Bam: File = new File("/seq/picard_aggregation/DEV-2930/Pond-233795,_Spin_Column_Single_Index_Midi_Plate/current/Pond-233795,_Spin_Column_Single_Index_Midi_Plate.bam")

  @Argument(shortName="sc", doc="set the number of the scattered jobs", required=false)
  var scatterCount: Int = 0

  @Argument(shortName="L", doc="An optional file with a list of intervals to proccess.",  required=false)
  var intervals: String = null

  @Argument(shortName="recalIndels", doc="If true, we will recalibrate indels, otherwise we'll just pass them straight through",  required=false)
  var recalIndels: Boolean = false

  val latestdbSNP = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_137.b37.vcf"  // Best Practices v4
  val hapmapSites = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/hapmap_3.3.b37.vcf"                       // Best Practices v4
  val omni_b37_sites = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_2141_samples.b37.vcf"    // Best Practices v4
  val training_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"  // from the MethodDevelopmentCallingPipeline scala script
  val indelGoldStandardCallset  = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/Mills_and_1000G_gold_standard.indels.b37.vcf" // Best Practices v4
  val dbSNP_129 = "/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.vcf"

  var nSamples: Int = 1

  trait BaseCommandArguments extends CommandLineGATK with RetryMemoryLimit {
    this.logging_level = "INFO"
    this.reference_sequence = qscript.referenceFile
    this.intervals :+= exomeIntervals
    if ( qscript.intervals != null ) {
      this.intervalsString :+= qscript.intervals
      this.interval_set_rule = org.broadinstitute.sting.utils.interval.IntervalSetRule.INTERSECTION
    }
    this.interval_padding = 50
    this.memoryLimit = 2
    this.javaGCThreads = 4
  }

  // 1a.) Call SNPs with UG
  class MyUnifiedGenotyper(name: String) extends UnifiedGenotyper with BaseCommandArguments {
    this.out = qscript.outputDir + "/" + name + ".UnifiedGenotyper.unfiltered.vcf"
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.analysisName = "CEUTrio_UGs"
    if (qscript.scatterCount > 0) this.scatterCount = qscript.scatterCount
    this.D = new File(latestdbSNP)
    this.onlyEmitSamples :+= "NA12878"
  }

  // 1.) HaplotypeCaller
  class MyHaplotypeCaller(name:String) extends HaplotypeCaller with BaseCommandArguments {
    this.out = qscript.outputDir + "/" + name + ".HaplotypeCaller.unfiltered.vcf"
    this.analysisName = "HaplotypeCaller"
    this.dbsnp = latestdbSNP
    this.memoryLimit = 4
    if (qscript.scatterCount > 0) this.scatterCount = qscript.scatterCount
  }

  /****************************************************************************************
    *                3.)   VariantRecalibrator                                              *
    *****************************************************************************************/

  class VQSRBase(vcf:File) extends VariantRecalibrator with BaseCommandArguments {
    this.input :+= vcf
    //this.nt = 4
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.8", "99.7", "99.5", "99.0", "98.5", "98.0", "97.0", "95.0", "90.0")
    this.memoryLimit = 4
    this.tranches_file = swapExt(outputDir, vcf, ".vcf", ".tranches")
    this.recal_file = swapExt(outputDir, vcf, ".vcf", ".recal")
    this.rscript_file = swapExt(outputDir, vcf, ".vcf", ".vqsr.R")

    this.use_annotation ++= List("FS", "ReadPosRankSum")
    //this.use_annotation ++= List("DP") // not appropriate for exomes
    if ( nSamples > 100 )
      this.use_annotation ++= List("InbreedingCoeff")
  }

  // 3a)
  class snpRecal(snpVCF: File) extends VQSRBase(snpVCF) with BaseCommandArguments{
    this.resource :+= new TaggedFile( hapmapSites, "known=false,training=true,truth=true,prior=15.0" )
    this.resource :+= new TaggedFile( omni_b37_sites, "known=false,training=true,truth=true,prior=12.0" ) // truth=false on the bast practices v4
    this.resource :+= new TaggedFile( training_1000G, "known=false,training=true,prior=10.0" )	// not part of the bast practices v4
    this.resource :+= new TaggedFile( dbSNP_129, "known=true,training=false,truth=false,prior=2.0" )    // prior=6.0 on the bast practices v4
    this.use_annotation ++= List("QD")
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
    this.analysisName = "CEUTrio_VQSRs"
  }

  // 3b)
  class indelRecal(indelVCF: String) extends VQSRBase(indelVCF) with BaseCommandArguments {
    this.resource :+= new TaggedFile( indelGoldStandardCallset, "known=false,training=true,truth=true,prior=12.0" ) // known=true on the bast practices v4
    this.resource :+= new TaggedFile( latestdbSNP, "known=true,prior=2.0" )  						// not part of the bast practices v4
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
    this.analysisName = "CEUTrio_VQSRi"
    this.maxGaussians = 4
  }

  /****************************************************************************************
    *               4.) Apply the recalibration table to the appropriate tranches         *
    ***************************************************************************************/

  class applyVQSRBase(vqsr: VariantRecalibrator) extends ApplyRecalibration with BaseCommandArguments  {
    val in = vqsr.input(0)
    this.input :+= in
    this.tranches_file = vqsr.tranches_file
    this.recal_file = vqsr.recal_file
    this.memoryLimit = 4
    this.out = swapExt(outputDir, in, ".unfiltered.vcf", ".recalibrated.vcf")
  }

  class applySnpVQSR(vqsr: VariantRecalibrator) extends applyVQSRBase(vqsr) with BaseCommandArguments {
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
    this.analysisName = "CEUTrio_AVQSRs"
    this.ts_filter_level = 99.5
  }

  class applyIndelVQSR(vqsr: VariantRecalibrator) extends applyVQSRBase(vqsr) with BaseCommandArguments {
    this.mode = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
    this.analysisName = "CEUTrio_AVQSRi"
    this.ts_filter_level = 99.0
  }

  def filterSNPsAndIndels(snpsAndIndelsVCF: File, hardFilter: Boolean): File = {
    val selectSNPs = new SelectVariants with BaseCommandArguments
    selectSNPs.V = snpsAndIndelsVCF
    selectSNPs.selectType = List(VariantContext.Type.SNP)
    selectSNPs.out = swapExt(outputDir, snpsAndIndelsVCF, ".vcf", ".snps.vcf")

    val selectIndels = new SelectVariants with BaseCommandArguments
    selectIndels.V = snpsAndIndelsVCF
    selectIndels.selectType = List(VariantContext.Type.INDEL, VariantContext.Type.MIXED, VariantContext.Type.MNP, VariantContext.Type.SYMBOLIC)
    selectIndels.out = swapExt(outputDir, snpsAndIndelsVCF, ".vcf", ".indels.vcf")

    add(selectSNPs, selectIndels)

    val filteredSNPs: File =
      if ( ! hardFilter ) {
        val snpRecalibrator = new snpRecal(selectSNPs.out)
        val snpApplyVQSR = new applySnpVQSR(snpRecalibrator)
        add(snpRecalibrator, snpApplyVQSR)
        snpApplyVQSR.out
      } else {
        filterVCF(selectSNPs.out, true)
      }

    val filteredIndels =
      if ( !hardFilter && recalIndels ) {
        val indelRecalibrator = new indelRecal(selectIndels.out)
        val indelApplyVQSR = new applyIndelVQSR(indelRecalibrator)
        add(indelRecalibrator, indelApplyVQSR)
        indelApplyVQSR.out
      } else {
        filterVCF(selectIndels.out, false)
      }

    val combineSNPsIndels = new CombineSNPsIndels(filteredSNPs, filteredIndels)
    combineSNPsIndels.out = swapExt(outputDir, snpsAndIndelsVCF, ".unfiltered.vcf", ".filtered.vcf")
    add(combineSNPsIndels)

    combineSNPsIndels.out
  }

  def filterVCF(vcf: File, useSNPFilters: Boolean): File = {
    val vf = new VariantFiltration with BaseCommandArguments
    vf.variant = vcf
    if ( useSNPFilters ) {
      vf.filterName = Seq("BAD_QD", "BAD_FS")
      vf.filterExpression = Seq("QD<2.0", "FS>60.0")
    } else {
      vf.filterName = Seq("BAD_FS")
      vf.filterExpression = Seq("FS>60.0")
    }
    vf.out = swapExt(outputDir, vcf, ".unfiltered.vcf", ".filtered.vcf")
    add(vf)

    vf.out
  }

  class CombineSNPsIndels(snpVCF:File, indelVCF:File) extends CombineVariants with BaseCommandArguments {
    this.variant :+= TaggedFile(new File(snpVCF), "snps")
    this.variant :+= TaggedFile(new File(indelVCF), "indels")
    this.filteredrecordsmergetype = FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    this.assumeIdenticalSamples = true
  }

  //------------------------------------------------------------------------------------ //
  //                      9) Assessing against NA12878 KB                                //
  //------------------------------------------------------------------------------------ //

  /****************************************************************************************
    *                script                                                                 *
    *****************************************************************************************/
  val callSets: scala.collection.mutable.Map[String, File] = scala.collection.mutable.Map[String,File]()

  def callSamples(name: String, hardFilter: Boolean, useHC: Boolean, BAMs: File*): File = {
    val callset = {
      if ( useHC ) {
        val hc = new MyHaplotypeCaller(name)
        for ( bam <- BAMs ) hc.input_file :+= bam
        add(hc)
        hc.out
      } else {
        val ug = new MyUnifiedGenotyper(name)
        for ( bam <- BAMs ) ug.input_file :+= bam
        add(ug)
        ug.out
      }
    }

    val finalCallset: File = filterSNPsAndIndels(callset, hardFilter)
    callSets += name -> finalCallset
    finalCallset
  }

  def script() {
    nSamples = io.Source.fromFile(jointBams).getLines().size + 1

    val na12878alone_HC = callSamples("NA12878Alone.HC", true, true, NA12878ReducedBam)
    val na12878alone_UG = callSamples("NA12878Alone.UG", true, false, NA12878ReducedBam)
    val na12878Fullalone_HC = callSamples("NA12878Alone.full.HC", true, true, NA12878FullBam)
    val na12878Fullalone_UG = callSamples("NA12878Alone.full.UG", true, false, NA12878FullBam)
    val clinicalNa12878alone_UG = callSamples("NA12878Alone.clinical.UG", true, false, ClinicalNA12878Bam)
    val clinicalNa12878alone_HC = callSamples("NA12878Alone.clinical.HC", true, true, ClinicalNA12878Bam)
    val fullJointCall = callSamples("fullJointCall", false, false, jointBams, NA12878ReducedBam)

    //val jointWithoutNA12878 = callSamples("jointWithoutNA12878", false, false, jointBams)

    for ( callsetToAssess <- List(fullJointCall, na12878alone_UG, na12878alone_HC,
      na12878Fullalone_HC, na12878Fullalone_UG, clinicalNa12878alone_HC, clinicalNa12878alone_UG)) {
      val assessKB = new AssessNA12878 with BaseCommandArguments
      assessKB.V :+= callsetToAssess
      assessKB.out = swapExt(outputDir, callsetToAssess, ".vcf", ".kb.gatkreport")
      assessKB.BAM = NA12878ReducedBam
      assessKB.detailed = true
      assessKB.badSites = swapExt(outputDir, callsetToAssess, ".vcf", ".kb.badsites.vcf")
      assessKB.ip = 0
      add(assessKB)
    }

    // write out a combined KB assessment
    def makeAssessment(ip: Int, ext: String) = {
      val assessKB = new AssessNA12878 with BaseCommandArguments
      assessKB.V :+= new TaggedFile(na12878alone_UG, "NA12878.called.alone.ug")
      assessKB.V :+= new TaggedFile(na12878alone_HC, "NA12878.called.alone.hc")
      assessKB.V :+= new TaggedFile(clinicalNa12878alone_HC, "NA12878.called.alone.clinical.hc")
      assessKB.V :+= new TaggedFile(clinicalNa12878alone_UG, "NA12878.called.alone.clinical.ug")
      assessKB.V :+= new TaggedFile(na12878Fullalone_HC, "NA12878.called.alone.full.hc")
      assessKB.V :+= new TaggedFile(na12878Fullalone_UG, "NA12878.called.alone.full.ug")
      assessKB.V :+= new TaggedFile(fullJointCall, "NA12878.called.jointly")
      assessKB.out = qscript.outputDir + "/NA12878.kb.assessment." + ext + ".gatkreport"
      assessKB.ip = ip
      add(assessKB)
    }
    makeAssessment(50, "ip50")
    makeAssessment(0, "ip0")

//    // write out a combined VCF for easy comparisons
//    val cv = new CombineVariants with BaseCommandArguments
//    cv.V = assessKB.V
//    cv.out = qscript.outputDir + "/NA12878.combined.vcf"
//    add(cv)
  }
}

