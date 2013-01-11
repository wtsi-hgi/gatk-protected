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

package org.broadinstitute.sting.queue.qscripts.variantQC

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.QScript
import org.apache.commons.io.FilenameUtils;

class distributedGATKPerformance extends QScript {
  qscript =>

  @Argument(shortName="gatk", doc="gatk jar file", required=true)
  var gatkJarFile: File = _

  @Argument(shortName="outputDir", doc="output directory", required=false)
  var outputDir: String = ""

  @Argument(shortName="dataset", doc="selects the datasets to run. If not provided, all datasets will be used", required=false)
  var datasets: List[String] = Nil

  @Argument(shortName="waysParallel", doc="selects the datasets to run. If not provided, all datasets will be used", required=false)
  var waysParallelArg: List[Int] = Nil

  @Argument(shortName="long", doc="runs long calculations", required=false)
  var long: Boolean = false

  @Argument(shortName="test", doc="runs long calculations", required=false)
  var test: Boolean = false

  @Argument(shortName="limitTo30Min", doc="runs long calculations", required=false)
  var limitTo30Min: Boolean = false

  @Argument(shortName="huge", doc="runs long calculations", required=false)
  var huge: Int = -1

  @Argument(shortName="justDist", doc="runs long calculations", required=false)
  var justDist: Boolean = false

  @Argument(shortName="justSG", doc="runs long calculations", required=false)
  var justSG: Boolean = false

  @Argument(shortName="trackerDir", doc="root directory for distributed tracker files", required=false)
  var trackerDir: String = "" // "/humgen/gsa-scr1/depristo/tmp/"

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK { logging_level = "DEBUG"; jarFile = gatkJarFile; memoryLimit = 2; }

  class Target(
          val baseName: String,
          val reference: File,
          val dbsnpFile: String,
          val hapmapFile: String,
          val maskFile: String,
          val bamList: File,
          val goldStandard_VCF: File,
          val intervals: String,
          val titvTarget: Double,
          val isLowpass: Boolean,
          val useBAQ: Boolean) {
    val name = qscript.outputDir + baseName
    val clusterFile = new File(name + ".clusters")
    def rawVCF(part: String) = new File(name + "." + part + ".raw.vcf")
    val filteredVCF = new File(name + ".filtered.vcf")
    val titvRecalibratedVCF = new File(name + ".titv.recalibrated.vcf")
    val tsRecalibratedVCF = new File(name + ".ts.recalibrated.vcf")
    val goldStandardName = qscript.outputDir + "goldStandard/" + baseName
    val goldStandardClusterFile = new File(goldStandardName + ".clusters")
  }

  val hg18 = new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
  val b36 = new File("/humgen/1kg/reference/human_b36_both.fasta")
  val b37 = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  val dbSNP_hg18 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_130_hg18.rod"
  val dbSNP_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_130_b36.rod"
  val dbSNP_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_132_b37.leftAligned.vcf"
  val hapmap_hg18 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.hg18_fwd.vcf"
  val hapmap_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b36_fwd.vcf"
  val hapmap_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
  val indelMask_b36 = "/humgen/1kg/processing/pipeline_test_bams/pilot1.dindel.mask.b36.bed"
  val indelMask_b37 = "/humgen/1kg/processing/pipeline_test_bams/pilot1.dindel.mask.b37.bed"

  // ToDos:
  // reduce the scope of the datasets so the script is more nimble
  // figure out how to give names to all the Queue-LSF logs (other than Q-1931@node1434-24.out) so that it is easier to find logs for certain steps
  // create gold standard BAQ'd bam files, no reason to always do it on the fly

  // Analysis to add at the end of the script:
  // auto generation of the cluster plots
  // spike in NA12878 to the exomes and to the lowpass, analysis of how much of her variants are being recovered compared to single sample exome or HiSeq calls
  // produce Kiran's Venn plots based on comparison between new VCF and gold standard produced VCF

  val lowPass: Boolean = true

  val targetDataSets: Map[String, Target] = Map(
    "HiSeq" -> new Target("NA12878.HiSeq", hg18, dbSNP_hg18, hapmap_hg18,
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.indels.10.mask",
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.HiSeq.WGS.bwa.cleaned.recal.bam"),
              new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.vcf"),
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/distributedGATK/whole_genome_chunked.hg18.intervals", 2.07, !lowPass, true),
    "FIN" -> new Target("FIN", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/FIN.79sample.Nov2010.chr20.bam"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),         // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/distributedGATK/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass, true),
    "WEx" -> new Target("NA12878.WEx", hg18, dbSNP_hg18, hapmap_hg18,
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.indels.10.mask",
              new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.WEx.cleaned.recal.bam"),
              new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.vcf"),
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.targets.interval_list", 2.6, !lowPass, true),
    "TGPWExGdA" -> new Target("1000G.WEx.GdA", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/distributedGATK/Barcoded_1000G_WEx_Reduced_Plate_1.20.cleaned.list"),        // BUGBUG: reduce from 60 to 20 people
              new File("/humgen/gsa-scr1/delangel/NewUG/calls/AugustRelease.filtered_Q50_QD5.0_SB0.0.allSamples.SNPs_hg19.WEx_UG_newUG_MQC.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 2.6, !lowPass, true),
    "LowPassN60" -> new Target("lowpass.N60", b36, dbSNP_b36, hapmap_b36, indelMask_b36,
              new File("/humgen/1kg/analysis/bamsForDataProcessingPapers/lowpass_b36/lowpass.chr20.cleaned.matefixed.bam"), // the bam list to call from
              new File("/home/radon01/depristo/work/oneOffProjects/VQSRCutByNRS/lowpass.N60.chr20.filtered.vcf"),           // the gold standard VCF file to run through the VQSR
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.b36.intervals", 2.3, lowPass,true),          // chunked interval list to use with Queue's scatter/gather functionality
    "LowPassAugust" -> new Target("ALL.august.v4", b37, dbSNP_b37, hapmap_b37, indelMask_b37,                               // BUGBUG: kill this, it is too large
              new File("/humgen/1kg/processing/allPopulations_chr20_august_release.cleaned.merged.bams/ALL.cleaned.merged.list"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),
              "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass, true),
    "LowPassEUR363Nov" -> new Target("EUR.nov2010", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
              new File("/humgen/1kg/processing/pipeline_test_bams/EUR.363sample.Nov2010.chr20.bam"),
              new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),         // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
              "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/distributedGATK/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass,false),
    "WExTrio" -> new Target("NA12878Trio.WEx", b37, dbSNP_b37, hapmap_b37, indelMask_b37,
        new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WEx.bwa.cleaned.recal.bams.list"),
        new File("/humgen/gsa-scr1/delangel/NewUG/calls/AugustRelease.filtered_Q50_QD5.0_SB0.0.allSamples.SNPs_hg19.WEx_UG_newUG_MQC.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
        "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 2.6, !lowPass, true)
  )

  def getTargetInterval(target: Target): List[String] = target.name match {
    case "NA12878.HiSeq" =>  List("chr1")
    case "FIN" => List("20")
    case "ALL.august.v4" => List("20")
    case "EUR.nov2010" => List("20")
    case _ => List(target.intervals)
  }

  def script = {

    // Selects the datasets in the -dataset argument and adds them to targets.
    var targets: List[Target] = List()
    if (!datasets.isEmpty)
      for (ds <- datasets)
        targets ::= targetDataSets(ds)                  // Could check if ds was mispelled, but this way an exception will be thrown, maybe it's better this way?
    else                                                // If -dataset is not specified, all datasets are used.
      for (targetDS <- targetDataSets.valuesIterator)   // for Scala 2.7 or older, use targetDataSets.values
        targets ::= targetDS

    val nWays = if ( test ) List(32) else { if ( long ) List(1,2,4,8) else if ( huge != -1 ) List(huge) else List(16,32,64,128) }
    //val nWays = List(2)

    for (target <- targets) {
      for ( scatterP <- if ( test ) List(false) else if ( justSG ) List(true) else if ( justDist ) List(false) else List(true, false) )
        for (nWaysParallel <- nWays ) {
          val aname = "ptype_%s.nways_%d".format(if ( scatterP ) "sg" else "dist", nWaysParallel)

          def addUG(ug: UnifiedGenotyper) = {
            if ( ! long )
              ug.jobLimitSeconds = 60 * 60 * 4
            if ( limitTo30Min )
              ug.jobLimitSeconds = 60 * 30
            add(ug);
          }

          // add scatter/gather or distributed parallelism
          if ( scatterP ) {
            var ug: UnifiedGenotyper = new UnifiedGenotyper(target, aname)
            ug.scatterCount = nWaysParallel
            ug.intervalsString ++= List(target.intervals)
            addUG(ug)
          } else {
            for ( part <- 1 to nWaysParallel) {
              var ug: UnifiedGenotyper = new UnifiedGenotyper(target, aname + ".part" + part)
              ug.intervalsString ++= getTargetInterval(target)
              ug.processingTracker = new File(trackerDir + target.name + "." + aname + ".distributed.txt")
              ug.processingTrackerID = part
              if ( part == 1 )
                ug.performanceLog = new File("%s.%s.pf.log".format(target.name, aname))
              ug.processingTrackerStatusFile = new File("%s.%s.%d.ptstatus.log".format(target.name, aname, part))
              addUG(ug)
            }
          }

        }
    }
  }

  // 1.) Call SNPs with UG
  class UnifiedGenotyper(t: Target, aname: String) extends org.broadinstitute.sting.queue.extensions.gatk.UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.reference_sequence = t.reference
    this.dcov =  if ( t.isLowpass ) { 50 } else { 250 }
    this.stand_call_conf =  if ( t.isLowpass ) { 4.0 } else { 30.0 }
    this.stand_emit_conf =  if ( t.isLowpass ) { 4.0 } else { 30.0 }
    this.input_file :+= t.bamList
    this.out = t.rawVCF(aname)
    this.baq =  if (t.useBAQ) {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE} else {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF}
    this.analysisName = t.name + "_UG." + aname
    if (t.dbsnpFile.endsWith(".rod"))
      this.DBSNP = new File(t.dbsnpFile)
    else if (t.dbsnpFile.endsWith(".vcf"))
      this.rodBind :+= RodBind("dbsnp", "VCF", t.dbsnpFile)
  }
}
