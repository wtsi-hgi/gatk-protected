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

import net.sf.picard.reference.FastaSequenceFile
import org.broadinstitute.sting.gatk.downsampling.DownsampleType
import org.broadinstitute.sting.pipeline.Pipeline
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.yaml.YamlUtils

class Phase1Calling extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="the chromosome to process", shortName="chr", required=false)
  var chr: Int = 20

  @Input(doc="output path", shortName="outputDir", required=false)
  var outputDir: String = "/humgen/1kg/processing/allPopulations_chr20_phase1_release/perPop.cleaned.BAQed.bams"

  @Input(doc="base output filename", shortName="baseName", required=false)
  var baseName: String = ""

  @Input(doc="path to tmp space for storing intermediate bam files", shortName="outputTmpDir", required=false)
  var outputTmpDir: String = "/humgen/1kg/processing/allPopulations_chr20_phase1_release/perPop.cleaned.BAQed.bams"

  private val tmpDir: File = new File("/broad/shptmp/rpoplin/tmp/")
  private val reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.rod")
  private val dindelPilotCalls: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg.pilot_release.merged.indels.sites.hg19.vcf"
  private val dindelAFRCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110111_august_dindel_indel_calls/AFR.dindel_august_release.20110110.sites.vcf.gz"
  private val dindelASNCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110111_august_dindel_indel_calls/ASN.dindel_august_release.20110110.sites.vcf.gz"
  private val dindelEURCalls: String = "/humgen/1kg/DCC/ftp/technical/working/20110111_august_dindel_indel_calls/EUR.dindel_august_release.20110110.sites.vcf.gz"
  private val dindelMask: String = "/humgen/1kg/processing/allPopulations_wholeGenome_august_release/pilot1.dindel.mask.bed"
  val hapmap = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
  val g1k = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg_pilot1_projectCalls/ALL.low_coverage.2010_07.hg19.vcf"
  val omni = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/764samples.deduped.b37.annot.vcf"
  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
  val populations = List("ASW","CEU","CHB","CHS","CLM","FIN","GBR","JPT","LWK","MXL","PUR","TSI","YRI")    
  //val populations = List("JPT","ASN","AMR")
  //val populations = List("EUR","AMR","ASN","AFR")
  //val populations = List("FIN", "LWK")
  private val intervals: String = "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals"
  //val populations = List("ZZZ") // small set used for debugging

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = 3
    this.jobTempDir = qscript.tmpDir
    this.DBSNP = qscript.dbSNP
  }

  def script = {
    callThisChunk() // using scatter/gather capabilities of Queue so no need to for loop over 1Mb chunks of the chromosome
  }

  def callThisChunk() = {

    val interval = "%d".format(qscript.chr)
    for( population <- qscript.populations ) {
      val baseName: String = qscript.outputDir + "/" + population + ".phase1.chr" + qscript.chr.toString
      var bamList: File = new File("/humgen/1kg/processing/allPopulations_chr20_phase1_release/perPop.cleaned.BAQed.bams/%s.phase1.chr%d.cleaned.bam".format(population, qscript.chr))
      if( population == "ASN" || population == "EUR" || population == "AFR" || population == "AMR" ) {
        bamList = new File("/humgen/1kg/processing/allPopulations_chr20_phase1_release/perPop.cleaned.BAQed.bams/%s.chr%d.cleaned.list".format(population, qscript.chr))
      }

      val rawCalls = new File(baseName + ".raw.vcf")
      val filteredCalls = new File(baseName + ".filtered.vcf")
      val clusterFile = new File(baseName + ".omni.clusters")
      val recalibratedCalls = new File(baseName + ".recal.vcf")
      val tranchesFile = new File(baseName + ".ts.omni.tranches")

      var call = new UnifiedGenotyper with CommandLineGATKArgs
      call.intervalsString ++= List(qscript.intervals)
      call.scatterCount = 63 // the smallest interval list has 63 intervals, one for each Mb on chr20
      call.dcov = 50
      call.stand_call_conf = 4.0
      call.stand_emit_conf = 4.0
      call.input_file :+= bamList
      call.out = rawCalls
      call.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
      call.analysisName = baseName + "_UG"

      var filter = new VariantFiltration with CommandLineGATKArgs
      filter.intervalsString ++= List(qscript.intervals)
      filter.scatterCount = 10
      filter.variantVCF = rawCalls
      filter.out = filteredCalls
      filter.filterName ++= List("HARD_TO_VALIDATE")
      filter.filterExpression ++= List("\"MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1\"")
      filter.analysisName = baseName + "_VF"
      //filter.rodBind :+= RodBind("mask", "Bed", qscript.dindelMask)
      //filter.maskName = "InDel"

      var gvc = new GenerateVariantClusters with CommandLineGATKArgs
      gvc.rodBind :+= RodBind("hapmap", "VCF", qscript.hapmap)
      gvc.rodBind :+= RodBind("1kg", "VCF", qscript.omni)
      gvc.rodBind :+= RodBind("input", "VCF", filteredCalls )
      gvc.clusterFile = clusterFile
      gvc.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
      gvc.analysisName = baseName + "_GVC"
      gvc.intervalsString ++= List(qscript.intervals)
      gvc.qual = 100 // clustering parameters to be updated soon pending new experimentation results
      gvc.std = 4.5
      gvc.mG = 6

      var vr = new VariantRecalibrator with CommandLineGATKArgs
      vr.rodBind :+= RodBind("1kg", "VCF", qscript.omni)
      vr.rodBind :+= RodBind("hapmap", "VCF", qscript.hapmap)
      vr.rodBind :+= RodBind("truthOmni", "VCF", qscript.omni)
      vr.rodBind :+= RodBind("truthHapMap", "VCF", qscript.hapmap)
      vr.rodBind :+= RodBind("input", "VCF", filteredCalls )
      vr.clusterFile = clusterFile
      vr.analysisName = baseName + "_VR"
      vr.intervalsString ++= List(qscript.intervals)
      vr.ignoreFilter ++= List("HARD_TO_VALIDATE")
      vr.target_titv = 2.3
      vr.sm = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibrator.SelectionMetricType.TRUTH_SENSITIVITY
      vr.tranche ++= List("0.1", "1.0", "2.0", "3.0", "5.0", "10.0", "100.0")
      vr.out = recalibratedCalls
      vr.priorDBSNP = 10.0
      vr.priorHapMap = 12.0
      vr.prior1KG = 12.0
      vr.tranchesFile = tranchesFile      

      add(call, filter, gvc, vr)
    }

  }
}