/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.queue.qscripts.calling

import htsjdk.samtools.reference.FastaSequenceFile
import org.broadinstitute.gatk.engine.downsampling.DownsampleType
import org.broadinstitute.gatk.utils.pipeline.Pipeline
import org.broadinstitute.gatk.tools.DownsampleType
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.samtools._
import org.broadinstitute.gatk.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.gatk.utils.yaml.YamlUtils

class ASHGcalling extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="the chromosome to process", shortName="chr", required=true)
  var chr: Int = _

  @Input(doc="output path", shortName="outputDir", required=false)
  var outputDir: String = "/humgen/1kg/processing/allPopulations_wholeGenome_august_release/calls/"

  @Input(doc="base output filename", shortName="baseName", required=false)
  var baseName: String = "ALL.august"

  @Input(doc="path to tmp space for storing intermediate bam files", shortName="outputTmpDir", required=false)
  var outputTmpDir: String = "/humgen/gsa-hpprojects/august_cleaned_bams"

  private val tmpDir: File = new File("/broad/shptmp/rpoplin/")
  private val reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  private val dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.rod")
  private val targetIntervals: File = new File("/humgen/1kg/processing/allPopulations_wholeGenome_august_release/knownIndels.intervals")
  private val dindelCalls: String = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg.pilot_release.merged.indels.sites.hg19.vcf"
  private val dindelMask: String = "/humgen/1kg/processing/allPopulations_wholeGenome_august_release/pilot1.dindel.mask.bed"
  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
  val populations = List("YRI","LWK","ASW","PUR","CEU","TSI","GBR","FIN","MXL","CHB","CHS","JPT")

  private var pipeline: Pipeline = _

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = 2
    this.DBSNP = qscript.dbSNP
    this.jobTempDir = qscript.tmpDir
  }

  class SamtoolsBaqFunction extends CommandLineFunction {
    @Input(doc="The input BAM file") var in_bam: File = _
    @Output(doc="The output BAM file") var out_bam: File = _
    def commandLine = "/humgen/gsa-scr1/rpoplin/samtools/samtools calmd -br %s %s > %s".format(in_bam.getAbsolutePath, qscript.reference, out_bam.getAbsolutePath)
  }

  class DeleteMeFunction extends CommandLineFunction {
    @Input(doc="The file to be deleted") var me: File = _
    @Input(doc="The file which must exist before we are allowed to delete") var trigger: File = _
    def commandLine = "rm -f %s".format(me.getAbsolutePath)
  }

  class DeleteMeAllFunction extends CommandLineFunction {
    @Input(doc="The file to be deleted") var me: File = _
    @Input(doc="The file which must exist before we are allowed to delete") var trigger: File = _
    def commandLine = "rm -f %s*".format(me.getAbsolutePath)
  }

  def script = {
    val basesPerJob: Int = 3000000
    val lastBase: Int = qscript.chromosomeLength(qscript.chr - 1)
    var start: Int = 1
    var stop: Int = start - 1 + basesPerJob
    if( stop > lastBase ) { stop = lastBase }
    var jobNumber: Int = 1
    while( jobNumber < (lastBase.toFloat / basesPerJob.toFloat) + 1.0) {
      callThisChunk("%d:%d-%d".format(qscript.chr, start, stop), jobNumber)
      start += basesPerJob
      stop += basesPerJob
      if( stop > lastBase ) { stop = lastBase }
      jobNumber += 1
    }

    /* CombineVariants parses the 800+ genotypes per record and is way too slow. Combine the vcf files together using grep, cat, and sortByRef.pl outside of Queue
    combineVariants = new CombineVariants with CommandLineGATKArgs
    combineVariants.rodBind = vcfChunks
    combineVariants.out = new TaggedFile(qscript.baseName + ".chr" + qscript.chr.toString + ".filtered.vcf", "vcf")
    combineVariants.variantmergeoption = org.broadinstitute.gatk.engine.contexts.variantcontext.VariantContextUtils.VariantMergeType.UNION
    combineVariants.genotypemergeoption = org.broadinstitute.gatk.engine.contexts.variantcontext.VariantContextUtils.GenotypeMergeType.UNSORTED
    combineVariants.setKey = "null"
    add(combineVariants)
    */
  }

  def callThisChunk(interval: String, jobNumber: Int) = {

    val baseName: String = qscript.outputDir + "/chr" + qscript.chr.toString + "/" + qscript.baseName + ".chr" + qscript.chr.toString + "." + jobNumber.toString +"."
    var call = new UnifiedGenotyperV2 with CommandLineGATKArgs
    val rawCalls = new File(baseName + "raw.vcf")

    for( population <- qscript.populations ) {
      val baseTmpName: String = qscript.outputTmpDir + "/chr" + qscript.chr.toString + "/" + population + ".august.chr" + qscript.chr.toString + "." + jobNumber.toString +"."
      val bamList: File = new File("/humgen/1kg/processing/allPopulations_wholeGenome_august_release/bamLists/%s.chr%d.bam.list".format(population, qscript.chr))

      // 1.) Clean at known indels
      var clean = new IndelRealigner with CommandLineGATKArgs
      val cleanedBam = new File(baseTmpName + "cleaned.bam")
      clean.memoryLimit = 4
      clean.input_file :+= bamList
      clean.intervalsString :+= interval
      clean.targetIntervals = qscript.targetIntervals
      clean.out = cleanedBam
      clean.rodBind :+= RodBind("indels", "VCF", qscript.dindelCalls)
      clean.knownsOnly = true
      clean.LOD = 1.0
      clean.sortInCoordinateOrderEvenThoughItIsHighlyUnsafe = true
      clean.compress = 2
      clean.jobName = baseName + population + ".clean"
      //clean.stripBam = true
      //clean.fileSystemUsage = "indium"

      // 2.) Apply BAQ calculation
      var baq = new SamtoolsBaqFunction
      val baqedBam = new File(baseTmpName + "cleaned.baq.bam")
      baq.memoryLimit = 4
      baq.in_bam = cleanedBam
      baq.out_bam = baqedBam
      baq.jobName = baseName + population + ".baq"
      //baq.fileSystemUsage = "iodine"

      // 3a.) Delete cleaned bam
      var deleteClean = new DeleteMeFunction
      deleteClean.me = cleanedBam
      deleteClean.trigger = baqedBam
      deleteClean.jobName = baseName + population + ".deleteClean"
      //deleteClean.fileSystemUsage = "iodine"

      // 3b.) Index BAQ'ed bam
      var index = new SamtoolsIndexFunction
      index.bamFile = baqedBam
      index.jobName = baseName + population + ".index"
      //index.fileSystemUsage = "iodine"

      // 5a.) Delete BAQ'ed bam and index
      //var deleteBaq = new DeleteMeAllFunction
      //deleteBaq.me = baqedBam
      //deleteBaq.trigger = rawCalls
      //deleteBaq.jobName = baseName + population + ".deleteBaq"
      //deleteBaq.fileSystemUsage = "iodine"

      call.input_file :+= baqedBam

      //add(clean, baq, deleteClean, index, deleteBaq)
      add(clean, baq, deleteClean, index)
    }

    // 4.) Call with UGv2
    call.memoryLimit = 4
    call.intervalsString :+= interval
    call.out = rawCalls
    call.dcov = 50
    call.standard_min_confidence_threshold_for_calling = 50
    call.standard_min_confidence_threshold_for_emitting = 30
    call.min_mapping_quality_score = 20
    call.min_base_quality_score = 20
    call.pnrm = org.broadinstitute.gatk.playground.gatk.walkers.genotyper.AlleleFrequencyCalculationModel.Model.GRID_SEARCH
    call.jobName = baseName + "call"
    //call.fileSystemUsage = "iodine"

    // 5b.) Filter near indels and HARD_TO_VALIDATE
    var filter = new VariantFiltration with CommandLineGATKArgs
    val filteredCalls = new File(baseName + "filtered.vcf")
    filter.memoryLimit = 1
    filter.out = filteredCalls
    filter.intervalsString :+= interval
    filter.variantVCF = rawCalls
    filter.rodBind :+= RodBind("mask", "Bed", qscript.dindelMask)
    filter.maskName = "InDel"
    filter.filterName ++= List("HARD_TO_VALIDATE")
    filter.filterExpression ++= List("\"MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1\"")
    filter.jobName = baseName + "filter"
    //filter.fileSystemUsage = "indium"

    // 6.) Delete raw calls and index
    var deleteRawCalls = new DeleteMeAllFunction
    deleteRawCalls.me = rawCalls
    deleteRawCalls.trigger = filteredCalls
    deleteRawCalls.jobName = baseName + "deleteRawCalls"
    //deleteRawCalls.fileSystemUsage = "indium"

    add(call, filter, deleteRawCalls)
  }
}
