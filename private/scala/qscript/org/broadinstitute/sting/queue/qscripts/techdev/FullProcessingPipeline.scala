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

package org.broadinstitute.sting.gatk.walkers.techdev

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard._
import net.sf.samtools.{SAMReadGroupRecord, SAMFileReader, SAMFileHeader}
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode
import collection.JavaConversions._
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException
import net.sf.samtools.SAMFileHeader.SortOrder

class FullProcessingPipeline extends QScript {
  qscript =>

  @Input(doc="input BAM files -- one per -I, no bam lists", fullName="input", shortName="I", required=true)
  var input: Seq[File] = _

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=false)
  var reference = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")

  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=false)
  var dbSNP = Seq(new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_137.b37.vcf"))

  @Argument(doc="Final root name of the BAM file", fullName="project", shortName="p", required=false)
  var project: String = ""

  @Argument(doc="hard clip illumina adapter sequence", fullName="hard_clip_adapter", shortName="adapter", required=false)
  var clipAdapter:Boolean = false

  @Hidden
  @Argument(doc="Skip reverting the bam", fullName = "skip_revert", shortName = "skip_revert", required=false)
  var skip_revert: Boolean = false

  @Hidden
  @Argument(doc="only process the minimum possible to get an analizable bam", shortName="fly", required=false)
  var flyThrough = false

  @Hidden
  @Argument(doc="Number of threads to use with BWA", shortName="t", required=false)
  var threads = 8

  @Hidden
  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = -1

  @Hidden
  @Argument(doc="Define the default platform for Count Covariates -- useful for techdev purposes only.", fullName="default_platform", shortName="dp", required=false)
  var defaultPlatform: String = ""



  val cleaningModel: ConsensusDeterminationModel = ConsensusDeterminationModel.USE_READS

  def script() {
    // keep a record of the number of contigs in the first bam file in the list
    if (qscript.nContigs < 0)
      qscript.nContigs = QScriptUtils.getNumberOfContigs(qscript.input(0))

    // assert that all bams are from the same sample.
    val sampleName = getSampleName(input)

    // BAM files generated by the pipeline
    val bams = align(convert_to_fastq(revert_bam(input)))
    val sortedBam  = new File(qscript.project + sampleName + ".sorted.bam")
    val cleanedBam = swapExt(sortedBam, ".bam", ".clean.bam")
    val recalBam   = swapExt(cleanedBam, ".bam", ".recal.bam")
    val reducedBam = swapExt(recalBam, ".bam", ".reduced.bam")

    // Accessory files
    val targetIntervals = swapExt(sampleName, ".bam", ".intervals")
    val recalFile       = swapExt(sampleName, ".bam", ".grp")


    add(
      merge(bams, sortedBam),
      target(sortedBam, targetIntervals),
      clean(sortedBam, targetIntervals, cleanedBam),
      bqsr(cleanedBam, recalFile),
      printreads(cleanedBam, recalFile, recalBam),
      reduce(recalBam, reducedBam)
    )
  }

  /****************************************************************************
    * Helper methods
    ****************************************************************************/
//  def splitIfMultipleReadGroups(bam: File): Seq[File] = {
//    val samReader = new SAMFileReader(bam)
//    val readGroups: List[SAMReadGroupRecord] = samReader.getFileHeader.getReadGroups
//    var outBAMs : Seq[File] = Seq(bam)
//    if (QScriptUtils.hasMultipleSamples(readGroups)) {
//      add(split(bam)) // split the multiple RGs into multiple files
//      outBAMs = Seq()
//      for (rg: SAMReadGroupRecord <- readGroups) {
//        val filename = rg.getSample + "_" + rg.getId() + ".bam"
//        outBAMs :+= new File(filename)
//      }
//    }
//    return outBAMs
//  }


  /**
   * Reverts all input bam files using Picard's RevertSAM
   *
   * @param bams input bam files
   * @return a sequence of reverted bam files
   */
  def revert_bam(bams:Seq[File]) : Seq[(File, String)] = {
    var revBams: Seq[(File, String)] = Seq()
    for (bam <- bams) {
      if (skip_revert) {
        revBams :+= (bam, getBWAReadGroupLine(bam))
      }
      else {
        val revertedBAM = swapExt(bam, ".bam", ".reverted.bam")
        add(revert(bam, revertedBAM))
        revBams :+= (revertedBAM, getBWAReadGroupLine(bam))
      }
    }
    revBams
  }

  /**
   * converts bam files to interleaved fastq clipping adapters if necessary
   *
   * @param bams the input bams (should be reverted)
   * @return the sequence of interleaved fastq files and their read group strings in a tuple
   */
  def convert_to_fastq(bams:Seq[(File, String)]) : Seq[(File, String)] = {
    var fqs: Seq[(File, String)] = Seq()
    for ( (bam, rg) <- bams) {
      if (clipAdapter) {
        fqs :+= (picard_fq(bam), rg)
      }
      else {
        fqs :+= (htslib_fq(bam), rg)
      }
    }
    fqs
  }

  /**
   * alignment of interleaved fastq files
   *
   * Aligns all interleaved fastqs
   *
   * @param fqs a sequence of tuples (fastq file, read group string)
   * @return list of aligned sam files
   */
  def align(fqs:Seq[(File, String)]) : Seq[File] = {
    var sams: Seq[File] = Seq()
    for ((fq, rg) <- fqs) {
      sams :+= align_fastq(fq, rg)
    }
    sams
  }


  /**
   * convert to fastq without clipping adapters (using htslib)
   *
   * @param bam input bam
   * @return output interleaved fastq
   */
  def htslib_fq(bam: File): File = {
    val fq         = swapExt(bam, ".bam", ".fq.gz")
    add(bam2fq(bam, fq))
    fq
  }

  /**
   * convert to fastq and clip adapter sequence (using picard)
   *
   * @param bam input bam
   * @return output interleaved fastq
   */
  def picard_fq(bam: File): File = {
    val queryBAM = swapExt(bam,".bam",".query.bam")
    val markedBAM = swapExt(bam,".bam",".adaptor_marked.bam")
    val fq = swapExt(bam,".bam",".clipped.fq")

    add (
      sortSam(bam, queryBAM, SortOrder.queryname),
      mark_adaptor(queryBAM, markedBAM),
      sam2fq(markedBAM, fq)
    )
    fq
  }

  /**
   * aligns an interleaved fastq file
   *
   * @param fq the fastq file
   * @param rg the read group string
   * @return the aligned sam file
   */
  def align_fastq(fq: File, rg: String): File = {
    val alignedSam = swapExt(fq, ".fq.gz", ".aligned.sam")
    add(bwa(fq, rg, alignedSam, qscript.threads))
    alignedSam
  }

  /**
   * parses the read group line from a bam file and returns the string format used by BWA
   *
   * @param bam the input bam file
   * @return the string used by BWA to reconstitute the read group information
   */
  def getBWAReadGroupLine (bam:File): String = {
    for (rg: SAMReadGroupRecord <- getReadGroupList(bam)) {  // freaking hate java / scala collection conversions!
      val id = "@RG\tID:" + rg.getId
      val lb = if (rg.getLibrary != null) {"\\tLB:" + rg.getLibrary} else {""}
      val pl = if (rg.getPlatform != null) {"\\tPL:" + rg.getPlatform} else {""}
      val pu = if (rg.getPlatformUnit != null) {"\\tPU:" + rg.getPlatformUnit} else {""}
      val sm = if (rg.getSample != null) {"\\tSM:" + rg.getSample} else {""}
      val cn = if (rg.getSequencingCenter != null) {"\\tCN:" + rg.getSequencingCenter} else {""}
      val ds = if (rg.getDescription != null) {"\\tDS:" + rg.getDescription} else {""}

      return "'" + id + lb + pl + pu + sm + cn + ds + "'" // we only do 1 readgroup per bam file
    }
    ""
  }

  /**
   * Extracts the sample name from all BAM files and checks that they all match
   * (in this pipeline, all bams must be from the same sample -- run multiple samples in parallel processess)
   *
   * @param bams all input bams
   * @return the sample name
   * @throws ReviewedStingException if the sample names are not the same across all files
   */
  def getSampleName (bams:Seq[File]): String = {
    val sample = getReadGroupList(bams(0))(0).getSample  // first read group of the first bam file in the list -- happy to blow up if this fails.
    for (bam <- bams) {
      for (rg: SAMReadGroupRecord <- getReadGroupList(bam)) {
        if (rg.getSample != sample)
          throw new ReviewedStingException(String.format("BAM file %s contains multiple samples, this pipeline cannot handle that", bam.getName))
      }
    }
    sample
  }

  /**
   * Gets the list of read groups from a BAM file (converts it neatly to a scala Sequence)
   *
   * @param bam the input bam file
   * @return a sequence of read group records
   */
  def getReadGroupList(bam:File): Seq[SAMReadGroupRecord] = {
    val samReader = new SAMFileReader(bam)
    val header = samReader.getFileHeader
    header.getReadGroups
  }

  /****************************************************************************
    * Classes (GATK Walkers)
    ****************************************************************************/

  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 4
    this.isIntermediate = true
  }

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
  }

  trait SAMargs extends PicardBamFunction with ExternalCommonArgs {
    this.maxRecordsInRam = 100000
  }

  case class bam2fq (@Input inBam: File, @Output outFQ: File) extends org.broadinstitute.sting.queue.function.CommandLineFunction with ExternalCommonArgs {
    def commandLine = "htscmd bamshuf -uOn 128 " + inBam + " " + outFQ + ".tmp" + " | htscmd bam2fq -a - | gzip > " + outFQ
    this.memoryLimit = 8
    this.analysisName = "FastQ"
  }

  case class bwa (@Input inFQ: File, inRG: String, @Output outSAM: File, threads: Int) extends org.broadinstitute.sting.queue.function.CommandLineFunction with ExternalCommonArgs{
    def commandLine = "bwa mem -p -M -t " + threads + " -R " + inRG + " " + reference + " " + inFQ + " > " + outSAM
    this.memoryLimit = 8
    this.nCoresRequest = threads
    this.analysisName = "BWA"
  }

  case class mark_adaptor(@Input inBAM: File, @Output outBAM: File) extends org.broadinstitute.sting.queue.function.CommandLineFunction with ExternalCommonArgs {
    def commandLine = "java -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar /seq/software/picard/current/bin/MarkIlluminaAdapters.jar INPUT=" + inBAM + " OUTPUT=" + outBAM + " PE=true ADAPTERS=PAIRED_END M=" + outBAM + ".adapter_metrics"
    this.memoryLimit = 4
    this.analysisName = "ADAPTER"
  }

  case class revert (inBam: File, outBam: File) extends RevertSam with SAMargs {
    this.output = outBam
    this.input :+= inBam
  }

  case class sam2fq (inBAM: File, outFQ: File) extends SamToFastq with SAMargs {
    this.input :+= inBAM
    this.fastq = outFQ
    this.interleave = true
    this.clippingAction = "X"
    this.clippingAttribute = "XT"
  }

  case class sortSam (inBam: File, outBam: File, sortOrderP: SortOrder) extends SortSam with SAMargs {
    this.input :+= inBam
    this.output = outBam
    this.sortOrder = sortOrderP
  }

  case class merge (inBams: Seq[File], outBam: File) extends MergeSamFiles with SAMargs {
    this.input = inBams
    this.output = outBam
    this.sortOrder = SAMFileHeader.SortOrder.coordinate
  }

  case class target (inBams: File, outIntervals: File) extends RealignerTargetCreator with CommandLineGATKArgs {
    this.input_file :+= inBams
    this.out = outIntervals
    this.mismatchFraction = 0.0
    this.known ++= qscript.dbSNP
    this.scatterCount = nContigs
  }

  case class clean (inBams: File, tIntervals: File, outBam: File) extends IndelRealigner with CommandLineGATKArgs {
    this.input_file :+= inBams
    this.targetIntervals = tIntervals
    this.out = outBam
    this.known ++= qscript.dbSNP
    this.consensusDeterminationModel = qscript.cleaningModel
    this.compress = 0
    this.read_filter :+= "BadCigar"
    this.scatterCount = nContigs
  }

  case class bqsr (inBam: File, outRecalFile: File) extends BaseRecalibrator with CommandLineGATKArgs {
    this.knownSites ++= qscript.dbSNP
    this.covariate ++= Seq("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "ContextCovariate")
    this.input_file :+= inBam
    this.out = outRecalFile
    if (!defaultPlatform.isEmpty) this.default_platform = defaultPlatform
    if (qscript.flyThrough) {this.intervalsString :+= "1"}
    this.scatterCount = nContigs
  }

  case class printreads (inBam: File, inRecalFile: File, outBam: File) extends PrintReads with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.BQSR = inRecalFile
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.out = outBam
    this.scatterCount = nContigs
    this.isIntermediate = false
  }

  case class reduce (inBam: File, outBam: File) extends ReduceReads with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.out = outBam
    this.scatterCount = nContigs
    this.isIntermediate = false
  }
}
