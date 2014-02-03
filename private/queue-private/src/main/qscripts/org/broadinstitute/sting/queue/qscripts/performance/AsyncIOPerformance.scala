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

package org.broadinstitute.sting.queue.qscripts.performance

import _root_.scala._
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.gatk.RodBind._
import org.broadinstitute.sting.queue.function.ListWriterFunction
import scala.io._
import org.broadinstitute.sting.queue.QScript

object IOType extends Enumeration {
  type IOType = Value
  val SYNC = Value("sync")
  val ASYNC = Value("async")
}
import IOType._

class AsyncIOPerformance extends QScript {
  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "bams", doc="BAMs", required=false)
  var bamListFile: File = new File("/humgen/gsa-pipeline/PA56S/T2DGO_Freeze2/wg/v3/T2DGO_Freeze2.bam.list");

  @Argument(fullName = "num_bams", shortName = "nb", doc="How many BAMs should be processed in one go.", required = false)
  var numBams: List[Int] = List(-1)

  @Argument(shortName = "logging_directory", doc="The target directory for log files", required=false)
  var projectDirectory: File = new File("/humgen/gsa-scr1/hanna/asyncio")

  @Argument(shortName = "dbsnp", doc="dbSNP file to use with RTC/IR",required = false)
  val dbsnp = "/humgen/gsa-pipeline/resources/b37/v2/dbsnp_132.b37.vcf"

  @Argument(shortName = "indels", doc="1kG indels file to use with IR",required = false)
  val indels = "/humgen/gsa-pipeline/resources/b37/v2/1000G_indels_for_realignment.b37.vcf"

  // Default range should be over chr20: first 10Mb, 10M-11M, near the centromere
  @Argument(shortName = "L",doc = "Subset of intervals to use",required = false)
  var intervals: List[String] = List("20:1-10000000","20:10000000-11000000","20:26000000-27000000")

  @Argument(shortName="nt", doc="Number of CPU threads.  Applicable to the UG joint test only.", required=false)
  var numThreads: List[Int] = Nil

  @Argument(fullName = "read_buffer_size", shortName = "rbs", doc="Number of reads for read shard to hold in memory.  Applicable to Indel Realigner only.", required = false)
  var readBufferSizes: List[Int] = List(-1)



  val outputDirectory = new File(projectDirectory,"output");
  val logDirectory = new File(projectDirectory,"logs")

  trait UNIVERSAL_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = referenceFile;
    this.dcov = 50
    this.memoryLimit = 4
  }

  def buildCallingJob(bamList: File, interval: String, synchronicity: IOType, uniquifier: String, numCPUThreads: Int = 1) : UnifiedGenotyper = {
    var caller = new UnifiedGenotyper with UNIVERSAL_ARGS
    caller.input_file = List(bamList)
    if(interval != null)
      caller.intervalsString = List(interval)
    caller.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP

    var outputBase = "ug." + synchronicity + ".calls." + uniquifier;
    if(numCPUThreads != 1) {
      caller.num_threads = numCPUThreads
      outputBase += ".nt" + numCPUThreads
    }

    caller.out = new File(outputDirectory,outputBase + ".vcf")
    caller.jobOutputFile = new File(logDirectory,outputBase + ".log")

    if(synchronicity == IOType.ASYNC)
      makeJobIOThreaded(caller)

    caller
  }

  def buildRTCJob(bamList: File, interval: String, synchronicity: IOType): RealignerTargetCreator =  {
    val rtc = new RealignerTargetCreator with UNIVERSAL_ARGS
    rtc.input_file = List(bamList)
    rtc.intervalsString = List(interval)
    rtc.mismatchFraction = 0.0
    rtc.known :+= dbsnp
    rtc.known :+= indels
    rtc.isIntermediate = true
    rtc.out = new File(projectDirectory,synchronicity + ".target." + interval + ".intervals")
    rtc.jobOutputFile = new File(logDirectory,"rtc." + synchronicity + ".joint." + interval + ".log")

    if(synchronicity == IOType.ASYNC) {
      makeJobIOThreaded(rtc)
    }

    rtc
  }

  def buildCleaningJob(bamList: File, numBams: Int, interval: String, dirtyIntervals: File, readBufferSize: Int, synchronicity: IOType): IndelRealigner =  {
    val clean = new IndelRealigner with UNIVERSAL_ARGS
    clean.input_file = List(bamList)
    clean.intervalsString = List(interval)
    clean.targetIntervals = dirtyIntervals
    clean.read_buffer_size = readBufferSize
    clean.known :+= dbsnp
    clean.known :+= indels
    clean.consensusDeterminationModel = org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel.USE_READS
    clean.simplifyBAM = true
    clean.bam_compression = 1
    var baseName: String = ""

    if(readBufferSize != -1)
      baseName = "clean." + synchronicity + ".joint." + interval + ".rbs" + readBufferSize
    else
      baseName = "clean." + synchronicity + ".joint." + interval
    if(numBams != -1)
      baseName += ".nb" + numBams

    clean.out = new File(outputDirectory,baseName + ".bam")
    clean.jobOutputFile = new File(logDirectory,baseName + ".log")

    if(synchronicity == IOType.ASYNC) {
      makeJobIOThreaded(clean)
    }

    clean
  }

  def makeJobIOThreaded(job: CommandLineGATK) {
    // add in the num threads and num io threads parameters.
    if(job.num_threads != None)
      job.num_threads += 2
    else
      job.num_threads = 2
    job.num_io_threads = 2
  }

  def script = {
    if(!projectDirectory.exists())
      projectDirectory.mkdirs();

    if(!outputDirectory.exists())
      outputDirectory.mkdir();

    if(!logDirectory.exists())
      logDirectory.mkdir();


    // Read in the list of BAM files, filtering out files that aren't present.
    var bamFiles = List.empty[File];
    for(line <- Source.fromFile(bamListFile).getLines) {
      val bamFile = new File(line);
      if(bamFile.exists()) {
        bamFiles :+= bamFile;
        if(bamFiles.size % 10 == 0)
          printf("PROGRESS: %d BAM files found%n",bamFiles.size);
      }
      else
        printf("Skipping BAM file %s because it doesn't exist on disk%n",bamFile.getAbsolutePath);
    }
    printf("A total of %d BAM files are prepared for processing%n",bamFiles.size);

    // Write a list file of processed BAMs.
    val writeBamList = new ListWriterFunction
    writeBamList.inputFiles = bamFiles
    writeBamList.listFile = new File(projectDirectory,"present_bams.list")
    writeBamList.run

    // Call everything without filters on each file individually.
    for(bamFile <- bamFiles) {
      var caller = buildCallingJob(bamFile,null,IOType.ASYNC,"single."+bamFile.getName)
      add(caller)
    }

    // Call everything on the entire fileset, both tagged and untagged.
    // Run RTC / IR on the entire merged dataset
    for(interval <- intervals) {
      // Add the single-CPU thread calling job
      add(buildCallingJob(writeBamList.listFile,interval,IOType.SYNC,"joint."+interval))
      add(buildCallingJob(writeBamList.listFile,interval,IOType.ASYNC,"joint."+interval))

      // If specified, add the multi-CPU threaded job
      for(threadCount <- numThreads) {
        add(buildCallingJob(writeBamList.listFile,interval,IOType.SYNC,"joint."+interval,threadCount))
        add(buildCallingJob(writeBamList.listFile,interval,IOType.ASYNC,"joint."+interval,threadCount))
      }

      var rtcSync = buildRTCJob(writeBamList.listFile,interval,IOType.SYNC)
      add(rtcSync)
      var rtcAsync = buildRTCJob(writeBamList.listFile,interval,IOType.ASYNC)
      add(rtcAsync)

      for(bamCount <- numBams; readBufferSize <- readBufferSizes) {
        println("bamCount = " + bamCount + ", readBufferSize = " + readBufferSize)
        var bamWorkingSet = writeBamList.listFile
        var workingSetSize = bamFiles.size

        if(bamCount != -1) {
          val writeBamListSubset = new ListWriterFunction
          writeBamListSubset.inputFiles = bamFiles.take(bamCount)
          writeBamListSubset.listFile = new File(projectDirectory,"present_bams."+ bamCount + ".list")
          writeBamListSubset.run

          bamWorkingSet = writeBamListSubset.listFile
          workingSetSize = bamCount
        }

        var cleanSync = buildCleaningJob(bamWorkingSet,workingSetSize,interval,rtcSync.out,readBufferSize,IOType.SYNC)
        add(cleanSync)
        var cleanAsync = buildCleaningJob(bamWorkingSet,workingSetSize,interval,rtcAsync.out,readBufferSize,IOType.ASYNC)
        add(cleanAsync)
      }
    }

  }

}

