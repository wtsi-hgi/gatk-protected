/*
 * Copyright (c) 2011, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

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

