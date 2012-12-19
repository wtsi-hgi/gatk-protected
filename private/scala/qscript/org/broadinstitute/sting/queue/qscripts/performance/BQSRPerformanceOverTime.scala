package org.broadinstitute.sting.queue.qscripts.performance

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.PathUtils
import org.broadinstitute.sting.commandline.ClassType

class BQSRPerformanceOverTime extends QScript {
  @Argument(shortName = "resources", doc="resources", required=true)
  val resourcesDir: String = ""

  @Argument(shortName = "myJarFile", doc="Path to the current GATK jar file", required=true)
  val myJarFile: File = null

  @Argument(shortName = "iterations", doc="it", required=false)
  val iterations: Int = 3

  @Argument(shortName = "ntTest", doc="For each value provided we will use -nt VALUE in the multi-threaded tests", required=false)
  @ClassType(classOf[Int])
  val ntTests: List[Int] = List(1, 2, 3, 4, 6, 8, 10, 12)

  val MY_TAG = "GATKPerformanceOverTime"
  val RECAL_BAM_FILENAME = "CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20GAV.8.bam"
  val RECAL_GATKREPORT_FILENAME = "CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20GAV.8.grp"
  val dbSNP_FILENAME = "dbsnp_132.b37.vcf"
  val b37_FILENAME = "human_g1k_v37.fasta"

  def makeResource(x: String): File = new File("%s/%s".format(resourcesDir, x))

  val GATK_RELEASE_DIR = new File("/humgen/gsa-hpprojects/GATK/bin/")
  val GATKs: Map[String, File] = Map(
    "v2.cur" -> myJarFile, // TODO -- how do I get this value?
    "v2.3" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "2.3"),
    "v2.2" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "2.2"))

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO"
    this.reference_sequence = makeResource(b37_FILENAME)
    this.memoryLimit = 8
    this.intervalsString = List("1")
  }

  def script() {
    // iterate over GATK's and data sets
    for ( iteration <- 0 until iterations ) {
      for ( (gatkName, gatkJar) <- GATKs ) {

        def makeBQSR(): BaseRecalibrator = {
          val BQSR = new BaseRecalibrator() with UNIVERSAL_GATK_ARGS
          BQSR.knownSites :+= makeResource(dbSNP_FILENAME)
          BQSR.input_file :+= makeResource(RECAL_BAM_FILENAME)
          BQSR.out = new File("/dev/null")
          BQSR.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName, "assessment" -> "BQSR"))
          BQSR.jarFile = gatkJar
          BQSR
        }
        addMultiThreadedTest(gatkName, makeBQSR)

        def makePrintReads(): PrintReads = {
          val PR = new PrintReads with UNIVERSAL_GATK_ARGS
          PR.input_file :+= makeResource(RECAL_BAM_FILENAME)
          PR.out = new File("/dev/null")
          PR.BQSR = makeResource(RECAL_GATKREPORT_FILENAME)
          PR.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName, "assessment" -> "PrintReads"))
          PR.jarFile = gatkJar
          PR
        }
        addMultiThreadedTest(gatkName, makePrintReads)
      }
    }
  }

  def addMultiThreadedTest(gatkName: String,
                           makeCommand: () => CommandLineGATK,
                           maxNT : Int = 1000) {
    if ( ntTests.size >= 1 ) {
      for ( nt <- ntTests ) {
        if ( nt <= maxNT ) {
          for ( useNT <- List(true, false) ) {
            val cmd = makeCommand()
//            Console.print(cmd)
            cmd.nct = nt
            cmd.addJobReportBinding("nct", nt)
            cmd.analysisName = cmd.analysisName + ".nct"
            add(cmd)
          }
        }
      }
    }
  }
}
