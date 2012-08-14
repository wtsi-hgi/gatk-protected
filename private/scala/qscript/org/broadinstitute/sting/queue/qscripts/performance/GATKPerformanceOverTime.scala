package org.broadinstitute.sting.queue.qscripts.performance

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import java.lang.Math
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode
import org.broadinstitute.sting.utils.PathUtils

class GATKPerformanceOverTime extends QScript {
  val STD_RESULTS_DIR = "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/gatkPerformanceOverTime"

  @Argument(shortName = "results", doc="results", required=false)
  val resultsDir: File = new File("runResults")

  @Argument(shortName = "test", doc="test", required=false)
  val TEST: Boolean = false

  @Argument(shortName = "justUG", doc="Just run UnifiedGenotyper tests", required=false)
  val justUG: Boolean = false

  @Argument(shortName = "resources", doc="resources", required=true)
  val resourcesDir: String = ""

  @Argument(shortName = "myJarFile", doc="Path to the current GATK jar file", required=true)
  val myJarFile: File = null

  @Argument(shortName = "iterations", doc="it", required=false)
  val iterations: Int = 3

  val nIterationsForSingleTestsPerIteration: Int = 3

  @Argument(shortName = "ntTest", doc="For each value provided we will use -nt VALUE in the multi-threaded tests", required=false)
  val ntTests: List[Int] = List(1, 2, 3, 4, 6, 8, 10, 12, 16, 24, 32)

  @Argument(shortName = "steps", doc="steps", required=false)
  val steps: Int = 10

  @Argument(shortName = "maxNSamples", doc="maxNSamples", required=false)
  val maxNSamples: Int = 1000000

  val RECAL_BAM_FILENAME = "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam"
  val dbSNP_FILENAME = "dbsnp_132.b37.vcf"
  val RECAL_FILENAME = "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.csv"  // TODO -- update to use recal table for BQSRv2
  val b37_FILENAME = "human_g1k_v37.fasta"

  def makeResource(x: String): File = new File("%s/%s".format(resourcesDir, x))
  def makeChunk(x: Int): File = makeResource("chunk_%d.vcf".format(x))
  def COMBINE_FILES: List[File] = Range(1,10).map(makeChunk).toList

  class AssessmentParameters(val name: String, val bamList: File, val intervals: File, val nSamples: Int, val dcov: Int)

  // TODO -- count the number of lines in the bam.list file
  val WGSAssessment = new AssessmentParameters("WGS", "wgs.bam.list.local.list", "wgs.bam.list.intervals", 1103, 50)
  val WGSDeepAssessment = new AssessmentParameters("WGS.deep", "wgs.deep.bam.list.local.list", "wgs.deep.bam.list.intervals", 1, 250)
  val WExAssessment = new AssessmentParameters("WEx", "wex.bam.list.local.list", "wex.bam.list.intervals", 140, 500)

  val assessments = List(WGSAssessment, WGSDeepAssessment, WExAssessment)

  val GATK_RELEASE_DIR = new File("/humgen/gsa-hpprojects/GATK/bin/")
  val GATKs: Map[String, File] = Map(
    "v2.cur" -> myJarFile, // TODO -- how do I get this value?
    "v1.6" -> findMostRecentGATKVersion("1.6"))

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO"
    this.reference_sequence = makeResource(b37_FILENAME)
    this.memoryLimit = 8
  }

  def script() {
    if ( ! resultsDir.exists ) resultsDir.mkdirs()

    for ( iteration <- 0 until iterations ) {
      for ( assess <- assessments ) {
        for (nSamples <- divideSamples(assess.nSamples) ) {
          val sublist = new SliceList(assess.name, nSamples, makeResource(assess.bamList))
          if ( iteration == 0 ) add(sublist) // todo - remove condition when Queue bug is fixed
          for ( (gatkName, gatkJar) <- GATKs ) {
            val name: String = "assess.%s_gatk.%s_iter.%d".format(assess.name, gatkName, iteration)

            trait VersionOverrides extends CommandLineGATK {
              this.jarFile = gatkJar
              this.dcov = assess.dcov

              // special handling of test intervals
              if ( TEST )
                this.intervalsString :+= "20:10,000,000-10,001,000"
              else
                this.intervals :+= makeResource(assess.intervals)

              this.configureJobReport(Map(
                "iteration" -> iteration,
                "gatk" -> gatkName,
                "nSamples" -> nSamples,
                "assessment" -> assess.name))
            }

            // SNP calling
            add(new Call(sublist.list, nSamples, name) with VersionOverrides)

            // CountLoci
            add(new MyCountLoci(sublist.list, nSamples, name) with VersionOverrides)

            if ( nSamples == assess.nSamples )
              addMultiThreadedTest(() => new Call(sublist.list, nSamples, name) with VersionOverrides)
          }
        }
      }

      if ( ! justUG ) {
        for ( subiteration <- 0 until nIterationsForSingleTestsPerIteration ) {
          for ( (gatkName, gatkJar) <- GATKs ) {
            { // Standard VCF tools
            trait VersionOverrides extends CommandLineGATK {
              this.jarFile = gatkJar
              this.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName))
            }

              val CV = new CombineVariants with UNIVERSAL_GATK_ARGS with VersionOverrides
              CV.variant = COMBINE_FILES
              CV.intervalsString = (if ( TEST ) List("1:10,000,000-10,010,000") else List("1", "2", "3", "4", "5"))
              CV.out = new File("/dev/null")
              add(CV)

              val SV = new SelectVariants with UNIVERSAL_GATK_ARGS with VersionOverrides
              SV.variant = makeResource("chunk_1.vcf")
              SV.sample_name = List("HG00096") // IMPORTANT THAT THIS SAMPLE BE IN CHUNK ONE
              if ( TEST ) SV.intervalsString = List("1:10,000,000-10,010,000")
              SV.out = new File("/dev/null")
              add(SV)

              def makeVE(): CommandLineGATK = {
                val VE = new VariantEval with UNIVERSAL_GATK_ARGS with VersionOverrides
                VE.eval :+= makeResource("chunk_1.vcf")
                if ( TEST ) VE.intervalsString = List("1:10,000,000-10,010,000")
                VE.out = new File("/dev/null")
                VE.comp :+= new TaggedFile(makeResource(dbSNP_FILENAME), "dbSNP")
                VE
              }

              add(makeVE())
              if ( subiteration == 0 )
                addMultiThreadedTest(makeVE)
              //            }
            }
          }
        }
      }
    }
  }

  def addMultiThreadedTest(makeCommand: () => CommandLineGATK) {
    if ( ntTests.size > 1 ) {
      for ( nt <- ntTests ) {
        val cmd = makeCommand()
        cmd.nt = nt
        cmd.addJobReportBinding("nt", nt)
        cmd.analysisName = cmd.analysisName + ".nt"
        add(cmd)
      }
    }
  }

  def divideSamples(nTotalSamples: Int): List[Int] = {
    val maxLog10: Double = Math.log10(Math.min(maxNSamples, nTotalSamples))
    val stepSize: Double = maxLog10 / steps
    val ten: Double = 10.0
    def deLog(x: Int): Int = Math.round(Math.pow(ten, stepSize * x)).toInt
    dedupe(Range(0, steps+1).map(deLog).toList)
  }

  class Call(@Input(doc="foo") bamList: File, n: Int, name: String) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.input_file :+= bamList
    this.stand_call_conf = 10.0
    this.o = outVCF
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE
    @Output(doc="foo") var outVCF: File = new File("/dev/null")
  }

  class MyCountLoci(@Input(doc="foo") bamList: File, n: Int, name: String) extends CountLoci with UNIVERSAL_GATK_ARGS {
    this.input_file :+= bamList
    @Output(doc="foo") var outFile: File = new File("/dev/null")
    this.o = outFile
  }

  class SliceList(prefix: String, n: Int, @Input bamList: File) extends CommandLineFunction {
    this.analysisName = "SliceList"
    @Output(doc="foo") var list: File = new File("%s/%s.bams.%d.list".format(resultsDir.getPath, prefix, n))
    def commandLine = "head -n %d %s | awk '{print \"%s/\" $1}' > %s".format(n, bamList, resourcesDir, list)
  }

  def dedupe(elements:List[Int]):List[Int] = {
    if (elements.isEmpty)
      elements
    else
      elements.head :: dedupe(for (x <- elements.tail if x != elements.head) yield x)
  }

  /**
   * Walk over the GATK released directories to find the most recent JAR files corresponding
   * to the version prefix.  For example, providing input "GenomeAnalysisTK-1.2" will
   * return the full path to the most recent GenomeAnalysisTK.jar in the GATK_RELEASE_DIR
   * in directories that match GATK_RELEASE_DIR/GenomeAnalysisTK-1.2*
   */
  def findMostRecentGATKVersion(version: String): File = {
    PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, version)
  }
}
