package org.broadinstitute.sting.queue.qscripts.performance

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import java.lang.Math
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode
import org.broadinstitute.sting.utils.PathUtils
import org.broadinstitute.sting.queue.function.QFunction
import org.broadinstitute.sting.queue.engine.JobRunInfo

class GATKPerformanceOverTime extends QScript {
  @Argument(shortName = "results", doc="results", required=false)
  val resultsDir: File = new File("runResults")

  @Argument(shortName = "resources", doc="resources", required=true)
  val resourcesDir: String = ""

  @Argument(shortName = "myJarFile", doc="Path to the current GATK jar file", required=true)
  val myJarFile: File = null

  @Argument(shortName = "iterations", doc="it", required=false)
  val iterations: Int = 2

  @Argument(shortName = "assessment", doc="Which assessments should we run?", required=false)
  val assessmentsArg: Set[String] = Assessment.values map(_.toString)

  val nIterationsForSingleTestsPerIteration: Int = 3

  @Argument(shortName = "ntTest", doc="For each value provided we will use -nt VALUE in the multi-threaded tests", required=false)
  val ntTests: List[Int] = List(1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 24)

  @Argument(shortName = "steps", doc="steps", required=false)
  val steps: Int = 10

  @Argument(shortName = "maxNSamples", doc="maxNSamples", required=false)
  val maxNSamples: Int = 1000000

  val MY_TAG = "GATKPerformanceOverTime"
  val RECAL_BAM_FILENAME = "wgs.deep.bam.list.local.list"
  val dbSNP_FILENAME = "dbsnp_132.b37.vcf"
  val BIG_VCF_WITH_GENOTYPES = "ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
  val RECAL_FILENAME = "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.csv"  // TODO -- update to use recal table for BQSRv2
  val b37_FILENAME = "human_g1k_v37.fasta"

  def makeResource(x: String): File = new File("%s/%s".format(resourcesDir, x))
  def makeChunk(x: Int): File = makeResource("chunk_%d.vcf".format(x))
  def COMBINE_FILES: List[File] = Range(1,10).map(makeChunk).toList

  class AssessmentParameters(val name: String,
                             val bamList: File,
                             val fullIntervals: String,
                             val shortIntervals: String,
                             val nSamples: Int,
                             val dcov: Int,
                             val baq: Boolean) {
    def addIntervals(gatkCmd : CommandLineGATK, useFull: Boolean): CommandLineGATK = {
      val intervals = if (useFull) fullIntervals else shortIntervals
      val maybeFile = makeResource(intervals)
      if ( maybeFile.exists() )
        gatkCmd.intervals :+= maybeFile
      else
        gatkCmd.intervalsString :+= intervals
      gatkCmd
    }
  }

  // TODO -- count the number of lines in the bam.list file
  val WGSAssessment = new AssessmentParameters("WGS.multiSample.4x", "wgs.bam.list.local.list", "wgs.bam.list.select.intervals", "20:10000000-11000000", 1103, 50, true)
  val WGSDeepAssessment = new AssessmentParameters("WGS.singleSample.60x", "wgs.deep.bam.list.local.list", "wgs.deep.bam.list.select.intervals", "1", 1, 250, true)
  val WExAssessment = new AssessmentParameters("WEx.multiSample.150x", "wex.bam.list.local.list", "wex.bam.list.select.intervals", "wex.bam.list.small.intervals", 140, 500, true)

  val dataSets = List(WGSAssessment, WGSDeepAssessment, WExAssessment)

  val GATK_RELEASE_DIR = new File("/humgen/gsa-hpprojects/GATK/bin/")
  val GATKs: Map[String, File] = Map(
    "v2.cur" -> myJarFile, // TODO -- how do I get this value?
    "v1.6" -> findMostRecentGATKVersion("1.6"))

  object Assessment extends Enumeration {
    type Assessment = Value
    val UG, UG_NT, CL, CL_NT, CV, CV_NT, VE, VE_NT, SV, BQSR_NT = Value
  }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO"
    this.reference_sequence = makeResource(b37_FILENAME)
    this.memoryLimit = 4
  }

  def script() {
    val assessments = assessmentsArg.map(Assessment.withName(_))

    if ( ! resultsDir.exists ) resultsDir.mkdirs()

    for ( iteration <- 0 until iterations ) {
      for ( assess <- dataSets ) {
        for (nSamples <- divideSamples(assess.nSamples) ) {
          val sublist = new SliceList(assess.name, nSamples, makeResource(assess.bamList))
          if ( iteration == 0 ) add(sublist) // todo - remove condition when Queue bug is fixed
          for ( (gatkName, gatkJar) <- GATKs ) {
            val name: String = "assess.%s_gatk.%s_iter.%d".format(assess.name, gatkName, iteration)

            trait VersionOverrides extends CommandLineGATK {
              this.jarFile = gatkJar
              this.dcov = assess.dcov

              this.configureJobReport(Map(
                "iteration" -> iteration,
                "gatk" -> gatkName,
                "nSamples" -> nSamples,
                "assessment" -> assess.name))
            }

            // SNP calling
            if ( assessments.contains(Assessment.UG) )
              addGATKCommand(assess.addIntervals(new Call(sublist.list, nSamples, name, assess.baq) with VersionOverrides, false))
            if ( assessments.contains(Assessment.UG_NT) && nSamples == assess.nSamples )
              addMultiThreadedTest(() => assess.addIntervals(new Call(sublist.list, nSamples, name, assess.baq) with VersionOverrides, false))

            // CountLoci
            if ( assessments.contains(Assessment.CL) )
              addGATKCommand(assess.addIntervals(new MyCountLoci(sublist.list, nSamples, name) with VersionOverrides, true))
            if ( assessments.contains(Assessment.CL_NT) && nSamples == assess.nSamples )
              addMultiThreadedTest(() => assess.addIntervals(new MyCountLoci(sublist.list, nSamples, name) with VersionOverrides, true))
          }
        }
      }

      // GATK v2 specific tests
      for ( iteration <- 0 until iterations ) {
        for ( (gatkName, gatkJar) <- GATKs ) {
          if ( gatkName.contains("v2") ) {
            if ( assessments.contains(Assessment.CV_NT) ) {
              for ( outputBCF <- List(true, false) ) {
                val outputName = if ( outputBCF ) "bcf" else "vcf"

                def makeCV(): CommandLineGATK = {
                  val CV = new CombineVariants with UNIVERSAL_GATK_ARGS
                  CV.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName, "assessment" -> outputName))
                  CV.jarFile = gatkJar
                  CV.intervalsString :+= "22"
                  CV.variant = List(makeResource(BIG_VCF_WITH_GENOTYPES))
                  CV.out = new File("/dev/null")
                  CV.bcf = outputBCF
                  CV
                }

                addMultiThreadedTest(makeCV)
              }
            }

            if ( assessments.contains(Assessment.BQSR_NT) ) {
              def makeBQSR(): CommandLineGATK = {
                val BQSR = new BaseRecalibrator with UNIVERSAL_GATK_ARGS
                BQSR.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName, "assessment" -> RECAL_BAM_FILENAME))
                BQSR.jarFile = gatkJar
                BQSR.intervalsString :+= "1:10,100,000-100,000,000"
                BQSR.knownSites :+= makeResource(dbSNP_FILENAME)
                //this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "ContextCovariate")
                BQSR.input_file :+= makeResource(RECAL_BAM_FILENAME)
                BQSR.out = new File("/dev/null")
                BQSR.no_plots = true
                BQSR.memoryLimit = 12
                BQSR
              }

              addMultiThreadedTest(makeBQSR, 8) // max nt until BQSR is performant
            }
          }
        }
      }

      for ( subiteration <- 0 until nIterationsForSingleTestsPerIteration ) {
        for ( (gatkName, gatkJar) <- GATKs ) {
          { // Standard VCF tools
          trait VersionOverrides extends CommandLineGATK {
            this.jarFile = gatkJar
            this.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName))
          }

            val CV = new CombineVariants with UNIVERSAL_GATK_ARGS with VersionOverrides
            CV.variant = COMBINE_FILES
            CV.intervalsString = List("1", "2", "3", "4", "5")
            CV.out = new File("/dev/null")
            if ( assessments.contains(Assessment.CV) )
              addGATKCommand(CV)

            val SV = new SelectVariants with UNIVERSAL_GATK_ARGS with VersionOverrides
            SV.variant = makeResource("chunk_1.vcf")
            SV.sample_name = List("HG00096") // IMPORTANT THAT THIS SAMPLE BE IN CHUNK ONE
            SV.out = new File("/dev/null")
            if ( assessments.contains(Assessment.SV) )
              addGATKCommand(SV)

            def makeVE(): CommandLineGATK = {
              val VE = new VariantEval with UNIVERSAL_GATK_ARGS with VersionOverrides
              VE.eval :+= makeResource("chunk_1.vcf")
              VE.out = new File("/dev/null")
              VE.comp :+= new TaggedFile(makeResource(dbSNP_FILENAME), "dbSNP")
              VE.addJobReportBinding("assessment", "chunk_1.vcf")
              VE
            }

            if ( assessments.contains(Assessment.VE) ) {
              addGATKCommand(makeVE())
            }

            if ( assessments.contains(Assessment.VE_NT) && subiteration == 0 )
              addMultiThreadedTest(makeVE)
          }
        }
      }
    }
  }

  def addMultiThreadedTest(makeCommand: () => CommandLineGATK, maxNT : Int = 1000) {
    if ( ntTests.size > 1 ) {
      for ( nt <- ntTests ) {
        if ( nt < maxNT ) {
          val cmd = makeCommand()
          cmd.nt = nt
          cmd.memoryLimit = cmd.memoryLimit * (if ( nt >= 8 ) (if (nt>=16) 4 else 2) else 1)
          cmd.addJobReportBinding("nt", nt)
          cmd.analysisName = cmd.analysisName + ".nt"
          addGATKCommand(cmd)
        }
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

  class Call(@Input(doc="foo") bamList: File, n: Int, name: String, useBaq: Boolean) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.input_file :+= bamList
    this.stand_call_conf = 10.0
    this.o = outVCF
    this.baq = if ( useBaq ) org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE else org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
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

  def addGATKCommand(gatkCmd: CommandLineGATK) {
    if ( gatkCmd.jarFile == null || ! gatkCmd.jarFile.getAbsolutePath.matches(".*-1.[0-9]*-.*") )
      gatkCmd.tag = MY_TAG
    add(gatkCmd)
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
