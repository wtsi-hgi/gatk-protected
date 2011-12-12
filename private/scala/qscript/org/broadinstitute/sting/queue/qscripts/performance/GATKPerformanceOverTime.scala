package org.broadinstitute.sting.queue.qscripts.performance

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.extensions.gatk._
import java.lang.Math

class GATKPerformanceOverTime extends QScript {
  val STD_RESULTS_DIR = "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/gatkPerformanceOverTime"

  @Argument(shortName = "results", doc="results", required=false)
  val resultsDir: File = new File("runResults");

  @Argument(shortName = "iterations", doc="it", required=false)
  val iterations: Int = 3;

  @Argument(shortName = "steps", doc="steps", required=false)
  val steps: Int = 10;

  val b37: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  class AssessmentParameters(val name: String, val bamList: File, val intervals: File, val nSamples: Int, val dcov: Int)

  // TODO -- count the number of lines in the bam.list file
  val WGSAssessment = new AssessmentParameters("WGS", "wgs.bam.list", "wgs.intervals", 1103, 50)
  val WExAssessment = new AssessmentParameters("WEx", "wex.bam.list", "wex.intervals", 140, 500)
  val assessments = List(WGSAssessment, WExAssessment)

  val GATK12 = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-1.2-65-ge4a583a/GenomeAnalysisTK.jar"
  val GATK13 = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-1.3-23-g13905c0/GenomeAnalysisTK.jar"
  val myJarFile = "/home/unix/depristo/dev/GenomeAnalysisTK/projects/gatkPerformance/dist/GenomeAnalysisTK.jar"
  val GATKs = Map(
    "v1.2" -> GATK12,
    "v1.3" -> GATK13,
    "this" -> myJarFile) // TODO -- how do I get this value?

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = b37;
    this.memoryLimit = 4
  }

  def script = {
    if ( ! resultsDir.exists ) resultsDir.mkdirs()

    for ( assess <- assessments ) {
      for (nSamples <- divideSamples(assess.nSamples) ) {
        val sublist = new SliceList(assess.name, nSamples, assess.bamList)
        add(sublist)
        for ( iteration <- 1 until iterations + 1 ) {
          for ( (gatkName, gatkJar) <- GATKs ) {
            val name: String = "assess.%s_gatk.%s_iter.%d".format(assess.name, gatkName, iteration)

            trait VersionOverrides extends CommandLineGATK {
              this.jarFile = gatkJar
              this.dcov = assess.dcov
              this.intervals :+= assess.intervals
              //this.analysisName = name
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
          }
        }

        //add(new CallWithSamtools(sublist.list, nSamples, name, assess.intervals))
      }
    }
  }

  def divideSamples(nTotalSamples: Int): List[Int] = {
    val maxLog10: Double = Math.log10(nTotalSamples)
    val stepSize: Double = maxLog10 / steps
    val ten: Double = 10.0
    def deLog(x: Int): Int = Math.round(Math.pow(ten, stepSize * x)).toInt
    dedupe(Range(0, steps+1).map(deLog).toList)
  }

  class Call(@Input(doc="foo") bamList: File, n: Int, name: String) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    @Output(doc="foo") var outVCF: File = new File("%s/%s.%d.%s.vcf".format(resultsDir.getPath, bamList.getName, n, name))
    this.input_file :+= bamList
    this.stand_call_conf = 10.0
    this.o = outVCF
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE
  }

  class MyCountLoci(@Input(doc="foo") bamList: File, n: Int, name: String) extends CountLoci with UNIVERSAL_GATK_ARGS {
    @Output(doc="foo") var outFile: File = new File("%s/%s.%d.%s.txt".format(resultsDir.getPath, bamList.getName, n, name))
    this.input_file :+= bamList
    this.o = outFile
  }

  class SliceList(prefix: String, n: Int, bamList: File) extends CommandLineFunction {
    this.analysisName = "SliceList"
    @Output(doc="foo") var list: File = new File("%s/%s.bams.%d.list".format(resultsDir.getPath, prefix, n))
    def commandLine = "head -n %d %s > %s".format(n, bamList, list)
  }

  class CallWithSamtools(@Input(doc="foo") bamList: File, n: Int, name: String, intervals: File) extends CommandLineFunction {
    this.analysisName = "samtools"
    @Output(doc="foo") var outVCF: File = new File("%s/%s.%d.%s.bcf".format(resultsDir.getPath, bamList.getName, n, name))
    def commandLine = "./samtools mpileup -ub %s -f %s -l %s | ./bcftools view -vcg -l %s - > %s".format(bamList, b37, intervals, intervals, outVCF)
  }

  def dedupe(elements:List[Int]):List[Int] = {
    if (elements.isEmpty)
      elements
    else
      elements.head :: dedupe(for (x <- elements.tail if x != elements.head) yield x)
  }
}
