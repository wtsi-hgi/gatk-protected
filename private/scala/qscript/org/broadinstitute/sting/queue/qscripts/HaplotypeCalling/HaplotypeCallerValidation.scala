package org.broadinstitute.sting.queue.qscripts.HaplotypeCalling

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.text.XReadLines
import collection.JavaConversions._
import org.broadinstitute.sting.utils.{GenomeLoc, GenomeLocParser}
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource
import org.broadinstitute.sting.queue.function.QFunction
import scala.math._

class HaplotypeCallerValidation extends QScript {
  qscript =>

  @Input(doc = "File listing run name, locus, and then samples to produce haplotype-based genotype calls", shortName = "I", required = true)
  var runsFile: File = _

  @Input(doc = "File mapping sample to BAM and SM tag in that BAM", shortName = "sample_bam_SM", required = true)
  var sample_bam_SM: File = _

  @Argument(doc = "gatk jar file", shortName = "J", required = true)
  var gatkJarFile: File = _

  @Argument(shortName = "R", doc = "ref", required = true)
  var referenceFile: File = _

  @Input(doc = "level of parallelism.   By default is set to 0 [no scattering].", shortName = "scatter", required = false)
  var scatterCountInput = 0

  @Input(doc = "Bases upstream and downstream to add when a single base locus is given", shortName = "extent", required = false)
  var defaultExtent = 100


  class BamSM(bamIn: File, SMin: String) {
    val bam = bamIn
    val SM = SMin
  }

  def script = {
    val sampleToBamSM = sampleToBAM_SMfromMapFile(sample_bam_SM)

    class HCrun(val name: String, val locus: GenomeLoc, val samples: List[String]) extends HaplotypeCaller with CommandLineGATKArgs {
      this.intervalsString = List(locus.toString)
      this.input_file = samples.reverse.map(s => {if (sampleToBamSM.contains(s)) sampleToBamSM(s).bam else throw new IllegalArgumentException("Sample " + s + " not found in " + sample_bam_SM)})
      this.out = name + ".HC.vcf"

      // Get full haplotypes (at the expense of missing out on the lower frequency variants with lots of samples):
      this.genotypeFullActiveRegion = true
    }

    class UGrun(val name: String, val locus: GenomeLoc, val samples: List[String]) extends UnifiedGenotyper with CommandLineGATKArgs {
      this.intervalsString = List(locus.toString)
      this.input_file = samples.reverse.map(s => {if (sampleToBamSM.contains(s)) sampleToBamSM(s).bam else throw new IllegalArgumentException("Sample " + s + " not found in " + sample_bam_SM)})
      this.out = name + ".UG.vcf"

      this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
      this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    }

    def createRuns(runsFile: File) : List[QFunction] = {
      var locParser = new GenomeLocParser(new ReferenceDataSource(referenceFile).getReference)
      var runs = List[QFunction]()

      var elems = asScalaIterator(new XReadLines(runsFile))
      while (elems.hasNext) {
        val line = elems.next
        val splitLine = line.split("\\s+")
        val name = splitLine(0)
        val locusStr = splitLine(1)
        var samples = List[String]()

        for (i <- 2 until splitLine.length) {
          samples ::= splitLine(i)
        }

        val splitLoc = locusStr.split("\\+")
        var locus = locParser.parseGenomeLoc(splitLoc(0))

        if (locus.getStart == locus.getStop) {
          var extent = defaultExtent
          if (splitLoc.length > 1) {
            extent = splitLoc(1).toInt
          }
          locus = locParser.setStart(locus, max(1, locus.getStart - extent))
          locus = locParser.setStop(locus, min(locParser.getContigInfo(locus.getContig).getSequenceLength, locus.getStop + extent))
        }

        runs ::= new HCrun(name, locus, samples)
        runs ::= new UGrun(name, locus, samples)
      }

      return runs
    }

    val runs = createRuns(runsFile)
    addAll(runs)
  }

  trait CommandLineGATKArgs extends CommandLineGATK with ScatterGatherableFunction {
    this.jarFile = qscript.gatkJarFile
    this.reference_sequence = qscript.referenceFile
    this.scatterCount = scatterCountInput
    this.memoryLimit = 2
    this.logging_level = "INFO"
  }

  def sampleToBAM_SMfromMapFile(sample_bam_SM_file: File) : scala.collection.mutable.Map[String, BamSM] = {
    var sampMap = scala.collection.mutable.Map.empty[String, BamSM]

    var elems = asScalaIterator(new XReadLines(sample_bam_SM_file))
    while (elems.hasNext) {
      val line = elems.next
      val splitLine = line.split("\\s+")
      val sample = splitLine(0)
      val bam = new File(splitLine(1))
      val SM = splitLine(2)

      sampMap.put(sample, new BamSM(bam, SM))
    }

    return sampMap
  }

}