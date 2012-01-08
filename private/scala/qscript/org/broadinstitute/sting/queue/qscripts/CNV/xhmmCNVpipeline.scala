package org.broadinstitute.sting.queue.qscripts.CNV

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.queue.util.VCF_BAM_utilities
import org.broadinstitute.sting.commandline.Hidden

class xhmmCNVpipeline extends QScript {
  qscript =>

  @Input(doc = "bam input, as .bam or as a list of files", shortName = "I", required = true)
  var bams: File = _

  @Argument(doc = "gatk jar file", shortName = "J", required = true)
  var gatkJarFile: File = _

  @Argument(doc = "xhmm executable file", shortName = "X", required = true)
  var xhmmExec: File = _

  @Argument(shortName = "R", doc = "ref", required = true)
  var referenceFile: File = _

  @Argument(shortName = "L", doc = "Intervals", required = false)
  var intervals: String = _

  @Input(doc = "level of parallelism for BAM DoC.   By default is set to 0 [no scattering].", shortName = "scatter", required = false)
  var scatterCountInput = 0

  @Input(doc = "Samples to run together for DoC.   By default is set to 1 [one job per sample].", shortName = "samplesPerJob", required = false)
  var samplesPerJob = 1

  @Output(doc = "Base name for files to output", shortName = "o", required = true)
  var outputBase: File = _

  @Input(doc = "Maximum depth (before GATK down-sampling kicks in...)", shortName = "MAX_DEPTH", required = false)
  var MAX_DEPTH = 20000

  @Hidden
  @Input(doc = "Number of read-depth bins", shortName = "NUM_BINS", required = false)
  var NUM_BINS = 200

  @Hidden
  @Input(doc = "Starting value of read-depth bins", shortName = "START_BIN", required = false)
  var START_BIN = 1

  @Input(doc = "Minimum read mapping quality", shortName = "MMQ", required = false)
  var minMappingQuality = 0

  @Input(doc = "Memory (in GB) required for storing the whole matrix in memory", shortName = "wholeMatrixMemory", required = false)
  var wholeMatrixMemory = -1

  @Argument(shortName = "rawFilters", doc = "xhmm command-line parameters to filter targets and samples from raw data", required = false)
  var targetSampleFiltersString: String = ""

  @Argument(shortName = "PCAnormalize", doc = "xhmm command-line parameters to Normalize data using PCA information", required = false)
  var PCAnormalizeMethodString: String = ""

  @Argument(shortName = "normalizedFilters", doc = "xhmm command-line parameters to filter targets and samples from PCA-normalized data", required = false)
  var targetSampleNormalizedFiltersString: String = ""

  @Argument(shortName = "xhmmParams", doc = "xhmm model parameters file", required = true)
  var xhmmParamsArg: File = _

  val DOC_OUTPUT_SUFFIX: String = ".sample_interval_summary"

  val RD_OUTPUT_SUFFIX: String = ".RD.txt"
  val FILTERED_TARGS_SUFFIX: String = ".filtered_targets.txt"
  val FILTERED_SAMPS_SUFFIX: String =  ".filtered_samples.txt"

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervalsString = List(qscript.intervals)
    this.jarFile = qscript.gatkJarFile
    this.reference_sequence = qscript.referenceFile
    this.logging_level = "INFO"
  }

  trait WholeMatrixMemoryLimit extends CommandLineFunction {
    // Since loading ALL of the data can take significant memory:
    if (wholeMatrixMemory < 0) {
      this.memoryLimit = 24
    }
    else {
      this.memoryLimit = wholeMatrixMemory
    }
  }

  // A target has a list of samples and bam files to use for DoC
  class Target(val name: String, val samples: List[String], val bams: List[File]) {
    // getName() just includes the file name WITHOUT the path:
    def DoC_output = new File(name + "." + outputBase.getName())

    override def toString(): String = String.format("[Target %s [%s] with samples %s against bams %s]", name, DoC_output, samples, bams)
  }

  def script = {
    val sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]] = VCF_BAM_utilities.getMapOfBAMsForSample(VCF_BAM_utilities.parseBAMsInput(bams))
    val samples: List[String] = sampleToBams.keys.toList
    Console.out.printf("Samples are %s%n", samples)

    val targets: List[Target] = buildTargets(samples, sampleToBams)
    for (target <- targets) {
      Console.out.printf("Target is %s%n", target)
      add(new DoC(target))
    }

    val mergeDepths = new MergeGATKdepths(targets.map(u => new File(u.DoC_output.getPath() + DOC_OUTPUT_SUFFIX)))
    add(mergeDepths)

    val filterCenterDepths = new FilterCenterRawMatrix(mergeDepths.mergedDoC)
    add(filterCenterDepths)

    val pca = new PCA(filterCenterDepths.filteredCentered)
    add(pca)

    val normalize = new Normalize(pca)
    add(normalize)

    val filterZscore = new FilterAndZscoreNormalized(normalize.normalized)
    add(filterZscore)

    val filterOriginal = new FilterOriginalData(mergeDepths.mergedDoC, filterCenterDepths, filterZscore)
    add(filterOriginal)

    val discover = new DiscoverCNVs(filterZscore.filteredZscored, filterOriginal.sameFiltered)
    add(discover)

    val genotype = new GenotypeCNVs(filterZscore.filteredZscored, discover.xcnv, filterOriginal.sameFiltered)
    add(genotype)
  }

  def buildTargets(samples: List[String], sampleToBams: scala.collection.mutable.Map[String, scala.collection.mutable.Set[File]]): List[Target] = {

    def buildTargetsHelper(samples: List[String], count: Int): List[Target] = (samples splitAt samplesPerJob) match {
      case (Nil, y) =>
        return Nil
      case (subsamples, remaining) =>
        return new Target("group" + count, subsamples, VCF_BAM_utilities.findBAMsForSamples(subsamples, sampleToBams)) ::
          buildTargetsHelper(remaining, count + 1)
    }

    return buildTargetsHelper(samples, 0)
  }

  class DoC(t: Target) extends CommandLineGATKArgs with ScatterGatherableFunction {
    this.analysis_type = "DepthOfCoverage"

    this.input_file = t.bams

    this.downsample_to_coverage = MAX_DEPTH
    this.downsampling_type = DownsampleType.BY_SAMPLE

    this.scatterCount = scatterCountInput
    this.scatterClass = classOf[IntervalScatterFunction]

    @Output
    @Gather(classOf[org.broadinstitute.sting.queue.function.scattergather.SimpleTextGatherFunction])
    var intervalSampleOut: File = new File(t.DoC_output.getPath() + DOC_OUTPUT_SUFFIX)

    // The HACK for DoC to work properly within Queue:
    val commandLineSuppliedOutputFilesPrefix = new File(intervalSampleOut.getParentFile(), t.DoC_output.getName())

    override def commandLine = super.commandLine +
      " --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 0 --minMappingQuality " + minMappingQuality +
      " --start " + START_BIN + " --stop " + MAX_DEPTH + " --nBins " + NUM_BINS +
      " -o " + commandLineSuppliedOutputFilesPrefix

    override def dotString = "DOC: " + t.DoC_output

    this.jobOutputFile = commandLineSuppliedOutputFilesPrefix + ".out"
  }

  class MergeGATKdepths(DoCsToCombine: List[File]) extends CommandLineFunction with WholeMatrixMemoryLimit {
    @Input(doc = "")
    var inputDoCfiles: List[File] = DoCsToCombine

    @Output
    val mergedDoC: File = new File(outputBase.getPath + RD_OUTPUT_SUFFIX)

    var command: String =
      xhmmExec + " --mergeGATKdepths" +
      " -o " + mergedDoC
    for (input <- inputDoCfiles) {
      command += " --GATKdepths " + input
    }
    def commandLine = command

    override def description = "Combines DoC outputs for multiple samples (at same loci): " + command
  }

  class FilterCenterRawMatrix(inputParam: File) extends CommandLineFunction with WholeMatrixMemoryLimit {
    @Input(doc = "")
    val input = inputParam

    @Output
    val filteredCentered: File = new File(outputBase.getPath + ".filtered_centered" + RD_OUTPUT_SUFFIX)
    @Output
    val filteredTargets: File = new File(filteredCentered.getPath + FILTERED_TARGS_SUFFIX)
    @Output
    val filteredSamples: File = new File(filteredCentered.getPath + FILTERED_SAMPS_SUFFIX)

    var command: String =
      xhmmExec + " --matrix" +
      " -r " + input +
      " --centerData --centerType target" +
      " -o " + filteredCentered +
      " --outputExcludedTargets " + filteredTargets +
      " --outputExcludedSamples " + filteredSamples
    if (targetSampleFiltersString != "")
      command += " " + targetSampleFiltersString

    def commandLine = command

    override def description = "Filters samples and targets and then mean-centers the targets: " + command
  }

  class PCA(inputParam: File) extends CommandLineFunction with WholeMatrixMemoryLimit {
    @Input(doc = "")
    val input = inputParam

    val PCAbase: String = outputBase.getPath + ".RD_PCA"

    @Output
    val outPC: File = new File(PCAbase + ".PC.txt")
    @Output
    val outPC_SD: File = new File(PCAbase + ".PC_SD.txt")
    @Output
    val outPC_LOADINGS: File = new File(PCAbase + ".PC_LOADINGS.txt")

    var command: String =
      xhmmExec + " --PCA" +
      " -r " + input +
      " --PCAfiles " + PCAbase

    def commandLine = command

    override def description = "Runs PCA on mean-centered data: " + command
  }

  class Normalize(pca: PCA) extends CommandLineFunction {
    @Input(doc = "")
    val input = pca.input

    @Input(doc = "")
    val inPC = pca.outPC

    @Input(doc = "")
    val inPC_SD = pca.outPC_SD

    @Input(doc = "")
    val inPC_LOADINGS = pca.outPC_LOADINGS

    @Output
    val normalized: File = new File(outputBase.getPath + ".PCA_normalized.txt")

    var command: String =
      xhmmExec + " --normalize" +
      " -r " + input +
      " --PCAfiles " + pca.PCAbase +
      " --normalizeOutput " + normalized
    if (PCAnormalizeMethodString != "")
      command += " " + PCAnormalizeMethodString

    def commandLine = command

    override def description = "Normalizes mean-centered data using PCA information: " + command
  }

  class FilterAndZscoreNormalized(inputParam: File) extends CommandLineFunction with WholeMatrixMemoryLimit {
    @Input(doc = "")
    val input = inputParam

    @Output
    val filteredZscored: File = new File(outputBase.getPath + ".PCA_normalized.filtered.sample_zscores" + RD_OUTPUT_SUFFIX)
    @Output
    val filteredTargets: File = new File(filteredZscored.getPath + FILTERED_TARGS_SUFFIX)
    @Output
    val filteredSamples: File = new File(filteredZscored.getPath + FILTERED_SAMPS_SUFFIX)

    var command: String =
      xhmmExec + " --matrix" +
      " -r " + input +
      " --centerData --centerType sample --zScoreData" +
      " -o " + filteredZscored +
      " --outputExcludedTargets " + filteredTargets +
      " --outputExcludedSamples " + filteredSamples
    if (targetSampleNormalizedFiltersString != "")
      command += " " + targetSampleNormalizedFiltersString

    def commandLine = command

    override def description = "Filters and z-score centers (by sample) the PCA-normalized data: " + command
  }

  class FilterOriginalData(inputParam: File, filt1: FilterCenterRawMatrix, filt2: FilterAndZscoreNormalized) extends CommandLineFunction with WholeMatrixMemoryLimit {
    @Input(doc = "")
    val input = inputParam

    @Input(doc = "")
    val targFilters: List[File] = List(filt1.filteredTargets, filt2.filteredTargets).map(u => new File(u))

    @Input(doc = "")
    val sampFilters: List[File] = List(filt1.filteredSamples, filt2.filteredSamples).map(u => new File(u))

    @Output
    val sameFiltered: File = new File(outputBase.getPath + ".same_filtered" + RD_OUTPUT_SUFFIX)

    var command: String =
      xhmmExec + " --matrix" +
      " -r " + input +
      targFilters.map(u => " --excludeTargets " + u).reduceLeft(_ + "" + _) +
      sampFilters.map(u => " --excludeSamples " + u).reduceLeft(_ + "" + _) +
      " -o " + sameFiltered

    def commandLine = command

    override def description = "Filters original read-depth data to be the same as filtered, normalized data: " + command
  }

  class DiscoverCNVs(inputParam: File, origRDParam: File) extends CommandLineFunction {
    @Input(doc = "")
    val input = inputParam

    @Input(doc = "")
    val xhmmParams = xhmmParamsArg

    @Input(doc = "")
    val origRD = origRDParam

    @Output
    val xcnv: File = new File(outputBase.getPath + ".xcnv")

    @Output
    val aux_xcnv: File = new File(outputBase.getPath + ".aux_xcnv")

    val posteriorsBase = outputBase.getPath

    @Output
    val dipPosteriors: File = new File(posteriorsBase + ".posteriors.DIP.txt")

    @Output
    val delPosteriors: File = new File(posteriorsBase + ".posteriors.DEL.txt")

    @Output
    val dupPosteriors: File = new File(posteriorsBase + ".posteriors.DUP.txt")

    var command: String =
      xhmmExec + " --discover" +
      " -p " + xhmmParams +
      " -r " + input +
      " -R " + origRD +
      " -c " + xcnv +
      " -a " + aux_xcnv +
      " -s " + posteriorsBase

    def commandLine = command

    override def description = "Discovers CNVs in normalized data: " + command
  }

  class GenotypeCNVs(inputParam: File, xcnv: File, origRDParam: File) extends CommandLineFunction {
    @Input(doc = "")
    val input = inputParam

    @Input(doc = "")
    val xhmmParams = xhmmParamsArg

    @Input(doc = "")
    val origRD = origRDParam

    @Input(doc = "")
    val inXcnv = xcnv

    @Output
    val vcf: File = new File(outputBase.getPath + ".vcf")

    var command: String =
      xhmmExec + " --genotype" +
      " -p " + xhmmParams +
      " -r " + input +
      " -g " + inXcnv +
      " -F " + referenceFile +
      " -R " + origRD +
      " -v " +  vcf

    def commandLine = command

    override def description = "Genotypes discovered CNVs in all samples: " + command
  }
}