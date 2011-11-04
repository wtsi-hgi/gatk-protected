package org.broadinstitute.sting.queue.qscripts.annotations

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.library.ipf.vcf.VCFExtractIntervals

class ScatteredFullVariantAnnotator extends QScript {
  qscript =>

  @Argument(doc = "gatk jar file", shortName = "J", required = true)
  var gatkJarFile: File = _

  @Argument(shortName = "R", doc = "ref", required = true)
  var referenceFile: File = _

  @Input(doc = "level of parallelism. By default is set to 0 [no scattering].", shortName = "scatter", required = false)
  var scatterCount = 0

  @Input(doc = "bam input, as .bam or as a list of files", shortName = "I", required = true)
  var bams: File = _

  @Input(doc = "variant calls to annotate", fullName = "variantVCF", shortName = "C", required = true)
  var variantVCF: File = _

  @Input(doc = "dbSNP annotations VCF file", fullName = "dbsnp", shortName = "D", required = false)
  var dbsnp: File = _

  @Output(doc = "annotated file to output", shortName = "o", required = true)
  var outputAnnotated: File = _

  @Output(doc = "Memory limit", fullName = "memoryLimit", shortName = "m", required = false)
  var memoryLimit = 3

  @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
  var annotation: List[String] = Nil

  @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
  var group: List[String] = Nil

  @Argument(fullName="requireExplicitAnnotations", shortName="requireExplicitAnnotations", doc="SUPPRESS the default option of using all annotations", required=false)
  var requireExplicitAnnotations: Boolean = false

  @Argument(fullName="downsample_to_coverage", shortName="dcov", doc="Per-sample downsampling to perform", required=false)
  var downsample_to_coverage: Int = 0

  def script = {
    var extractIntervals : VCFExtractIntervals = new VCFExtractIntervals(qscript.variantVCF, swapExt(qscript.variantVCF, ".vcf", ".intervals.list"), true)
    add(extractIntervals)

    trait CommandLineGATKArgs extends CommandLineGATK {
      this.intervals :+= extractIntervals.listOut

      this.jarFile = qscript.gatkJarFile
      this.reference_sequence = qscript.referenceFile
      this.input_file = List(qscript.bams)

      this.memoryLimit = qscript.memoryLimit
      this.logging_level = "INFO"
    }

    class ScatteredFullVariantAnnotator() extends org.broadinstitute.sting.queue.extensions.gatk.VariantAnnotator with CommandLineGATKArgs {
      this.scatterCount = qscript.scatterCount
      this.variant = qscript.variantVCF

      this.useAllAnnotations = !qscript.requireExplicitAnnotations
      this.annotation = qscript.annotation
      this.group = qscript.group

      this.dbsnp = qscript.dbsnp

      this.out = qscript.outputAnnotated

      if (qscript.downsample_to_coverage > 0) {
        this.downsample_to_coverage = qscript.downsample_to_coverage
        this.downsampling_type = org.broadinstitute.sting.gatk.DownsampleType.BY_SAMPLE
      }
    }

    add(new ScatteredFullVariantAnnotator())
  }
}