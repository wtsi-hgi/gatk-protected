package org.broadinstitute.sting.queue.batchMerges.BatchMerging

import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.library.ipf.vcf.{VCFSimpleMerge, VCFExtractSites,VCFExtractIntervals}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.baq.BAQ
import org.broadinstitute.sting.utils.text.XReadLines
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.gatk.walkers.genotyper.{ExactAFCalculationModel, GenotypeLikelihoodsCalculationModel, UnifiedGenotyperEngine}

class BatchedCallUnionMerger extends QScript {
  batchMerge =>

  @Argument(doc = "gatk jar file", shortName = "J", required = true)
  var gatkJarFile: File = _

  @Argument(shortName = "R", doc = "ref", required = true)
  var referenceFile: File = _

  @Argument(doc="VCF list",shortName="vcfs")
  var vcfList: File = _

  @Input(doc = "bam input, as .bam or as a list of files", shortName = "I", required = true)
  var bamList: File = _

  @Argument(doc="batched output",shortName="o")
  var batchOut: File = _

  //@Argument(doc="read UG settings from header",shortName="ugh") var ugSettingsFromHeader : Boolean = false

  @Hidden
  @Argument(doc="Min base q",shortName="mbq",required=false)
  var mbq : Int = 17

  @Hidden
  @Argument(doc="baq gap open penalty, using sets baq to calc when necessary",shortName="baqp",required=false)
  var baq : Int = -1

  @Argument(fullName="downsample_to_coverage", shortName="dcov", doc="Per-sample downsampling to perform", required=false)
  var downsample_to_coverage: Int = 0

  @Argument(fullName="annotate", shortName="annotate", doc="Annotate output VCF", required=false)
  var annotate: Boolean = false

  @Input(doc = "level of parallelism. By default is set to 0 [no scattering].", shortName = "scatter", required = false)
  var scatterCount = 0

  @Input(doc = "dbSNP annotations VCF file", fullName = "dbsnp", shortName = "D", required = false)
  var dbsnp: File = _

  @Output(doc = "Memory limit", fullName = "memoryLimit", shortName = "m", required = false)
  var memoryLimit = 3

  @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
  var annotation: List[String] = Nil

  @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", required=false)
  var excludeAnnotation: List[String] = Nil

  @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
  var group: List[String] = Nil

  @Argument(fullName="requireExplicitAnnotations", shortName="requireExplicitAnnotations", doc="SUPPRESS the default option of using all annotations", required=false)
  var requireExplicitAnnotations: Boolean = false

  def script = {

    var vcfs : List[File] = extractFileEntries(vcfList)
    var bams : List[File] = extractFileEntries(bamList)

    trait ExtractArgs extends VCFExtractSites {
      this.keepFilters = true
      this.keepInfo = true
      this.keepQual = true
    }

    trait CombineVariantsArgs extends CombineVariants {
      this.reference_sequence = batchMerge.referenceFile
      this.jarFile = batchMerge.gatkJarFile
      this.scatterCount = 10
      this.memoryLimit = 4
      this.filteredrecordsmergetype = VariantContextUtils.FilteredRecordMergeType.KEEP_UNCONDITIONAL;
      this.multipleAllelesMergeType = VariantContextUtils.MultipleAllelesMergeType.MIX_TYPES;
      this.setKey = "set";
    }

    var combine : CombineVariants = new CombineVariants with CombineVariantsArgs
    combine.out = swapExt(batchOut,".vcf",".variant.combined.vcf")
    combine.variant ++= vcfs.map( u => new TaggedFile(u, "VCF") )
    add(combine)

    var getVariantAlleles : List[VCFExtractSites] = vcfs.map( u => new VCFExtractSites(u, swapExt(batchOut.getParent,u,".vcf",".alleles.vcf")) with ExtractArgs)
    var combineVCFs : VCFSimpleMerge = new VCFSimpleMerge
    combineVCFs.vcfs = getVariantAlleles.map(u => u.outVCF)
    combineVCFs.fai = new File(referenceFile.getAbsolutePath+".fai")
    combineVCFs.outVCF = swapExt(batchOut,".vcf",".pf.alleles.vcf")

    var extractIntervals : VCFExtractIntervals = new VCFExtractIntervals(combine.out,swapExt(combine.out,".vcf",".intervals.list"),true)
    //addAll(getVariantAlleles)
    //add(combineVCFs,extractIntervals)
    add(extractIntervals)

    trait CalcLikelihoodArgs extends UGCalcLikelihoods {
      this.reference_sequence = batchMerge.referenceFile
      this.min_base_quality_score = batchMerge.mbq
      if ( batchMerge.baq >= 0 ) {
        this.baqGapOpenPenalty = batchMerge.baq
        this.baq = BAQ.CalculationMode.CALCULATE_AS_NECESSARY
      }
      this.intervals :+= extractIntervals.listOut
      this.alleles = new TaggedFile(combine.out, "VCF")
      this.jarFile = batchMerge.gatkJarFile
      this.memoryLimit = 4
      this.scatterCount = 60
      this.output_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
      this.genotyping_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES

      if (batchMerge.downsample_to_coverage > 0) {
        this.downsample_to_coverage = batchMerge.downsample_to_coverage
        this.downsampling_type = org.broadinstitute.sting.gatk.DownsampleType.BY_SAMPLE
      }

      this.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.BOTH;
    }

    def newUGCL( bams: (List[File],Int) ) : UGCalcLikelihoods = {
      var ugcl = new UGCalcLikelihoods with CalcLikelihoodArgs
      ugcl.input_file ++= bams._1
      ugcl.out = new File("MBatch%d.likelihoods.vcf".format(bams._2))
      return ugcl
    }

    var calcs: List[UGCalcLikelihoods] = bams.grouped(20).toList.zipWithIndex.map(u => newUGCL(u))
    addAll(calcs)

    trait CallVariantsArgs extends UGCallVariants {
      this.reference_sequence = batchMerge.referenceFile
      this.intervals :+= extractIntervals.listOut
      this.jarFile = batchMerge.gatkJarFile
      this.scatterCount = 30
      this.output_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
      this.genotyping_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES

      if (batchMerge.downsample_to_coverage > 0) {
        this.downsample_to_coverage = batchMerge.downsample_to_coverage
        this.downsampling_type = org.broadinstitute.sting.gatk.DownsampleType.BY_SAMPLE
      }

      this.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.BOTH;

      // The memory-intensive part is limited by the number of ALT alleles:
      this.memoryLimit = 16
      this.max_alternate_alleles = 2
    }

    var callVariants : UGCallVariants = new UGCallVariants with CallVariantsArgs
    callVariants.variant ++= calcs.map( a => new TaggedFile(a.out, "VCF,custom=variant" + a.out.getName.replace(".vcf","")) )
    callVariants.alleles = new TaggedFile(combine.out, "VCF")

    callVariants.out = batchOut
    if (batchMerge.annotate) {
      callVariants.out = new File("unannotated." + batchOut.getName)
    }

    add(callVariants)

    if (batchMerge.annotate) {
      trait CommandLineGATKArgs extends CommandLineGATK {
        this.intervals :+= extractIntervals.listOut

        this.jarFile = batchMerge.gatkJarFile
        this.reference_sequence = batchMerge.referenceFile
        this.input_file = List(batchMerge.bamList)

        this.memoryLimit = batchMerge.memoryLimit
        this.logging_level = "INFO"
      }

      class ScatteredFullVariantAnnotator(inputParam: File) extends org.broadinstitute.sting.queue.extensions.gatk.VariantAnnotator with CommandLineGATKArgs {
        this.scatterCount = batchMerge.scatterCount
        this.variant = inputParam

        this.useAllAnnotations = !batchMerge.requireExplicitAnnotations
        this.annotation = batchMerge.annotation
        this.excludeAnnotation = batchMerge.excludeAnnotation
        this.group = batchMerge.group

        this.dbsnp = batchMerge.dbsnp

        this.out = batchMerge.batchOut

        if (batchMerge.downsample_to_coverage > 0) {
          this.downsample_to_coverage = batchMerge.downsample_to_coverage
          this.downsampling_type = org.broadinstitute.sting.gatk.DownsampleType.BY_SAMPLE
        }
      }
      
      var annotateVariants = new ScatteredFullVariantAnnotator(callVariants.out)
      add(annotateVariants)
    }
  }

  def extractFileEntries(in: File): List[File] = {
    return (new XReadLines(in)).readLines.toList.map( new File(_) )
  }
}
