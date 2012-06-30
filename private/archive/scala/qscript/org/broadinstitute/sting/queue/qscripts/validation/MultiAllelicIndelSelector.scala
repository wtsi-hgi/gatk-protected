package org.broadinstitute.sting.queue.qscripts.validation

import org.broadinstitute.sting.pipeline.Pipeline
import org.broadinstitute.sting.queue.{QException, QScript}
import org.broadinstitute.sting.queue.extensions.gatk._


class MultiAllelicIndelSelector extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="Sample file", shortName="sampleFile", required=false)
  var sampleFile: File = new File("/humgen/gsa-hpprojects/dev/largeScaleValidation/validationSamples.txt")

  @Input(doc="Sample file", shortName="inputFile", required=false)
  var inputFile: File = new File("/humgen/1kg/processing/official_release/phase1/projectConsensus/ALL.wgs.VQSR_V2_GLs_polarized.20101123.indels.genotypes.vcf.gz")

  @Input(doc="Number of Sites", shortName="numSites", required=false)
  var numSites: Int = 2000

  @Input(doc="output path", shortName="outputDir", required=false)
  var outputDir: String =  "/humgen/gsa-hpprojects/dev/largeScaleValidation/outputVCFs/"

  private val tmpDir: File = new File("/broad/shptmp/delangel/tmp/")
  private val reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  private var pipeline: Pipeline = _
  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = 2
    this.jobTempDir = qscript.tmpDir
    this.jobQueue = "gsa";
  }

  def script = {

    // bind as input all Phase 1 autosome files
    var chrList =  List(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

    var selector = new ValidationSiteSelector with CommandLineGATKArgs
    selector.out = qscript.outputDir + "/"+"ALL.wgs.%d_validation_sites_multiAllelicIndels.sites.vcf".format(numSites)
    selector.jobOutputFile = selector.out + ".out"
    selector.sampleMode = org.broadinstitute.sting.gatk.walkers.ValidationSiteSelector.ValidationSiteSelectorWalker.SAMPLE_SELECTION_MODE.POLY_BASED_ON_GT
    selector.freqMode = org.broadinstitute.sting.gatk.walkers.ValidationSiteSelector.ValidationSiteSelectorWalker.AF_COMPUTATION_MODE.UNIFORM
    selector.numSites = qscript.numSites
    for(chr <- chrList) {
      // no X chr in phase 1 official release yet
      //  if (chr == 23)
      //    chrStr = "X"
      var vs = new SelectVariants with CommandLineGATKArgs
      vs.intervalsString :+= chr.toString
      vs.variant = qscript.inputFile
      vs.sample_file :+= qscript.sampleFile
      vs.ef = true
      vs.env = true
      vs.selectType:+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.INDEL
      vs.out  = qscript.tmpDir + "/"+"ALL.chr%d.%d_multiAllelicIndels.sites.vcf".format(chr, selector.numSites)
      vs.restrictAllelesTo = org.broadinstitute.sting.gatk.walkers.variantutils.SelectVariants.NumberAlleleRestriction.MULTIALLELIC

      selector.variant :+= vs.out


      add(vs)

    }
    add(selector)

  }

}
