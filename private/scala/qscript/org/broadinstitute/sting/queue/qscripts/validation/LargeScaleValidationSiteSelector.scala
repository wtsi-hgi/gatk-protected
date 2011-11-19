package org.broadinstitute.sting.queue.qscripts.validation


import org.broadinstitute.sting.pipeline.Pipeline
import org.broadinstitute.sting.queue.{QException, QScript}
import org.broadinstitute.sting.queue.extensions.gatk._

class LargeScaleValidationSiteSelector extends QScript {
  qscript =>

  @Input(doc="path to GATK jar", shortName="gatk", required=true)
  var gatkJar: File = _

  @Input(doc="Sample file", shortName="sampleFile", required=false)
  var sampleFile: File = new File("/humgen/gsa-hpprojects/dev/largeScaleValidation/validationSamples.txt")

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

    // SNP selector: 8000 based on AF
    selectSites(8000, "KEEP_AF_SPECTRUM","SNP", false)

    // SNP Selector: 8000 sampled uniformly, polymorphic in 8 validation samples
    selectSites(8000, "UNIFORM","SNP", true)

    // Indel selector: 5000 based on AF
    selectSites(5000, "KEEP_AF_SPECTRUM","INDEL", false)

    // Indel Selector: 5000 sampled uniformly, polymorphic in 8 validation samples
    selectSites(5000, "UNIFORM","INDEL", true)



  }

  def selectSites( numSites: Int, afMode: String, varType: String, useSampleFile: Boolean) = {

    // bind as input all Phase 1 autosome files
    var chrList =  List(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

    var genomeLength: Double = 0.0;

    for (len <- qscript.chromosomeLength) {
      genomeLength += len
    }

    // determine how many sites to validate per chromosome so that total length matches requested sites
    var sitesToValidatePerChr = new Array[Int](22)

    var k:Int = 0
    var totalSites:Int = 0
    var lens = qscript.chromosomeLength.toArray

    //logger.debug("numSites:%4.1f".format(numSites.toDouble))
    for(chr <- chrList) {
      //logger.debug("lens(k)=%d, (numSites.toDouble * lens(k).toDouble )/genomeLength.toDouble=%4.1f".format(lens(k),(numSites.toDouble * lens(k).toDouble)/genomeLength.toDouble))
      sitesToValidatePerChr(k) = math.round(  (numSites.toDouble * lens(k).toDouble )/genomeLength.toDouble).toInt
      //logger.debug("sitesToValidatePerChr(k)=%d".format(sitesToValidatePerChr(k)))
      totalSites +=  sitesToValidatePerChr(k)
      k=k+1
    }
   //logger.debug("totalSites=%d".format(totalSites))
    // tmp hack: whatever rounding difference we have gets binned in a single random chromosome
    val rnd = new scala.util.Random

    sitesToValidatePerChr(rnd.nextInt(sitesToValidatePerChr.length)) += (numSites - totalSites)

    var combined = new CombineVariants with CommandLineGATKArgs
    combined.out = qscript.outputDir + "/"+"ALL.wgs.%d_validation_sites_%s.%s.sites.vcf".format(numSites, afMode, varType)
    combined.jobOutputFile = combined.out + ".out"

    for(chr <- chrList) {
      // no X chr in phase 1 official release yet
      //  if (chr == 23)
      //    chrStr = "X"
      var selector = new ValidationSiteSelector with CommandLineGATKArgs

      selector.numSites = sitesToValidatePerChr(chr-1)

      if (varType.matches("INDEL"))
      { selector.selectType:+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.INDEL}
      else {
        selector.selectType:+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.SNP
      }

      if (afMode.matches("UNIFORM"))
        selector.freqMode = org.broadinstitute.sting.gatk.walkers.ValidationSiteSelector.ValidationSiteSelectorWalker.AF_COMPUTATION_MODE.UNIFORM
      else
        selector.freqMode = org.broadinstitute.sting.gatk.walkers.ValidationSiteSelector.ValidationSiteSelectorWalker.AF_COMPUTATION_MODE.KEEP_AF_SPECTRUM

      if (useSampleFile)
        selector.sample_file:+= qscript.sampleFile

      selector.intervalsString :+= chr.toString
      selector.variant :+= new File("/humgen/1kg/DCC/ftp/release/20110521/ALL.chr" +chr.toString + ".merged_beagle_mach.20101123.snps_indels_svs.genotypes.vcf.gz")
      selector.sampleMode = org.broadinstitute.sting.gatk.walkers.ValidationSiteSelector.ValidationSiteSelectorWalker.SAMPLE_SELECTION_MODE.POLY_BASED_ON_GT

      selector.out = qscript.tmpDir + "/"+"ALL.chr%d.%d_%s_%s.sites.vcf".format(chr, selector.numSites, varType, afMode)
      selector.jobOutputFile = selector.out + ".out"
      combined.variant :+= selector.out

      add(selector)

    }
    add(combined)

  }

}
