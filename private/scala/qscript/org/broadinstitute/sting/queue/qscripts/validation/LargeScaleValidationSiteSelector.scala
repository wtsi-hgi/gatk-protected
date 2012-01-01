package org.broadinstitute.sting.queue.qscripts.validation


import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.{QException, QScript}

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

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.reference
    this.memoryLimit = 32
    this.jobTempDir = qscript.tmpDir
    this.jobQueue = "gsa";
    /*    this.intervalsString =List("/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list")
     if ( !qscript.intervals.isEmpty)
       this.intervalsString :+= (qscript.intervals)
    */
    //this.DBSNP = qscript.dbSNP
  }
  /*
  trait ExpandedIntervals extends CommandLineGATK {
    if (qscript.expandIntervals > 0)
      this.intervals :+= flankIntervals
  }
     */
  def script = {

    // SNP selector: 8000 based on AF
    var AFSnps = new SiteSelector
    AFSnps.selectType:+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.SNP
    AFSnps.freqMode = org.broadinstitute.sting.gatk.walkers.ValidationSiteSelector.ValidationSiteSelectorWalker.AF_COMPUTATION_MODE.KEEP_AF_SPECTRUM
    AFSnps.numSites = 8000
    AFSnps.out = qscript.outputDir + "/"+"ALL.wgs.%d_validation_sites_AF_distributed.snp.sites.vcf".format(AFSnps.numSites)
    AFSnps.jobOutputFile = AFSnps.out + ".out"
    //add(AFSnps)

    // SNP Selector: 8000 sampled uniformly, polymorphic in 8 validation samples
    var USnps = new SiteSelector
    AFSnps.selectType:+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.SNP
    USnps.numSites = 8000
    USnps.freqMode = org.broadinstitute.sting.gatk.walkers.ValidationSiteSelector.ValidationSiteSelectorWalker.AF_COMPUTATION_MODE.UNIFORM
    USnps.sample_file:+= qscript.sampleFile
    USnps.out = qscript.outputDir + "/"+"ALL.wgs.%d_validation_sites_Uniformly_distributed.snp.sites.vcf".format(USnps.numSites)
    USnps.jobOutputFile = USnps.out + ".out"
    //add(USnps)

    // Indel selector: 5000 based on AF
    var AFIndels = new SiteSelector
    AFIndels.selectType:+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.INDEL
    AFIndels.freqMode = org.broadinstitute.sting.gatk.walkers.ValidationSiteSelector.ValidationSiteSelectorWalker.AF_COMPUTATION_MODE.KEEP_AF_SPECTRUM
    AFIndels.numSites = 5000
    AFIndels.out = qscript.outputDir + "/"+"ALL.wgs.%d_validation_sites_AF_distributed.indels.sites.vcf".format(AFIndels.numSites)
    AFIndels.jobOutputFile = AFIndels.out + ".out"
    AFIndels.numBins = 50
    add(AFIndels)

    // Indel Selector: 5000 sampled uniformly, polymorphic in 8 validation samples
    var UIndels = new SiteSelector
    UIndels.selectType:+= org.broadinstitute.sting.utils.variantcontext.VariantContext.Type.INDEL
    UIndels.numSites = 5000
    UIndels.freqMode = org.broadinstitute.sting.gatk.walkers.ValidationSiteSelector.ValidationSiteSelectorWalker.AF_COMPUTATION_MODE.UNIFORM
    UIndels.sample_file:+= qscript.sampleFile
    UIndels.out = qscript.outputDir + "/"+"ALL.wgs.%d_validation_sites_Uniformly_distributed.indels.sites.vcf".format(UIndels.numSites)
    UIndels.jobOutputFile = UIndels.out + ".out"
    UIndels.numBins = 50
    add(UIndels)




  }
  class SiteSelector() extends ValidationSiteSelector with CommandLineGATKArgs {

    // bind as input all Phase 1 autosome files
    var chrList =  List(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

    for(chr <- chrList) {
      // no X chr in phase 1 official release yet
      //  if (chr == 23)
      //    chrStr = "X"

      this.variant :+= new File("/humgen/1kg/DCC/ftp/release/20110521/ALL.chr" +chr.toString + ".merged_beagle_mach.20101123.snps_indels_svs.genotypes.vcf.gz")
    }
    this.sampleMode = org.broadinstitute.sting.gatk.walkers.ValidationSiteSelector.ValidationSiteSelectorWalker.SAMPLE_SELECTION_MODE.POLY_BASED_ON_GT

  }

}
