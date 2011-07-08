package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport


  // ToDos:
  // reduce the scope of the datasets so the script is more nimble
  // create gold standard BAQ'd bam files, no reason to always do it on the fly

  // Analysis to add at the end of the script:
  // auto generation of the cluster plots
  // spike in NA12878 to the exomes and to the lowpass, analysis of how much of her variants are being recovered compared to single sample exome or HiSeq calls
  // produce Kiran's Venn plots based on comparison between new VCF and gold standard produced VCF


class MethodsDevelopmentCallingPipeline extends QScript {
  qscript =>

  @Argument(shortName="gatk", doc="gatk jar file", required=true)
  var gatkJarFile: File = _
  

  val hg19 = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")  
  val hg18 = new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
  val b36 = new File("/humgen/1kg/reference/human_b36_both.fasta")
  val b37 = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  val dbSNP_hg18_129 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_hg18.rod"            // Special case for NA12878 collections that can't use 132 because they are part of it.
  val dbSNP_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b36.rod"
  val dbSNP_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_132_b37.leftAligned.vcf"
  val dbSNP_b37_129 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf"              // Special case for NA12878 collections that can't use 132 because they are part of it.
  val hapmap_hg18 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.hg18_fwd.vcf"
  val hapmap_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b36_fwd.vcf"
  val hapmap_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
  val training_hapmap_b37 = "/humgen/1kg/processing/pipeline_test_bams/hapmap3.3_training_chr20.vcf"
  val omni_b36 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b36.vcf"
  val omni_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf"
  val indelMask_b36 = "/humgen/1kg/processing/pipeline_test_bams/pilot1.dindel.mask.b36.bed"
  val indelMask_b37 = "/humgen/1kg/processing/pipeline_test_bams/pilot1.dindel.mask.b37.bed"
  val exome_targets = "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"
  val training_1kg_high_quality = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"
  val training_1kg_terrible = "/humgen/1kg/processing/production_wgs_phase1/consensus_wgs/v2b/recal/phase1.wgs.projectConsensus.v2b.recal.terrible.vcf"

  val queueLogDir = ".qlog/"

  def script = {

    for (n <- List(1,5,10,20,30,40,50,60,70,80,96)) {
        var t = new File(n + ".tranches")
        var r = new File(n + ".recal")
        add(new VQSR(n, t, r))
        add(new applyVQSR(n, t, r))
  }
}

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    jarFile = gatkJarFile;
    memoryLimit = 4;
  }

  // 3.) Variant Quality Score Recalibration - Generate Recalibration table
  class VQSR(n: Int, t: File, r: File) extends VariantRecalibrator with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 4
    this.reference_sequence = hg19
    this.intervalsString ++= List(exome_targets)
    this.rodBind :+= RodBind("input", "VCF", new File("/broad/shptmp/ebanks/WEx1000G_forRyan/Barcoded_1000G_WEx_Plate_1.cleaned.annotated.snps." + n + "sample" + (if(n>1) {"s"} else {""}) +".vcf") )
    this.rodBind :+= RodBind("hapmap", "VCF", hapmap_b37, "known=false,training=true,truth=true,prior=15.0")
    this.rodBind :+= RodBind("omni", "VCF", omni_b37, "known=false,training=true,truth=false,prior=12.0")
    this.rodBind :+= RodBind("1kg", "VCF", training_1kg_high_quality, "known=false,training=true,truth=false,prior=12.0")
    this.rodBind :+= RodBind("badSites", "VCF", training_1kg_terrible, "known=false,training=false,truth=false,prior=2.0,bad=true")
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP_b37, "known=true,training=false,truth=false,prior=4.0")
    this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS")
    if(n >= 10) {
      this.use_annotation ++= List("InbreedingCoeff")
    }
    this.tranches_file = t
    this.recal_file = r
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
    this.analysisName = n + "_VQSR"
    this.jobName =  queueLogDir + n + ".VQSR"
    this.mG = 5
    this.std = 14
    this.percentBad = 0.04
    this.nt = 5
}

  // 4.) Apply the recalibration table to the appropriate tranches
  class applyVQSR (n: Int, t: File, r: File) extends ApplyRecalibration with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 4
    this.reference_sequence = hg19
    this.intervalsString ++= List(exome_targets)
    this.rodBind :+= RodBind("input", "VCF", new File("/broad/shptmp/ebanks/WEx1000G_forRyan/Barcoded_1000G_WEx_Plate_1.cleaned.annotated.snps." + n + "sample" + (if(n>1) {"s"} else {""}) +".vcf") )
    this.tranches_file = t
    this.recal_file = r
    this.ts_filter_level = 99.0
    this.out = new File("snps." + n + "sample" + (if(n>1) {"s"} else {""}) +".vcf")
    this.analysisName = n + "_AVQSR"
    this.jobName =  queueLogDir + n + ".applyVQSR"
  }

}
