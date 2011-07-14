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
  val exome_flanks = "/humgen/gsa-hpprojects/dev/rpoplin/exome_t2d/GoT2D_exomes_batch_005.cleaned.flanks.interval_list"
  val training_1kg_high_quality = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"
  val training_1kg_terrible = "/humgen/1kg/processing/production_wgs_phase1/consensus_wgs/v2b/recal/phase1.wgs.projectConsensus.v2b.recal.terrible.vcf"

  val queueLogDir = ".qlog/"

  def script = {

    for (n <- List(1,5,10,20,30,40,50,60,70,80,96)) {
        var v = new File("t2d." + n + ".vcf")
        var b = new File("lists/t2d." + n + ".bam.list")
        var t = new File("t2d." + n + ".tranches")
        var r = new File("t2d." + n + ".recal")
        var o = new File("t2d.recal." + n + ".vcf")
        var s = new File("t2d.recal.selected." + n + ".vcf")
        add(new UG(n, b, v, t, r, o, s))
        add(new VQSR(n, b, v, t, r, o, s))
        add(new ApplyVQSR(n, b, v, t, r, o, s))
        add(new Select(n, b, v, t, r, o, s))

        v = new File("t2d.1kg." + n + ".vcf")
        b = new File("lists/t2d.1kg." + n + ".bam.list")
        t = new File("t2d.1kg." + n + ".tranches")
        r = new File("t2d.1kg." + n + ".recal")
        o = new File("t2d.1kg.recal." + n + ".vcf")
        s = new File("t2d.1kg.recal.selected." + n + ".vcf")
        add(new UG(n, b, v, t, r, o, s))
        add(new VQSR(n, b, v, t, r, o, s))
        add(new ApplyVQSR(n, b, v, t, r, o, s))
        add(new Select(n, b, v, t, r, o, s))
  }
}

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    jarFile = gatkJarFile;
    reference_sequence = hg19
    intervalsString ++= List(exome_targets, exome_flanks)
  }

  // 1.) Unified Genotyper
  class UG(n: Int, b: File, v: File, t: File, r: File, o: File, s: File) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.scatterCount = 120
    this.dcov = 250
    this.stand_call_conf = 30.0
    this.stand_emit_conf = 30.0
    this.input_file :+= b
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP_b37)
    this.out = v
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
  }

  // 3.) Variant Quality Score Recalibration - Generate Recalibration table
  class VQSR(n: Int, b: File, v: File, t: File, r: File, o: File, s: File) extends VariantRecalibrator with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 4
    this.rodBind :+= RodBind("input", "VCF", v)
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
    this.mG = 6
    this.std = 14
    this.percentBad = 0.03
    this.nt = 2
}

  // 4.) Apply the recalibration table to the appropriate tranches
  class ApplyVQSR(n: Int, b: File, v: File, t: File, r: File, o: File, s: File) extends ApplyRecalibration with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 4
    this.rodBind :+= RodBind("input", "VCF", v)
    this.tranches_file = t
    this.recal_file = r
    this.ts_filter_level = 99.0
    this.out = o
  }

  class Select(n: Int, b: File, v: File, t: File, r: File, o: File, s: File) extends SelectVariants with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.rodBind :+= RodBind("variant", "VCF", o)
    this.out = s
    this.sf ++= List( new File("lists/t2d." + n + ".sample.list") )
  }
}
