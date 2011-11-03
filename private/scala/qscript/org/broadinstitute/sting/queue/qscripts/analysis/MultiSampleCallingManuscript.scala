/*
 * Copyright (c) 2011, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.queue.qscripts.analysis

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport
import java.io.FileWriter
import org.broadinstitute.sting.utils.exceptions.UserException

class MultiSampleCallingManuscript extends QScript {
  qscript =>

  @Argument(shortName="outputDir", doc="output directory", required=false)
  var outputDir: String = "multiSampleManuscriptResults"

  @Argument(shortName="extraSamples", doc="BAM list of extra samples", required=true)
  var extraSamples: File = _

  @Argument(shortName="sample", doc="Samples to include in Variant Eval", required=false)
  var samples: List[String] = List("NA12878")

  @Argument(shortName = "L", fullName = "intervals", doc="intervals", required=false)
  val myIntervals: List[String] = null;

  @Argument(shortName = "sc", fullName = "scatterCount", doc = "Scatter/Gather jobs to use", required=false)
  val myScatterCount: Int = 1;

//  val b37_decoy = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  val bundle = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/")
  val b37 = new File(bundle.getPath + "/human_g1k_v37.fasta")
  val dbSNP_b37 = new File(bundle.getPath + "/dbsnp_132.b37.vcf")
  val dbSNP_b37_129 = new File(bundle.getPath + "/dbsnp_132.b37.excluding_sites_after_129.vcf")
  val hapmap_b37 = new File(bundle.getPath + "/hapmap_3.3.b37.vcf")
  val omni_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf"
  val training_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf"
  val badSites_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.terrible.vcf"
  val projectConsensus_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/ALL.wgs.projectConsensus_v2b.20101123.snps.sites.vcf"
  val queueLogDir = ".qlog/"

  val TITV_EXPECTED = 2.1
  // todo -- should really not be decoy but waiting for newest bwa results
  val NA12878BAM = new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam")
  val NA12878_AVERAGE_DEPTH = 64

  class Target(
          val baseName: String,
          val bamList: File,
          val titvTarget: Double,
          val trancheTarget: Double,
          val nSamples: Int,
          val na12878Depth: Int,
          val nAdditionalSamples: Int,
          val additionalSamples: File) {
    val name = "%s%s_d%d_w%d".format(qscript.outputDir, baseName, na12878Depth, nAdditionalSamples)
    val clusterFile = new File(name + ".clusters")
    val rawVCF = new File(name + ".raw.vcf")
    val recalibratedVCF = new File(name + ".recalibrated.vcf")
    val tranchesFile = new File(name + ".tranches")
    val vqsrRscript = name + ".vqsr.r"
    val recalFile = new File(name + ".tranches.recal")
    val evalFile = new File(name + ".snp.eval")
    val evalIndelFile = new File(name + ".indel.eval")
    def useBAQ = true
    def needToDownsampleFraction = na12878Depth < NA12878_AVERAGE_DEPTH
    def downsampleFraction = na12878Depth / (1.0*NA12878_AVERAGE_DEPTH)
  }

  def createTargets: List[Target] = {
    var targets: List[Target] = List()
    for ( nAdditionalSamples <- List(0, 10, 100, 300) ) {
      for ( na12878Depth <- List(4, 30, 64) ) {
        for ( sensitivityThreshold <- List(90.0, 99.0) ) {
          val file = if ( nAdditionalSamples > 0 ) createAdditionalSampleFile(nAdditionalSamples) else null
          val t = new Target("NA12878", NA12878BAM, TITV_EXPECTED, sensitivityThreshold, 1, na12878Depth, nAdditionalSamples, file)
          targets :+= t
        }
      }
    }

    targets
  }

  def createAdditionalSampleFile(nAdditionalSamples: Int): File = {
    val lines = scala.io.Source.fromFile(extraSamples).getLines().slice(0, nAdditionalSamples).toList

    if ( nAdditionalSamples > lines.size )
      throw new UserException("not enough extra samples provided")

    val destFile = new File(outputDir + "/additional_samples_%d.bam.list".format(nAdditionalSamples))
    val fw = new FileWriter(destFile)
    for ( line <- lines ) fw.write("%s%n".format(line))
    fw.close()

    destFile
  }


  def script = {
    if ( ! outputDir.exists() ) outputDir.mkdirs()
    val targets = createTargets

    for (target <- targets) {
      add(new callVariants(target))
      //add(new indelFilter(target), new indelEvaluation(target))
      add(new VQSR(target))
      add(new applyVQSR(target))
      add(new evalVariants(target))
    }
  }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    logging_level = "INFO";
    memoryLimit = 4;
    reference_sequence = b37
    intervalsString = myIntervals
  }

  // 1.) Unified Genotyper Base
  class callVariants (t: Target) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.scatterCount = myScatterCount
    //this.nt = 2
    this.dcov = 250
    this.stand_call_conf = 30.0
    this.stand_emit_conf = 30.0
    if ( t.needToDownsampleFraction ) this.downsample_to_fraction = t.downsampleFraction
    this.input_file :+= t.bamList
    this.input_file :+= t.additionalSamples
    this.D = dbSNP_b37
    this.out = t.rawVCF
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.baq = if (t.useBAQ) {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY} else {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF}
    this.analysisName = t.name + "_UGs"
    this.jobName =  queueLogDir + t.name + ".snpcall"
  }

//  // todo -- I'd like this to be one process -> snp + indel -> snp + indel (filtered) -> snp (filtered) + indel (filtered)
//  // 2.) Hard Filtering for indels
//  class indelFilter (t: Target) extends VariantFiltration with UNIVERSAL_GATK_ARGS {
//    this.memoryLimit = 2
//    this.scatterCount = 10
//    this.V = t.rawIndelVCF
//    this.out = t.filteredIndelVCF
//    this.filterName ++= List("IndelQD", "IndelReadPosRankSum", "IndelFS")
//    this.filterExpression ++= List("\"QD < 2.0\"", "\"ReadPosRankSum < -20.0\"", "\"FS > 200.0\"")
//    if (t.nSamples >= 10) {
//        this.filterName ++= List("IndelInbreedingCoeff")
//        this.filterExpression ++= List("\"InbreedingCoeff < -0.8\"")
//    }
//    this.analysisName = t.name + "_VF"
//    this.jobName =  queueLogDir + t.name + ".indelfilter"
//  }

  // 3.) Variant Quality Score Recalibration - Generate Recalibration table
  class VQSR(t: Target) extends VariantRecalibrator with UNIVERSAL_GATK_ARGS {
    this.nt = 2
    this.input :+= t.rawVCF
    this.resource :+= new TaggedFile( hapmap_b37, "training=true,truth=true,prior=15.0" )
    this.resource :+= new TaggedFile( omni_b37, "training=true,truth=true,prior=12.0" )
    this.resource :+= new TaggedFile( training_1000G, "training=true,prior=10.0" )
    this.resource :+= new TaggedFile( dbSNP_b37, "known=true,prior=2.0" )
    this.resource :+= new TaggedFile( projectConsensus_1000G, "prior=8.0" )
    this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS")
    if(t.nSamples >= 10) { // InbreedingCoeff is a population-wide statistic that requires at least 10 samples to calculate
        this.use_annotation ++= List("InbreedingCoeff")
    }
    this.use_annotation ++= List("DP")
    //this.resource :+= new TaggedFile( badSites_1000G, "bad=true,prior=2.0" )
    this.tranches_file = t.tranchesFile
    this.recal_file = t.recalFile
    this.allPoly = true // TODO -- what does this do?
    this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
    this.rscript_file = t.vqsrRscript
    this.analysisName = t.name + "_VQSR"
    this.jobName = queueLogDir + t.name + ".VQSR"
  }

  // 4.) Apply the recalibration table to the appropriate tranches
  class applyVQSR (t: Target) extends ApplyRecalibration with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 6
    this.input :+= t.rawVCF
    this.tranches_file = t.tranchesFile
    this.recal_file = t.recalFile
    this.ts_filter_level = t.trancheTarget
    this.out = t.recalibratedVCF
    this.analysisName = t.name + "_AVQSR"
    this.jobName = queueLogDir + t.name + ".applyVQSR"
  }

  // 5.) Variant Evaluation Base(OPTIONAL)
  class evalVariants(t: Target) extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.comp :+= new TaggedFile(hapmap_b37, "hapmap" )
    this.D = new File(dbSNP_b37_129)
    this.sample = samples
    this.comp :+= new TaggedFile( omni_b37, "omni" )
    this.eval :+= t.recalibratedVCF
    this.out =  t.evalFile
    this.analysisName = t.name + "_VEs"
    this.evalModule :+= "IndelStatistics"
    this.jobName = queueLogDir + t.name + ".snp.eval"
  }
}
