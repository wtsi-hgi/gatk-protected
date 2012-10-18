/**
 * Created with IntelliJ IDEA.
 * User: kcibul
 * Date: 9/26/12
 */

package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript

import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.cancer.MuTect

class CMICancerPipeline extends QScript {
  qscript =>


  /****************************************************************************
    * Required executable locations
    ****************************************************************************/
  @Input(doc="The path to the binary of MuTect", fullName="mutect_jar", shortName="mj", required=true)
  var mutectJar: File = _

  @Input(doc="The path to the binary of IndelGenotyper", fullName="indel_jar", shortName="ij", required=true)
  var indelGenotyperJar: File = _

  @Input(doc="The path to the indel filter script", fullName="indel_filter", shortName="if", required=true)
  var indelFilterPath: File = _

  @Input(doc="The path to the germlin indel filter script", fullName="indel_germline_filter", shortName="igf", required=true)
  var indelGermlineFilterPath: File = _

  @Input(doc="The path to the indel caller refgene file", fullName="indel_refgene", shortName="ir", required=true)
  var indelRefGenePath: File = _

  @Input(doc="The path to the indel caller germline event file", fullName="indel_germline_database", shortName="igdb", required=true)
  var indelGermlineDatabase: File = _

  @Input(doc="The path to the indel caller maflite conversion script", fullName="indel_to_maflite", shortName="itm", required=true)
  var indelToMaflitePath: File = _

  /****************************************************************************
    * Required Parameters
    ****************************************************************************/
  @Input(doc="a BAM list of one or more raw BAMs", fullName="bam_list", shortName="bl", required=true)
  var bamList: File = _

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=true)
  var reference: File = _

  @Input(doc="DBSNP or known callset to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=true)
  var dbSNP: Seq[File] = Seq()

  @Input(doc="Intervals to process", fullName="intervals", shortName="L", required=true)
  var intervals: Seq[File] = Seq()

  @Input(doc="List of indices of tumor BAMs", fullName="tumor", shortName="t", required=true)
  var tumorBamIndicies: Seq[Int] = _


  /****************************************************************************
    * Optional Input Parameters
    ****************************************************************************/
  @Input(doc="Panel Of Normals or known artifact sites to use (must be in VCF format)", fullName="panel_of_normals", shortName="pon", required=false)
  var pon: Seq[File] = Seq()

  @Input(doc="COSMIC sites to use (must be in VCF format)", fullName="cosmic", shortName="C", required=false)
  var cosmic: Seq[File] = Seq()

  /****************************************************************************
    * Optional Input Parameters with sensible defaults
    ****************************************************************************/
  @Input(doc="Base filter for SomaticIndelDetector", fullName="indel_base_filter", shortName = "ibf", required=false)
  var indelCallerBaseFilter: String = "T_COV<6||N_COV<4||(T_INDEL_F<=0.3&&T_CONS_CNT<7)||T_INDEL_CF<=0.7"


  /****************************************************************************
    * Hidden Parameters
    ****************************************************************************/
  @Hidden
  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = 0

  @Hidden
  @Argument(doc="Run the pipeline in test mode only", fullName = "test_mode", shortName = "test", required=false)
  var testMode: Boolean = false


  /****************************************************************************
    * Global Variables
    ****************************************************************************/


  /****************************************************************************
    * Main script
    ****************************************************************************/

  def script() {

    val bams: Seq[File] = QScriptUtils.createSeqFromFile(bamList)


    // TODO: actually implement this!!  For now assume tumor is (0) and normal is (1);
    val tumorBams: Seq[File] = Seq(bams(0))
    val normalBams: Seq[File] = Seq(bams(1))

    assert(bams.length <= 5, "Current implementation is limited to 5 bams. See source code for details on the n-way out indel cleaning procedure")

    print("Indicies:" + tumorBamIndicies + "\n")
    print("Tumor:" + tumorBams + "\n")
    print("Normal:" + normalBams + "\n")

    for (tumorBam <- tumorBams) {
      for (normalBam <- normalBams) {
        val tumorName = swapExt(tumorBams(0).getName, ".bam", "")
        val normalName = swapExt(normalBams(0).getName, ".bam", "")

        val outPrefix = tumorName + "-vs-" + normalName // TODO: use CMI ids here
        cancer(tumorName, tumorBam, normalName, normalBam, 0.01f, outPrefix)
      }
    }
  }




  /****************************************************************************
    * Classes (GATK Walkers)
    ****************************************************************************/



  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 4
    this.isIntermediate = false
  }

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
  }

  def cancer (tumorName : String, tumorBam : File, normalName : String, normalBam : File, tumorFractionContamination : Float, outPrefix : String) {
    val rawMutations = outPrefix + ".call_stats.txt"
    val rawVcf = outPrefix + ".vcf"
    val rawCoverage = outPrefix + ".wig.txt"

    add(mutect(tumorName, tumorBam, normalName, normalBam, tumorFractionContamination, rawMutations, rawVcf, rawCoverage))
    indels(tumorName, tumorBam, normalName, normalBam, outPrefix)

  }

  def indels (tumorName : String, tumorBam : File, normalName : String, normalBam : File, outPrefix : String) {
    val rawIndels = outPrefix + ".raw.indels.txt"
    val nFilteredIndels = outPrefix + ".nfilter.indels.txt"
    val ntFilteredIndels = outPrefix + ".nfilter.tfilter.indels.txt"
    val ntgFilteredIndels = outPrefix + ".nfilter.tfilter.gfilter.indels.txt"
    val mafliteIndels = outPrefix + ".filtered.maflite"

    add(callIndels(tumorBam, normalBam, rawIndels))
    add(filterIndelsNormal(rawIndels, nFilteredIndels))
    add(filterIndelsTumor(nFilteredIndels, ntFilteredIndels))
    add(filterIndelsGermline(ntFilteredIndels, ntgFilteredIndels))
    add(filterIndelsGermline(ntFilteredIndels, ntgFilteredIndels))
    add(convertIndelCallsToMaflite(tumorName, normalName, ntgFilteredIndels, mafliteIndels))
  }

  case class mutect (tumorName : String, tumorBam : File, normalName : String, normalBam : File, tumorFractionContamination : Float, outMutations : File,  outVcf : File, outCoverage : File) extends MuTect with CommandLineGATKArgs {
    this.scatterCount = 1
    this.memoryLimit = 4
    this.jarFile = qscript.mutectJar
    this.intervals = qscript.intervals

    this.dbsnp = qscript.dbSNP
    this.cosmic = qscript.cosmic
    this.normal_panel = qscript.pon

    this.only_passing_calls = true
    this.enable_extended_output = true
    this.downsample_to_coverage = 1000 // TODO: how deep should this be?
    this.fraction_contamination = Some(tumorFractionContamination)

    this.input_file :+= new TaggedFile(tumorBam, "tumor")
    this.tumor_sample_name = tumorName

    this.input_file :+= new TaggedFile(normalBam, "normal")
    this.normal_sample_name = normalName

    this.out = outMutations
    this.coverage_file = outCoverage
    this.vcf = outVcf


    this.analysisName = tumorBam.toString + ".mutect"
    this.jobName = this.analysisName

    print ("MuTect CMD:" + this.commandLine)
  }

  case class callIndels (tumorBam : File, normalBam : File, rawIndels : File) extends SomaticIndelDetector with CommandLineGATKArgs {
    this.scatterCount = 1
    this.memoryLimit = 4
    this.jarFile = qscript.indelGenotyperJar

    val baseFilter = qscript.indelCallerBaseFilter

    this.input_file :+= new TaggedFile(tumorBam, "tumor")
    this.input_file :+= new TaggedFile(normalBam, "normal")

    this.verbose = rawIndels
    this.refseq = indelRefGenePath
    this.window_size = 400 // why?
    this.maxNumberOfReads = 8000 // why?
    this.filter = Seq(baseFilter)
  }

  /****************************************************************************
    * Classes (non-GATK programs)
    ****************************************************************************/


  case class filterIndelsNormal(originalIndelCalls : File, filteredIndelCalls : File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="indel call file to be filtered") var inCalls = originalIndelCalls
    @Output(doc="filtered indel calls file") var outCalls = filteredIndelCalls
    def commandLine = "perl " + indelFilterPath +
      " --calls " + inCalls +
      " --prefix N_" +
      " --max_cons_av_mm 1000 " +
      " --max_ref_av_mm 1000 " +
      " --max_cons_nqs_av_mm 1000 " +
      " --min_ref_nqs_av_qual 0 " +
      " --min_cons_nqs_av_qual 0 " +
      " --min_cons_count 0 " +
      " --min_readpos_median 10 " +
      " --min_readpos_mad 3 " +
      " --mode ANNOTATE " +
      " --output " + filteredIndelCalls

      this.analysisName = outCalls + ".filterIndelsNormal"
    this.jobName = this.analysisName
  }

  case class filterIndelsTumor(originalIndelCalls : File, filteredIndelCalls : File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="indel call file to be filtered") var inCalls = originalIndelCalls
    @Output(doc="filtered indel calls file") var outCalls = filteredIndelCalls
    def commandLine = "perl " + indelFilterPath +
      " --calls " + inCalls +
      " --prefix T_" +
      " --max_cons_av_mm 4 " +
      " --max_ref_av_mm 4 " +
      " --max_cons_nqs_av_mm 0.3 " +
      " --min_ref_nqs_av_qual 15 " +
      " --min_cons_nqs_av_qual 15 " +
      " --min_cons_count 0 " +
      " --min_readpos_median 10 " +
      " --min_readpos_mad 3 " +
      " --mode ANNOTATE " +
      " --output " + filteredIndelCalls

    this.analysisName = outCalls + ".filterIndelsTumor"
    this.jobName = this.analysisName
  }

  case class filterIndelsGermline(originalIndelCalls : File, filteredIndelCalls : File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="indel call file to be filtered") var inCalls = originalIndelCalls
    @Output(doc="filtered indel calls file") var outCalls = filteredIndelCalls

    def commandLine = "perl " + indelGermlineFilterPath +
      " --calls " + inCalls +
      " --filter " + indelGermlineDatabase +
      " --window 10 " +
      " --mode ANNOTATE " +
      " --pos_column 1" +
      " --output " + filteredIndelCalls

    this.analysisName = outCalls + ".filterIndelsGermline"
    this.jobName = this.analysisName
  }

  case class convertIndelCallsToMaflite(tumorName : String, normalName : String, originalIndelCalls : File, maflite : File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="indel call file to be filtered") var inCalls = originalIndelCalls
    @Output(doc="filtered indel calls in maflite format") var outCalls = maflite

    def commandLine = "perl " + indelToMaflitePath +
      " --build 37 " +
      " " + tumorName +
      " " + normalName +
      " " + originalIndelCalls +
      " " + maflite +
      " tumor_f,t_ref_count,t_alt_count "

    this.analysisName = outCalls + ".convertIndelsToMaflite"
    this.jobName = this.analysisName
  }
}
