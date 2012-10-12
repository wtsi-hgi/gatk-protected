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

  @Input(doc="The path to the indel caller refgene file", fullName="indel_refgene", shortName="ir", required=true)
  var indelRefGenePath: File = _


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
    var tumorBams: Seq[File] = Seq(bams(0))
    var normalBams: Seq[File] = Seq(bams(1))

    assert(bams.length <= 5, "Current implementation is limited to 5 bams. See source code for details on the n-way out indel cleaning procedure")

    print("Indicies:" + tumorBamIndicies + "\n")
    print("Tumor:" + tumorBams + "\n")
    print("Normal:" + normalBams + "\n")

    for (tumorBam <- tumorBams) {
      for (normalBam <- normalBams) {
        val outPrefix = tumorBams(0).getName + "-vs-" + normalBams(0).getName // TODO: use CMI ids here
        cancer(tumorBam, normalBam, 0.01f, outPrefix)
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
    this.intervals = qscript.intervals
  }

  def cancer (tumorBam : File, normalBam : File, tumorFractionContamination : Float, outPrefix : String) {
    val rawMutations = outPrefix + ".call_stats.txt"
    val rawVcf = outPrefix + ".vcf"
    val rawCoverage = outPrefix + ".wig.txt"

    add(mutect(tumorBam, normalBam, tumorFractionContamination, rawMutations, rawVcf, rawCoverage))
    indels(tumorBam, normalBam, outPrefix)

  }

  def indels (tumorBam : File, normalBam : File, outPrefix : String) {
    val rawIndels = outPrefix + ".raw.indels.txt"
    val filteredIndels = outPrefix + ".filtered.indels.txt"

    add(callIndels(tumorBam, normalBam, rawIndels))
 //   add(filterIndelsStep1(rawIndels, filteredIndels))
  }

  case class mutect (tumorBam : File, normalBam : File, tumorFractionContamination : Float, outMutations : File,  outVcf : File, outCoverage : File) extends MuTect with CommandLineGATKArgs {
    this.scatterCount = 1
    this.memoryLimit = 4
    this.jarFile = qscript.mutectJar

    this.dbsnp = qscript.dbSNP
    this.cosmic = qscript.cosmic
    this.normal_panel = qscript.pon

    this.only_passing_calls = true
    this.enable_extended_output = true
    this.downsample_to_coverage = 1000 // TODO: how deep should this be?
    this.fraction_contamination = Some(tumorFractionContamination)

    this.input_file :+= new TaggedFile(tumorBam, "tumor")
    this.input_file :+= new TaggedFile(normalBam, "normal")

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

    val baseFilter = qscript.indelCallerBaseFilter;

    this.input_file :+= new TaggedFile(tumorBam, "tumor")
    this.input_file :+= new TaggedFile(normalBam, "normal")

    this.verbose = rawIndels
    this.refseq = indelRefGenePath
    this.window_size = 400 // why?
    this.maxNumberOfReads = 8000 // why?
    this.filter = Seq(baseFilter)
  }

  case class filterIndelsStep1(originalIndelCalls : File, filteredIndelCalls : File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="indel call file to be filtered") var inCalls = originalIndelCalls
    @Output(doc="filtered indel calls file") var outCalls = filteredIndelCalls
    def commandLine = indelFilterPath +
      " --calls " + inCalls + " --prefix ?? --params <params.file> " +
      " --max_cons_av_mm <max.cons.av.mm> --max_ref_av_mm <max.ref.av.mm> "
      " --max_cons_nqs_av_mm <max.cons.nqs.av.mm> --min_ref_nqs_av_qual <min.ref.nqs.av.qual> " +
      " --min_cons_nqs_av_qual <min.cons.nqs.av.qual> --min_cons_count <min.cons.count> " +
      " --min_readpos_median <min.readpos.median> --min_readpos_mad <min.readpos.mad> " +
      " --mode ANNOTATE --output " + filteredIndelCalls

    this.analysisName = outCalls + ".filterIndelsStep1"
    this.jobName = this.analysisName
  }

  /****************************************************************************
    * Classes (non-GATK programs)
    ****************************************************************************/


}

