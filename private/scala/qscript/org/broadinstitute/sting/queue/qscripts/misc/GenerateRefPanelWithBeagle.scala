package org.broadinstitute.sting.queue.qscripts.misc

import org.broadinstitute.sting.queue.QScript
import scala.collection.JavaConversions._
import org.broadinstitute.sting.utils.text.XReadLines
import org.broadinstitute.sting.commandline.Argument
import org.broadinstitute.sting.queue.qscripts.misc.GenerateRefPanelWithBeagle.RefineGenotypesWithBeagle
import org.broadinstitute.sting.queue.extensions.gatk.{BeagleOutputToVCF, ProduceBeagleInput, MultiplyLikelihoods}
import java.io.PrintStream
import org.broadinstitute.sting.queue.qscripts.misc.GenerateRefPanelWithBeagle.MyVCFGather

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/20/12
 * Time: 11:36 AM
 * To change this template use File | Settings | File Templates.
 */

class GenerateRefPanelWithBeagle extends QScript {
  @Input(shortName="LP",fullName="lowPass",doc="VCF from low-pass calls and genotypes, containing likelihoods",required=true)
  val lowPassCalls : File = _
  @Input(shortName="EX",fullName="exome",doc="VCF from exome calls and genotypes, containing likelihoods",required=true)
  val exomeCalls : File = _
  @Input(shortName="CHIP",fullName="genotypeChip",doc="VCF from genotype chip genotypes, not necessarily containing likelihoods",required=true)
  val chipCalls : File = _
  @Input(shortName="L",fullName="chunkIntervals",doc="Interval list containing the chunking strategy for running this analysis",required=true)
  val chunks : File = _
  @Input(shortName="R",fullName="referenceFile",doc="The reference fasta. Should be obvious.",required=false)
  val ref : File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  @Input(shortName="B",fullName="beagleJar",doc="Path to the beagle jar",required=false)
  val beagleJar = new File("/humgen/gsa-hpprojects/software/beagle/beagle.jar")

  val BEAGLE_MEM_IN_GB : Int = 8
  val TMPDIR : String = System.getProperty("java.io.tmpdir");


  class RefineGenotypesWithBeagle(@Input beagleInput: File, @Argument outputBase : String) extends CommandLineFunction {
    this.memoryLimit = BEAGLE_MEM_IN_GB

    // Note: These get set
    @Output val beaglePhasedFile: File = new File(outputBase +".phased.gz")
    @Output val beagleLikelihoods: File = new File(outputBase +".gprobs.gz")
    @Output val beagleRSquared: File = new File(outputBase +".r2")

    def commandLine = "java -Djava.io.tmpdir=%s -Xmx%dg -jar %s out=%s like=%s niterations=50".format(
      TMPDIR, BEAGLE_MEM_IN_GB, beagleJar,outputBase,beagleInput.getAbsolutePath)
  }

  def script = {

    // for each chunk
    var idx : Int = 0
    var panelChunks : List[File] = Nil
    asScalaIterator(new XReadLines(chunks)).foreach( chunk => {
      // 1: multiply together the likelihoods
      val mLik = new MultiplyLikelihoods
      mLik.reference_sequence = ref
      mLik.Variants :+= lowPassCalls
      mLik.Variants :+= exomeCalls
      mLik.Variants :+= chipCalls
      mLik.intervalsString :+= chunk
      mLik.out = new File("Mult_Likelihoods.chunk%d.vcf".format(idx))
      add(mLik)

      // 2: create the input for beagle
      val beagOut = new ProduceBeagleInput
      beagOut.reference_sequence = ref
      beagOut.variant = mLik.out
      beagOut.out = new File("Mult_Likelihoods.chunk%d.beagle".format(idx))
      beagOut.intervalsString :+= chunk
      beagOut.memoryLimit = Some(BEAGLE_MEM_IN_GB)
      add(beagOut)

      // 3: refine the genotypes (and phase)
      val beagRefine = new RefineGenotypesWithBeagle(beagOut.out,"Beag_Raw_Out.chunk%d".format(idx))
      add(beagRefine)

      // 4: convert beagle output back to VCF
      val beag2vcf = new BeagleOutputToVCF
      beag2vcf.reference_sequence = ref
      beag2vcf.memoryLimit = Some(4)
      beag2vcf.out = new File("Ref_Panel.chunk%d.vcf".format(idx))
      beag2vcf.keep_monomorphic = true
      beag2vcf.beaglePhased = beagRefine.beaglePhasedFile
      beag2vcf.beagleProbs = beagRefine.beagleLikelihoods
      beag2vcf.beagleR2 = beagRefine.beagleRSquared
      beag2vcf.intervalsString :+= chunk
      add(beag2vcf)
      panelChunks :+= beag2vcf.out
      idx.+
    })

    val gather = new MyVCFGather
    gather.vcfs ++= panelChunks
    gather.outVCF = new File("Reference_Panel_Final.vcf")
  }

    class MyVCFGather extends InProcessFunction {
    @Input(doc="VCFs to be merged") var vcfs: List[File] = Nil
    @Output(doc="The final VCF to write to") var outVCF : File = _

    def run : Unit = {
      var stream : PrintStream = new PrintStream(outVCF)
      var first : XReadLines = new XReadLines(vcfs(0))
      var line : String = first.next()
      while ( line != null && line.startsWith("#") ) {
        stream.printf("%s%n",line)
        if ( ! first.hasNext ) {
          line = null
        } else {
          line = first.next()
        }
      }
      first.close()
      vcfs.map(
        u => asScalaIterator[String](new XReadLines(u)) ).foreach(
        x => x.filter( z => ! z.startsWith("#") ).foreach(
         s => stream.printf("%s%n",s)))
      stream.close()
    }
  }

}