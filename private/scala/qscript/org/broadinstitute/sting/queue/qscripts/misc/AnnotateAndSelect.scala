package org.broadinstitute.sting.queue.qscripts.misc

import org.broadinstitute.sting.queue.QScript
import scala.collection.JavaConversions._
import org.broadinstitute.sting.utils.text.XReadLines
import org.broadinstitute.sting.queue.extensions.gatk.{SelectVariants, VariantAnnotator}
import java.io.PrintStream

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/24/12
 * Time: 2:47 PM
 * To change this template use File | Settings | File Templates.
 */

class AnnotateAndSelect extends QScript {

  @Input(fullName="omni",shortName="omni",doc="input omni vcf",required=true)
  var omni : File = _

  @Input(fullName="chunk",shortName = "chunk",doc="chunk interval list",required=true)
  var chunk : File = _

  @Output(fullName="mono",shortName = "mono",doc="output mono file",required=true)
  var monoVCF : File = _

  @Output(fullName="poly",shortName="poly",doc="output poly file",required=true)
  var polyVCF : File = _

  val ref : File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  val TMPDIR : String = System.getProperty("java.io.tmpdir");

  def script {

    var monoOut : List[File] = Nil
    var polyOut : List[File] = Nil
    var idx : Int = 1
    asScalaIterator(new XReadLines(chunk)).foreach( u => {
      var annotate = new VariantAnnotator
      annotate.reference_sequence = ref
      annotate.intervalsString :+= u
      annotate.variant = omni
      annotate.A :+= "ChromosomeCounts"
      annotate.out = new File(TMPDIR,"omni_annotate_%d.vcf".format(idx))
      annotate.memoryLimit = 2
      add(annotate)
      var selPoly = new SelectVariants
      selPoly.reference_sequence = ref
      selPoly.intervalsString :+= u
      selPoly.variant = annotate.out
      selPoly.select_expressions :+= "AC>0"
      selPoly.out = new File(TMPDIR,"omni_select_poly_%d.vcf".format(idx))
      selPoly.memoryLimit = 2
      add(selPoly)
      var selMono = new SelectVariants
      selMono.reference_sequence = ref
      selMono.variant = annotate.out
      selMono.intervalsString :+= u
      selMono.excludeIntervals :+= selPoly.out
      selMono.out = new File(TMPDIR,"omni_select_mono_%d.vcf".format(idx))
      selMono.memoryLimit = 2
      add(selMono)
      monoOut :+= selMono.out
      polyOut :+= selPoly.out
      idx = idx + 1
    })

    val monoGather = new MyVCFGather
    monoGather.vcfs = monoOut
    monoGather.outVCF = monoVCF
    add(monoGather)

    val polyGather = new MyVCFGather
    polyGather.vcfs = polyOut
    polyGather.outVCF = polyVCF
    add(polyGather)
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