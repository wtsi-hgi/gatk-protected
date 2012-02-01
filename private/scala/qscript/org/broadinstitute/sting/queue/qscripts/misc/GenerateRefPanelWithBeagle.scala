package org.broadinstitute.sting.queue.qscripts.misc

import org.broadinstitute.sting.queue.QScript
import scala.collection.JavaConversions._
import org.broadinstitute.sting.utils.text.XReadLines
import org.broadinstitute.sting.commandline.Argument
import java.io.PrintStream
import org.broadinstitute.sting.queue.extensions.gatk._
import net.sf.picard.reference.FastaSequenceIndex
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource
import collection.mutable.HashSet

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/20/12
 * Time: 11:36 AM
 * To change this template use File | Settings | File Templates.
 */

class GenerateRefPanelWithBeagle extends QScript {
  @Input(shortName="LP",fullName="lowPass",doc="VCF from low-pass calls and genotypes, containing likelihoods",required=true)
  var lowPassCalls : File = _
  @Input(shortName="EX",fullName="exome",doc="VCF from exome calls and genotypes, containing likelihoods",required=true)
  var exomeCalls : File = _
  @Input(shortName="CHIP",fullName="genotypeChip",doc="VCF from genotype chip genotypes, not necessarily containing likelihoods",required=true)
  var chipCalls : File = _
  @Input(shortName="CI",fullName="chunkIntervals",doc="Interval list containing the chunking strategy for running this analysis",required=true)
  var chunks : File = _
  @Input(shortName="R",fullName="referenceFile",doc="The reference fasta. Should be obvious.",required=false)
  var ref : File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  @Input(shortName="B",fullName="beagleJar",doc="Path to the beagle jar",required=false)
  var beagleJar = new File("/humgen/gsa-hpprojects/software/beagle/beagle.jar")
  @Input(shortName="PC",fullName="phasedComparison",doc="A file phased with beagle for comparison of phasing and genotypes. Usually lowpassed phased on its own.",required=true)
  var phasedComp : File = _
  @Input(shortName="X",fullName="excludeIdx",doc="Exclude these job indeces (useful for beagle failing when no variants present)",required=false)
  var exclude : File = _

  val BEAGLE_MEM_IN_GB : Int = 8
  val TMPDIR : String = System.getProperty("java.io.tmpdir");


  class RefineGenotypesWithBeagle(@Input beagleInput: File) extends CommandLineFunction {
    this.memoryLimit = BEAGLE_MEM_IN_GB

    // Note: These get set
    @Output val beaglePhasedFile: File = new File(beagleInput.getAbsolutePath +".phased.gz")
    @Output val beagleLikelihoods: File = new File(beagleInput.getAbsolutePath +".gprobs.gz")
    @Output val beagleRSquared: File = new File(beagleInput.getAbsolutePath +".r2")

    def commandLine = "java -Djava.io.tmpdir=%s -Xmx%dg -jar %s out=foo like=%s niterations=50 omitprefix=true".format(
      TMPDIR, BEAGLE_MEM_IN_GB, beagleJar,beagleInput.getAbsolutePath)
  }

  class Gunzip(@Input zippedFile : File ) extends CommandLineFunction {
    @Output val outputFile : File = new File(zippedFile.getAbsolutePath.replace(".gz",""))

    def commandLine = "gunzip -c %s > %s".format(zippedFile.getAbsolutePath,outputFile.getAbsolutePath)
  }

  def script = {

    // identify what to exclude
    var exIdx : HashSet[Int] = new HashSet[Int]
    if ( exclude != null ) {
      exIdx.addAll(asJavaCollection(asScalaIterable(new XReadLines(exclude)).map(u => u.toInt)))
    }
    // for each chunk
    var idx : Int = 0
    var panelChunks : List[File] = Nil
    asScalaIterator(new XReadLines(chunks)).foreach( chunk => {
      if ( ! exIdx.contains(idx) ) {
        // calculate the flanks (150kb)
        val chr = chunk.split(":")(0)
        val start = Integer.parseInt(chunk.split(":")(1).split("-")(0))
        val stop = Integer.parseInt(chunk.split(":")(1).split("-")(1))
        var flanks : List[String] = Nil
        var leftHalfFlank : String = null
        if ( start > 150000 ) {
          flanks :+= "%s:%d-%d".format(chr,start-150000,start)
          leftHalfFlank = "%s:%d-%d".format(chr,start-75000,start)
        }
        var rightHalfFlank : String = null
        val dif : Int = (new ReferenceDataSource(ref)).getReference.getSequenceDictionary.getSequence(chr).getSequenceLength - stop
        if ( dif > 0 ) {
          val ept = scala.math.min(dif,150000)
          flanks :+= "%s:%d-%d".format(chr,stop,stop+ept)
          rightHalfFlank = "%s:%d-%d".format(chr,stop,stop+ept/2)
        }
        // 1: multiply together the likelihoods
        val mLik = new MultiplyLikelihoods
        mLik.reference_sequence = ref
        mLik.Variants :+= lowPassCalls
        mLik.Variants :+= exomeCalls
        mLik.Variants :+= chipCalls
        mLik.intervalsString :+= chunk
        mLik.intervalsString ++= flanks
        mLik.out = new File(TMPDIR,"Mult_Likelihoods.chunk%d.vcf".format(idx))
        mLik.memoryLimit = 2
        add(mLik)

        // 2: create the input for beagle
        val beagOut = new ProduceBeagleInput
        beagOut.reference_sequence = ref
        beagOut.variant = mLik.out
        beagOut.out = new File(TMPDIR,"Mult_Likelihoods.chunk%d.beagle".format(idx))
        beagOut.intervalsString :+= chunk
        beagOut.intervalsString ++= flanks
        beagOut.memoryLimit = 2
        add(beagOut)

        // 3: refine the genotypes (and phase)
        val beagRefine = new RefineGenotypesWithBeagle(beagOut.out)
        add(beagRefine)

        // 3a: gunzip this stuff
        val phase = new Gunzip(beagRefine.beaglePhasedFile)
        add(phase)
        val like = new Gunzip(beagRefine.beagleLikelihoods)
        add(like)

        // 4: convert beagle output back to VCF
        val beag2vcf = new BeagleOutputToVCF
        beag2vcf.reference_sequence = ref
        beag2vcf.memoryLimit = Some(4)
        beag2vcf.out = new File(TMPDIR,"Ref_Panel.chunk%d.vcf".format(idx))
        beag2vcf.keep_monomorphic = true
        beag2vcf.beaglePhased = new TaggedFile(phase.outputFile,"ph,BEAGLE")
        beag2vcf.beagleProbs = new TaggedFile(like.outputFile,"pr,BEAGLE")
        beag2vcf.beagleR2 = new TaggedFile(beagRefine.beagleRSquared,"r2,BEAGLE")
        beag2vcf.intervalsString :+= chunk
        if ( leftHalfFlank != null ) {
          beag2vcf.intervalsString :+= leftHalfFlank
        }
        if ( rightHalfFlank != null ) {
          beag2vcf.intervalsString :+= rightHalfFlank
        }
        beag2vcf.variant = mLik.out
        panelChunks :+= beag2vcf.out
        beag2vcf.memoryLimit = 2
        add(beag2vcf)
      }
      idx = 1 + idx;
    })

    val gather = new MyVCFGather
    gather.vcfs ++= panelChunks
    gather.outVCF = new File("Reference_Panel_Final.vcf")
    add(gather)

    // now we want to compare this with a pre-phased file to ensure no truly gross errors
    val gtEval = new VariantEval
    gtEval.reference_sequence = ref
    gtEval.eval :+= new TaggedFile(gather.outVCF,"Consensus,VCF")
    gtEval.comp :+= new TaggedFile(phasedComp,"PhasedComp,VCF")
    gtEval.EV :+= "GenotypeConcordance"
    gtEval.EV :+= "GenotypePhasingEvaluator"
    gtEval.out = new File("Reference_Panel_Eval.gatk")
    gtEval.intervals :+= chunks
    gtEval.memoryLimit = 4
    add(gtEval)
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