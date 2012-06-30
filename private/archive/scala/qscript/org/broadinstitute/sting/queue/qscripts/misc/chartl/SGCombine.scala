package org.broadinstitute.sting.queue.qscripts.misc.chartl

import org.broadinstitute.sting.queue.QScript
import java.io.PrintStream
import org.broadinstitute.sting.utils.text.XReadLines
import collection.JavaConversions._
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.GenotypeMergeType
import org.broadinstitute.sting.queue.extensions.gatk.{TaggedFile, SelectVariants, CombineVariants}

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 2/16/12
 * Time: 7:54 PM
 * To change this template use File | Settings | File Templates.
 */

class SGCombine extends QScript {
  @Input(doc="The input VCFs you want to merge without subset",shortName="V",fullName="V",required=true)
  var vcfNoSubsetInput : List[File] = Nil
  @Input(doc="The input VCFs you want to merge with subset",shortName="VS",fullName="VS",required=true)
  var vcfWithSubsetInput : List[File] = Nil
  @Input(doc="The samples you want to subset to",shortName="sf",fullName="sf",required=true)
  var samplesFile : File = _
  @Output(doc="The output VCF to write to",shortName="o",fullName="o",required=true)
  var outFile : File = _
  @Argument(doc="The number of variants per running job",fullName="c",shortName="c",required=false)
  var maxChunkSize : Int = 50000;
  @Input(doc="The interval list to run over",required=true,shortName="L",fullName="L")
  var intervals : File = _
  @Input(fullName="preserveChromosomes",doc="Restrict chunks to one chromosome (smaller chunk at end of chromosome)",required=false)
  var preserve : Boolean = false
  @Input(fullName="reference",doc="The reference file",required=false)
  var ref : File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  @Input(doc="Genotype prioritization (by filename, strip off .vcf)",shortName="P",fullName="P",required=true)
  var priority : String = _


  val tmpdir : File = System.getProperty("java.io.tmpdir")

  def script = {
    var vcfNoSubset = vcfNoSubsetInput.map( u => new TaggedFile(u,u.getName.stripSuffix(".gz").stripSuffix(".vcf")))
    var vcfWithSubset = vcfWithSubsetInput.map(u => new TaggedFile(u,u.getName.stripSuffix(".gz").stripSuffix(".vcf")))
    var chunkNum = 1
    var numLinesInChunk = 0
    var chromosome : String = asScalaIterator(new XReadLines(intervals)).next().split(":")(0)
    var chunkFile : File = new File(tmpdir,"ChunkVCF.chunk%d.intervals.list".format(chunkNum))
    var chunkWriter = new PrintStream(chunkFile)
    var combinedChunks : List[File] = Nil
    asScalaIterator(new XReadLines(intervals)).foreach( int => {
      // check new chromosome or full chunk
      if ( ( preserve && ! int.split(":")(0).equals(chromosome) ) || numLinesInChunk > maxChunkSize ) {
        chunkWriter.close()
        // first, subset the VCFs that need be subset over this interval
        val chunkSelect : List[SelectVariants] = vcfWithSubset.map( u => {
          val v : SelectVariants = new SelectVariants
          v.reference_sequence = ref
          v.memoryLimit = 2
          v.sample_file :+= samplesFile
          v.intervals :+= chunkFile
          v.variant = u
          v.out = swapExt(tmpdir,u,".vcf.gz",".chunk%d.selected.vcf".format(chunkNum))
          v
        })
        addAll(chunkSelect)
        // second, combine the VCFs over the region
        val chunkCombine : CombineVariants = new CombineVariants
        chunkCombine.reference_sequence = ref
        chunkCombine.memoryLimit = 4
        chunkCombine.intervals :+= chunkFile
        chunkCombine.variant ++= chunkSelect.map(u => new TaggedFile(u.out,u.out.getName.stripSuffix(".gz").stripSuffix(".vcf")))
        chunkCombine.variant ++= vcfNoSubset.map(u => u.asInstanceOf[File] )
        chunkCombine.out = swapExt(tmpdir,outFile,".vcf",".chunk%d.combined.vcf".format(chunkNum))
        chunkCombine.genotypeMergeOptions = GenotypeMergeType.PRIORITIZE
        chunkCombine.priority = parsePriority(priority,chunkCombine.variant)
        add(chunkCombine)
        combinedChunks :+= chunkCombine.out
        chunkNum += 1
        numLinesInChunk = 0
        chromosome = int.split(":")(0)
        chunkFile = new File(tmpdir,"ChunkVCF.chunk%d.intervals.list".format(chunkNum))
        chunkWriter = new PrintStream(chunkFile)
      }
      chunkWriter.printf("%s%n",int)
      numLinesInChunk += 1
    })
    // last chunk
    if ( numLinesInChunk > 0 ) {
      // some work to do
      val chunkSelect : List[SelectVariants] = vcfWithSubset.map( u => {
        val v : SelectVariants = new SelectVariants
        v.reference_sequence = ref
        v.memoryLimit = 2
        v.sample_file :+= samplesFile
        v.intervals :+= chunkFile
        v.variant = u
        v.out = swapExt(tmpdir,u,".vcf.gz",".chunk%d.vcf".format(chunkNum))
        v
      })
      chunkWriter.close()
      addAll(chunkSelect)
      val chunkCombine : CombineVariants = new CombineVariants
      chunkCombine.reference_sequence = ref
      chunkCombine.memoryLimit = 4
      chunkCombine.intervals :+= chunkFile
      chunkCombine.variant ++= chunkSelect.map(u => u.out)
      chunkCombine.variant ++= vcfNoSubset
      chunkCombine.out = swapExt(tmpdir,outFile,".vcf.gz",".chunk%d.combined.vcf".format(chunkNum))
      add(chunkCombine)
      combinedChunks :+= chunkCombine.out
    }

    var gather : MyVCFGather = new MyVCFGather
    gather.vcfs ++= combinedChunks
    gather.outVCF = outFile
    add(gather)
  }

  def parsePriority(baseP:String, binds : Seq[File] ) : String = {
    logger.debug(baseP)
    logger.debug(binds.map(u => u.getName).reduceLeft(_+","+_))
    logger.debug(binds.map(t => t.getName.startsWith(baseP.split(",").head)).map(v => v.toString).reduceLeft(_+","+_))
    baseP.split(",").map(u =>  binds.filter( t => t.getName.startsWith(u)).head ).map(v => v.getName.stripSuffix(".gz").stripSuffix(".vcf")).reduceLeft(_ + "," + _)
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

