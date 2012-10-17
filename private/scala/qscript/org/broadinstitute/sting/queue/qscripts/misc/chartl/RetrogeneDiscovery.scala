package org.broadinstitute.sting.queue.qscripts.analysis

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.baq.BAQ
import org.broadinstitute.sting.queue.function.scattergather._
import org.broadinstitute.sting.queue.function._
import java.io.File
import org.broadinstitute.sting.utils.text.XReadLines
import java.io.PrintStream
import scala.collection.JavaConversions._
import org.broadinstitute.sting.utils.interval.IntervalSetRule
import org.broadinstitute.sting.commandline.{Output, Input, ArgumentSource}
import org.broadinstitute.sting.queue.extensions.gatk.{InsertSizeDistribution, ExonJunctionHypothesisGenerator, CommandLineGATK, TaggedFile, ExonJunctionGenotyper }

class RetrogeneDiscovery extends QScript {

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "r", doc="refSeq", required=false)
  var rsFile: File = new File("/humgen/gsa-hpprojects/GATK/data/refGene_b37.sorted.txt")

  @Argument(shortName = "hb", doc="Bam list for hypothesis generation", required=true)
  var bam: File = null;

  @Argument(shortName = "eb", doc="Bam list for hypothesis evaluation",required=false)
  var evalBam : File = null;

  @Argument(shortName = "o", doc = "output",required=true)
  var output: File = null;

  @Argument(shortName= "gl", doc="gene list",required=true)
  var geneIntervalsToRunOver: File = new File("/humgen/gsa-hphome1/chartl/projects/oneoffs/targets/all_genes.b37.interval_list")

  @Argument(shortName = "L", required = false, doc="foo")
  var captureIntervalList : File = _

  @Argument(shortName = "H", required = false, doc = "If you've preprocessed insert size interval lists, you can supply the list file here")
  var insertDistributionFileList : File = _

  @Argument(shortName = "intermediate", required=true, doc="The writing directory for intermediate files")
  var intermediateDir : String = _

  @Argument(fullName="histogramOnly",required=false,doc="Make histograms only")
  var histogramOnly : Boolean = false

  var refSeqFile : TaggedFile = _

   trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.memoryLimit = 2;
    this.reference_sequence = referenceFile
  }

  def bai(bam: File) = new File(bam + ".bai")

  def hypGen(bamList:File, ival: File, outFile: File) : ExonJunctionHypothesisGenerator = {
    var generator = new ExonJunctionHypothesisGenerator with UNIVERSAL_GATK_ARGS
    generator.refSeq = refSeqFile;
    generator.out = outFile
    generator.intervals :+= ival
    generator.input_file :+= bamList
    if ( captureIntervalList != null ){
        generator.intervals :+= captureIntervalList
        generator.interval_set_rule = IntervalSetRule.INTERSECTION
      }
    generator.memoryLimit = 4
    generator
  }

  def isdTable(bam: File) : InsertSizeDistribution = {
    var isd = new InsertSizeDistribution
    isd.reference_sequence = referenceFile
    isd.intervalsString ++= (new Range(1,23,1)).map(_.toString())
    // note -- probably good enough just to do chr1 and chr2
    isd.memoryLimit = Some(2)
    isd.out = swapExt(intermediateDir,bam,".bam",".isd.table")
    isd.input_file :+= bam
    isd.maxInsertSize = 2000
    isd
  }

  def hypEval(ival: File, outFile: File, tableList: List[File], hypothesis: File) : ExonJunctionGenotyper = {
    var gt = new ExonJunctionGenotyper with UNIVERSAL_GATK_ARGS
    gt.hypothesis = hypothesis
    gt.insertHistogram ++= tableList
    gt.intervals :+= ival
    gt.refSeq = refSeqFile
    gt.out = outFile
    gt.input_file :+= evalBam
    if ( captureIntervalList != null ) {
        gt.intervals :+= captureIntervalList
        gt.interval_set_rule = IntervalSetRule.INTERSECTION
      }
    gt
  }

  def script = {
    if ( evalBam == null ) {
      evalBam = bam
    }

    refSeqFile = new TaggedFile(rsFile,"REFSEQ")
    var histograms : List[File] = Nil
    if ( ! bam.getName().endsWith(".list") ) {
      println("Script not implemented for single bam files. Please supply a bam list, ending with '.list'")
    } else if ( insertDistributionFileList == null ){
      for ( s <- asScalaIterator(new XReadLines(evalBam)) ) {
        var insertSize : InsertSizeDistribution = isdTable(s)
        histograms :+= insertSize.out
        add(insertSize)
      }
    } else {
      // scripter's note: if you supply an ISD histogram file on the cmd line, making histograms will be omitted
      histograms :+= insertDistributionFileList
    }

    if ( ! geneIntervalsToRunOver.getName.endsWith(".list")) {
      println("Script not implemented for single intervals. Please supply an interval list to break up over, ending with '.list'")
    } else {
      var intervalCounter = 1
      var vcfIntermediates :  List[File] = Nil
      for ( s <- (asScalaIterator(new XReadLines(geneIntervalsToRunOver))).grouped(10) ) {
        var intermediateHypothesis = new File(intermediateDir,"RJG_Hyp_%d.tbl".format(intervalCounter))
        var intermediateVCF = new File(intermediateDir,"RJG_Hyp_%d.vcf".format(intervalCounter))
        var intermediateInterval = new File(intermediateDir,"RJG_Hyp_%d.intervals.list".format(intervalCounter))
        printIntervals(s.toList,intermediateInterval)
        var hyp = hypGen(bam,intermediateInterval,intermediateHypothesis)
        var typ = hypEval(intermediateInterval,intermediateVCF,histograms,new TaggedFile(intermediateHypothesis,"TABLE"))
        vcfIntermediates :+= intermediateVCF
        if ( ! histogramOnly )
         add(hyp,typ)
        intervalCounter = intervalCounter + 1
      }
      // want to gather the intermediates through an in process function
      var gather = new MyVCFGather
      gather.vcfs = vcfIntermediates
      gather.outVCF = output
      if ( ! histogramOnly )
       add(gather)
    }
  }

  def printIntervals(ival : List[String], iFile : File) {
    var stream : PrintStream = new PrintStream(iFile)
    ival.foreach(u => stream.printf("%s%n",u))
    stream.close()
  }
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

  class MyGatherFunction extends GatherFunction with InProcessFunction {

    def run : Unit = {

      var stream : PrintStream = new PrintStream(originalOutput)
      var first : XReadLines = new XReadLines(gatherParts(0))
      var line : String = first.next()
      while ( line.startsWith("#") ) {
        stream.printf("%s%n",line)
        line = first.next()
      }
      first.close()
      for ( f <- gatherParts ) {
        var xrl : XReadLines = new XReadLines(f)
        for ( s <- asScalaIterator(xrl) ) {
          if ( ! s.startsWith("#") ) {
            stream.printf("%s%n",s)
          }
        }
        xrl.close()
      }

      // originalOutput
      // gatherParts
    }
  }
