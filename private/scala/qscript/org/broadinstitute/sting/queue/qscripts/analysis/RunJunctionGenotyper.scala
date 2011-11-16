package org.broadinstitute.sting.queue.qscripts.analysis

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.baq.BAQ
import org.broadinstitute.sting.queue.function.scattergather._
import org.broadinstitute.sting.queue.function._
import org.broadinstitute.sting.commandline.ArgumentSource
import java.io.File
import org.broadinstitute.sting.utils.text.XReadLines
import java.io.PrintStream
import scala.collection.JavaConversions._
import org.broadinstitute.sting.queue.extensions.gatk._



class RunJunctionGenotyper extends QScript {

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "r", doc="refSeq", required=false)
  var refSeqFile: File = new File("/humgen/gsa-hpprojects/GATK/data/refGene_b37.sorted.txt")

  @Argument(shortName = "bam", doc="BAM", required=true)
  val bam: File = null;

  @Argument(shortName = "o", doc = "output",required=true)
  val output: File = null;

  @Argument(shortName= "L", doc="interval list",required=false)
  val intervals: File = new File("/humgen/gsa-hphome1/chartl/projects/oneoffs/targets/all_genes.b37.interval_list")

   trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    logging_level = "INFO";
    memoryLimit = 6;
    this.reference_sequence = referenceFile
    this.input_file :+= bam
  }

  def bai(bam: File) = new File(bam + ".bai")

  // 1.) Unified Genotyper Base
  class GenotyperBase (t: File) extends IntronLossGenotyperV2 with UNIVERSAL_GATK_ARGS {
    this.gatherClass = { case argumentSource : ArgumentSource if ( argumentSource.field.getName=="out" ) => classOf [MyGatherFunction] }
    this.intervals :+= t
    this.scatterCount = 20
    this.refSeq = new TaggedFile(refSeqFile,"REFSEQ")
  }

  class InsertSizeDistributionBase extends InsertSizeDistribution {
    this.intervalsString ++= (new Range(1,23,1)).map( _.toString() )
    this.reference_sequence = referenceFile
  }

  def script = {
    if ( bam.getName().endsWith(".bam") ) {
      val insertSize : InsertSizeDistribution = new InsertSizeDistributionBase()
      insertSize.out = swapExt(bam,".bam",".isd.table")
      insertSize.input_file :+= bam
      insertSize.memoryLimit = Some(4)
      add(insertSize)
      val genotype : IntronLossGenotyperV2 = new GenotyperBase(intervals)
      genotype.insertHistogram :+= insertSize.out
      genotype.out = output
      add(genotype)
    } else {
      val genotype : IntronLossGenotyperV2 = new GenotyperBase(intervals)
      genotype.out = output
      for ( s <- asScalaIterator(new XReadLines(bam)) ) {
        var insertSize : InsertSizeDistribution = new InsertSizeDistributionBase()
        insertSize.out = swapExt(s,".bam",".isd.table")
        insertSize.reference_sequence = referenceFile
        insertSize.input_file :+= new File(s)
        insertSize.memoryLimit = Some(4)
        genotype.insertHistogram :+= insertSize.out
        add(insertSize)
      }
      add(genotype)
    }
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
