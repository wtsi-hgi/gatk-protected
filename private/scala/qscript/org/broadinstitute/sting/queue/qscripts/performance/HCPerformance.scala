/*
 * Copyright (c) 2012 The Broad Institute
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.queue.qscripts.performance

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.PathUtils
import org.broadinstitute.sting.commandline.ClassType
import org.broadinstitute.sting.gatk.downsampling.DownsampleType

class HCPerformance extends QScript {
  @Argument(shortName = "resources", doc="resources", required=true)
  val resourcesDir: String = ""

  @Argument(shortName = "longRun", doc="Run for 5x the normal time", required=false)
  val longRun : Boolean = false

  @Argument(shortName = "iterations", doc="it", required=false)
  val iterations: Int = 1

  val b37_FILENAME = "human_g1k_v37.fasta"

  def makeResource(x: String): File = new File("%s/%s".format(resourcesDir, x))

  class Data(val name: String,
                   val bamList: File,
                   val intervalsFile: File,
                   val intervalsString: List[String]) {
    def addIntervals(gatkCmd : CommandLineGATK): CommandLineGATK = {
      if ( intervalsFile != null )
        gatkCmd.intervals :+= intervalsFile
      gatkCmd.intervalsString = intervalsString
      gatkCmd.interval_set_rule = org.broadinstitute.sting.utils.interval.IntervalSetRule.INTERSECTION
      gatkCmd
    }
  }

  val assessments: List[(String, () => CommandLineGATK)] = List(
    ("CountReads", makeCountReads),
    ("CountLoci", makeCountLoci),
    ("CountReads.ART", makeCountReadsART),
    ("HC.justProfile.noDS", makeHC(ds = false,  justProfile = true)),
    ("HC.justProfile",      makeHC(ds = true,   justProfile = true)),
    ("HC.noDS",             makeHC(ds = false,  justProfile = false)),
    ("HC",                  makeHC(ds = true,   justProfile = false))
  )

  def makeCountReadsART(): CommandLineGATK = {
    new CountReadsInActiveRegions() with UNIVERSAL_GATK_ARGS
  }

  def makeCountLoci(): CommandLineGATK = {
    new CountLoci() with UNIVERSAL_GATK_ARGS
  }

  def makeCountReads(): CommandLineGATK = {
    new CountReads() with UNIVERSAL_GATK_ARGS
  }


  def makeHC(ds: Boolean, justProfile: Boolean): () => CommandLineGATK = {
    def maker(): CommandLineGATK = {
      val HC = new HaplotypeCaller() with UNIVERSAL_GATK_ARGS
      HC.out = new File("/dev/null")
      HC.justDetermineActiveRegions = justProfile
      if ( ! ds )
        HC.downsampling_type = org.broadinstitute.sting.gatk.downsampling.DownsampleType.NONE
      HC
    }
    maker
  }

  val data: List[Data] = List(
    new Data("NA12878.WGS.problem",
      new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam"),
      null,
      List("1:142683008-144886616")),
      //List("1:142683008-148886616")),
    new Data("NA12878.WGS",
      new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam"),
      null,
      List("1:10,000,000-11,000,000")),
    new Data("NA12878.WEx",
      new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.bam"),
      new File("/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"),
      List("20")),
    new Data("1KG.1100samples",
      new File("/humgen/1kg/DCC/ftp/20120522.chrom20.bam.list"),
      null,
      List("20:10,000,000-10,010,000")),
    new Data("1KG.300samples",
      new File("20120522.chrom20.bam.300samples.list"),
      null,
      List("20:10,000,000-10,100,000"))
    )

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO"
    this.reference_sequence = makeResource(b37_FILENAME)
    this.memoryLimit = 16
  }

  def script() {
    // iterate over GATK's and data sets
    for ( iteration <- 1 to iterations ) {
      for ( datum <- data ) {
        logger.info("Datum %s %s %s %s".format(datum.name, datum.intervalsFile, datum.intervalsString, datum.bamList))
        for ( (name, maker) <- assessments) {
          val ARTWalker = maker()
          ARTWalker.input_file :+= datum.bamList
          ARTWalker.configureJobReport(Map( "iteration" -> iteration, "analysisName" -> name, "data" -> datum.name))
          datum.addIntervals(ARTWalker)
          add(ARTWalker)
        }
      }
    }
  }
}