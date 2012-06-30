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

package org.broadinstitute.sting.queue.qscripts.performance

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * A simple script that runs a basic GATK tool (CountRODs) to
 * count records in dbSNP that can be used to contrast performance
 * of farm nodes using the QueueJobReport
 */
class AssessFarmPerformance extends QScript {
  @Argument(shortName = "BUNDLE", doc = "Directory holding all of our data files", required=false)
  val BUNDLE_DIR: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current")

  @Argument(shortName = "iterations", doc = "Number of iterations we should execute", required=false)
  val iterations: Int = 1000;

  def withBundle(filename: String): File = {
    new File(BUNDLE_DIR.getAbsolutePath + "/b37/" + filename)
  }

  def referenceFile = withBundle("human_g1k_v37.fasta")
  def dbsnp = withBundle("dbsnp_132.b37.vcf")
  def bam = withBundle("NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam")

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = referenceFile;
    this.memoryLimit = 2
  }

  def script = {
    for ( iteration <- 1 until iterations + 1 ) {
      val cr = new CountRODs() with UNIVERSAL_GATK_ARGS
      cr.configureJobReport(Map("iteration" -> iteration))
      cr.analysisName = "CountingDBSNPRecords"
      cr.rod :+= dbsnp
      //cr.intervalsString = List("1")
      add(cr)
    }
  }
}

