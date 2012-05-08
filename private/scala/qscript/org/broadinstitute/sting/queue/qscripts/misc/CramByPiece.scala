/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.queue.qscripts.misc

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.QFunction
import org.broadinstitute.sting.queue.engine.JobRunInfo
import org.broadinstitute.sting.gatk.report.GATKReport
import java.io.{PrintStream, FileOutputStream}
import org.broadinstitute.sting.commandline.Argument._

class CramByPiece extends QScript {
  qscript =>

  @Argument(shortName = "R", doc = "ref", required = false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "cramtools", doc = "cramtools", required = false)
  val cramtools: File = "/Users/depristo/Dropbox/bin/cramtools.jar"

  @Argument(shortName = "I", doc = "I", required = true)
  val BAM: File = null

  @Argument(shortName = "start", doc = "start", required = false)
  val originalStart: Int = 10000000

  @Argument(shortName = "stop", doc = "start", required = false)
  val fullStop: Int = 11000000

  class Datum(val BAM: File, val CRAM: File, val start: Int, val step: Int) {
    def originalSize:java.lang.Long = BAM.length()
    def cramSize:java.lang.Long = CRAM.length()
    def compressionRatio:java.lang.Double = BAM.length() / CRAM.length()
  }

  var data: List[Datum] = List()

  def script() {
    for ( step <- List(1000000, 100000, 10000)) {
      var start = originalStart
      while ( start < fullStop ) {
        val stop = start + step
        val partialBAM = swapExt(BAM, ".bam", ".start_%d.step_%d.bam".format(start, step))
        val CRAM = swapExt(partialBAM, ".bam", ".cram")
        val myIntervals = "%s:%d-%d".format("20", start, stop)
        add(new PrintRegion(BAM, partialBAM, myIntervals))
        add(new CRAM(partialBAM, CRAM))
        start += step
        data :+= new Datum(partialBAM, CRAM, start, step)
      }
    }
  }

  override def onExecutionDone(jobs: Map[QFunction, JobRunInfo], success: Boolean) {
    val report = GATKReport.newSimpleReport("CompressionSizes", "start", "step", "original.bam.size", "cram.size", "compression.ratio", "original.bam", "cram");
    for ( d: Datum <- data ) {
      System.out.printf("%s %s %d %d%n".format(d.BAM, d.CRAM, d.start, d.step))
      report.addRow(new Integer(d.start), new Integer(d.step), d.originalSize, d.cramSize, d.compressionRatio, d.BAM, d.CRAM)
    }
    report.print(new PrintStream(new FileOutputStream(new File("sizes.gatkreport.txt"))))
  }

  class PrintRegion(@Input fullBAM: File, @Input partialBAM: File, myIntervals: String) extends PrintReads {
    this.reference_sequence = referenceFile
    this.intervalsString +:= myIntervals
    this.input_file +:= fullBAM
    this.out = partialBAM
  }

  class CRAM(@Input var BAM: File, @Output var CRAM: File) extends CommandLineFunction {
    def commandLine = "java -Xmx1g -jar %s cram --input-bam-file %s --reference-fasta-file %s --output-cram-file %s --capture-all-quality-scores --include-unmapped-reads".format(cramtools, BAM, referenceFile, CRAM)
  }
}
