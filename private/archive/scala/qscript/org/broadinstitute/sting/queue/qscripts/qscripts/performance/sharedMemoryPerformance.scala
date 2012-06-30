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
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.walkers.recalibration.CountCovariatesWalker

class SharedMemoryPerformance extends QScript {
  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "bam", doc = "foo", required=true)
  val BAM: File = null;

  @Argument(shortName = "dbsnp", doc = "foo", required=true)
  val dbsnp: File = null;

  @Argument(shortName = "intervals", doc="intervals", required=false)
  val TARGET_INTERVAL: String = null;

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    //this.jarFile = gatkJarFile;
    if ( TARGET_INTERVAL != null )
      this.intervalsString = List(new File(TARGET_INTERVAL));
    this.reference_sequence = referenceFile;
    this.mem oryLimit = 4
  }

  def script = {
    for ( usedbsnp <- List(true, false))
    for (nt <- List(1, 2, 3, 4, 6, 8, 10, 12, 14, 16)) {
      val cc = new CountCovariates() with UNIVERSAL_GATK_ARGS
      cc.input_file :+= BAM
      cc.o = "nt" + nt + "_dbsnp" + usedbsnp + ".csv"
      cc.nt = nt
      if ( usedbsnp ) {
        cc.knownSites :+= dbsnp
      } else {
        cc.run_without_dbsnp_potentially_ruining_quality = true;
      }
      add(cc)
    }
  }
}

