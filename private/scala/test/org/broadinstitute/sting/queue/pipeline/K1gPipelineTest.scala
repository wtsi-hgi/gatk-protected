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

package org.broadinstitute.sting.queue.pipeline

import org.broadinstitute.sting.BaseTest
import org.broadinstitute.sting.queue.util.Logging
import java.io.{FileNotFoundException, File}
import org.broadinstitute.sting.pipeline.PicardAggregationUtils

object K1gPipelineTest extends BaseTest with Logging {

  case class K1gBam(project: String, sample: String) {
    override val toString = project + "/" + sample
  }

  /** 1000G BAMs used for validation */
  val k1gBams = List(
    new K1gBam("C474", "NA19651"),
    new K1gBam("C474", "NA19655"),
    new K1gBam("C474", "NA19669"),
    new K1gBam("C454", "NA19834"),
    new K1gBam("C460", "HG01440"),
    new K1gBam("C456", "NA12342"),
    new K1gBam("C456", "NA12748"),
    new K1gBam("C474", "NA19649"),
    new K1gBam("C474", "NA19652"),
    new K1gBam("C474", "NA19654"))

  validateK1gBams()

   /**
   * Throws an exception if any of the 1000G bams do not exist.
   */
  private def validateK1gBams() {
    var missingBams: List[String] = Nil
    for (k1gBam <- k1gBams) {
      try {
        getPicardBam(k1gBam)
      } catch {
        case e: FileNotFoundException =>
          missingBams :+= k1gBam.toString
      }
    }
    if (missingBams.size > 0) {
      val nl = "%n".format()
      throw new FileNotFoundException("The following 1000G bam files are missing.%n%s".format(missingBams.mkString(nl)))
    }
  }

  private def getPicardBam(k1gBam: K1gBam): File =
    new File(PicardAggregationUtils.getSampleBam(k1gBam.project, k1gBam.sample))
}
