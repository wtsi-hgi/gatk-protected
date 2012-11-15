/**
 * User: dvoet
 */

package org.broadinstitute.sting.queue.qscripts

import us.countmein.queueext._
import org.broadinstitute.sting.queue.extensions.gatk._

import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction
import java.io.{PrintWriter, File}

class TrivialTask extends CmiScript {
  qscript =>


  /****************************************************************************
    * Output Parameters Parameters
    ****************************************************************************/
  @Hidden @Output(fullName="test_output", shortName = "to", doc="test output file", required = false)
  var testOutput: File = _

  def script() {
    testOutput = new File("hello.txt")
    val writer = new PrintWriter(testOutput)
    writer.println("hello world!")
    writer.close()
  }

}
