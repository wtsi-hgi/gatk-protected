package org.broadinstitute.sting.queue.qscripts.misc.chartl.ResolveBamsFromProject

import org.broadinstitute.sting.commandline.Input
import org.apache.commons.io.FilenameUtils
import org.broadinstitute.sting.pipeline.PicardAggregationUtils
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.queue.QScript
import scala.collection.JavaConversions._

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 2/13/12
 * Time: 2:53 PM
 * To change this template use File | Settings | File Templates.
 */

class ResolveBamPathsFromProject extends QScript {
   @Input(doc="Tab separated squid projects and samples. Name must be <projectName>.tsv", shortName="tsv", required=true)
  var projectSampleTsv: File = _

  var projectName: String = null

  def script = {
      projectName = FilenameUtils.getBaseName(projectSampleTsv)

      val picardSamples = PicardAggregationUtils.parseSamples(projectSampleTsv)

      val writeBamList = new ListWriterFunction
      writeBamList.inputFiles = asScalaIterable(PicardAggregationUtils.getSampleBams(picardSamples)).toSeq
      writeBamList.listFile = projectName +".bam.list"
      writeBamList.jobOutputFile = writeBamList.listFile + ".out"
      add(writeBamList)
  }

}