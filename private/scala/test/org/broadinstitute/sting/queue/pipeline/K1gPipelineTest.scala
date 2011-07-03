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
import collection.JavaConversions._
import org.broadinstitute.sting.pipeline.{PicardAggregationUtils, Pipeline, PipelineProject, PipelineSample}

object K1gPipelineTest extends BaseTest with Logging {

  case class K1gBam(squidId: String, sampleId: String, version: Int)

  /** 1000G BAMs used for validation */
  val k1gBams = List(
    new K1gBam("C474", "NA19651", 2),
    new K1gBam("C474", "NA19655", 2),
    new K1gBam("C474", "NA19669", 2),
    new K1gBam("C454", "NA19834", 2),
    new K1gBam("C460", "HG01440", 2),
    new K1gBam("C456", "NA12342", 2),
    new K1gBam("C456", "NA12748", 2),
    new K1gBam("C474", "NA19649", 2),
    new K1gBam("C474", "NA19652", 2),
    new K1gBam("C474", "NA19654", 2))

  validateK1gBams()

  /**
   * Creates a new pipeline from a project.
   * @param project Pipeline project info.
   * @param samples List of samples.
   * @return a new pipeline project.
   */
  def createPipeline(project: PipelineProject, samples: List[PipelineSample]) = {
    val pipeline = new Pipeline
    pipeline.setProject(project)
    pipeline.setSamples(samples)
    pipeline
  }

  /**
   * Creates a new pipeline project for hg19 with b37 132 dbsnp for genotyping, and b37 129 dbsnp for eval.
   * @param projectName Name of the project.
   * @param intervals The intervals file to use.
   * @return a new pipeline project.
   */
  def createHg19Project(projectName: String, intervals: String) = {
    val project = new PipelineProject
    project.setName(projectName)
    project.setReferenceFile(new File(BaseTest.hg19Reference))
    project.setGenotypeDbsnp(new File(BaseTest.b37dbSNP132))
    project.setEvalDbsnp(new File(BaseTest.b37dbSNP129))
    project.setRefseqTable(new File(BaseTest.hg19Refseq))
    project.setIntervalList(new File(intervals))
    project
  }

  /**
   * Creates a 1000G pipeline sample from one of the bams.
   * @param idPrefix Text to prepend to the sample name.
   * @param k1gBam bam to create the sample for.
   * @return the created pipeline sample.
   */
  def createK1gSample(idPrefix: String, k1gBam: K1gBam) = {
    val sample = new PipelineSample
    sample.setId(idPrefix + "_" + k1gBam.sampleId)
    sample.setBamFiles(Map("cleaned" -> getPicardBam(k1gBam)))
    sample
  }

   /**
   * Throws an exception if any of the 1000G bams do not exist and warns if they are out of date.
   */
  private def validateK1gBams() {
    var missingBams = List.empty[File]
    for (k1gBam <- k1gBams) {
      val latest = getLatestVersion(k1gBam)
      val bam = getPicardBam(k1gBam)
      if (k1gBam.version != latest)
        logger.warn("1000G bam is not the latest version %d: %s".format(latest, k1gBam))
      if (!bam.exists)
        missingBams :+= bam
    }
    if (missingBams.size > 0) {
      val nl = "%n".format()
      throw new FileNotFoundException("The following 1000G bam files are missing.%n%s".format(missingBams.mkString(nl)))
    }
  }

  private def getPicardBam(k1gBam: K1gBam): File =
    new File(PicardAggregationUtils.getSampleBam(k1gBam.squidId, k1gBam.sampleId, k1gBam.version))

  private def getLatestVersion(k1gBam: K1gBam): Int =
    PicardAggregationUtils.getLatestVersion(k1gBam.squidId, k1gBam.sampleId, k1gBam.version)
}
