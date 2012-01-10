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

import org.testng.annotations.{DataProvider, Test}
import org.broadinstitute.sting.BaseTest
import org.apache.commons.io.FileUtils
import org.broadinstitute.sting.pipeline.PicardAggregationUtils
import collection.JavaConversions._

class HybridSelectionPipelineTest {
  def datasets = List(k1gChr20Dataset)

  val k1gChr20Dataset = {
    val dataset = newK1gDataset("Barcoded_1000G_WEx_chr20", BaseTest.hg19Chr20Intervals)

    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.all.all.all", "nCalledLoci", 1464)
    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.all.known.all", "nCalledLoci", 1196)
    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.all.novel.all", "nCalledLoci", 268)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.all.all.all", "tiTvRatio", 3.56)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.all.known.all", "tiTvRatio", 3.80)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.all.novel.all", "tiTvRatio", 2.72)

    dataset
  }

  def newK1gDataset(projectName: String, intervals: String) = {
    val bams = K1gPipelineTest.k1gBams.map(ps => PicardAggregationUtils.getSampleBam(ps.project, ps.sample))
    new PipelineDataset(projectName, intervals, bams)
  }

  @DataProvider(name="datasets")//, parallel=true)
  final def convertDatasets: Array[Array[AnyRef]] =
    datasets.map(dataset => Array(dataset.asInstanceOf[AnyRef])).toArray

  @Test(dataProvider="datasets")
  def testHybridSelectionPipeline(dataset: PipelineDataset) {
    val testName = "HybridSelectionPipeline-" + dataset.projectName
    val bamList = writeBamList(dataset.projectName + ".bam.list", dataset.bams)

    // Run the pipeline with the expected inputs.
    val pipelineCommand =
      ("-retry 1" +
        " -S private/scala/qscript/org/broadinstitute/sting/queue/qscripts/pipeline/HybridSelectionPipeline.scala" +
        " -I %s" +
        " -L %s" +
        " -varFilter HARD")
        .format(bamList, dataset.intervals)

    val pipelineSpec = new PipelineTestSpec
    pipelineSpec.name = testName
    pipelineSpec.args = pipelineCommand
    pipelineSpec.jobQueue = dataset.jobQueue

    pipelineSpec.evalSpec = new PipelineTestEvalSpec
    pipelineSpec.evalSpec.evalReport = dataset.projectName + ".by_sample.eval"
    pipelineSpec.evalSpec.validations = dataset.validations

    PipelineTest.executeTest(pipelineSpec)
  }

  private def writeBamList(fileName: String, bams: List[String]) = {
    val bamList = BaseTest.createNetworkTempFile(fileName)
    FileUtils.writeLines(bamList, bams)
    bamList
  }

  class PipelineDataset(var projectName: String,
                        var intervals: String,
                        var bams: List[String],
                        var validations: List[PipelineValidation[_]] = Nil,
                        var jobQueue: String = null) {
    override def toString = projectName
  }
}
