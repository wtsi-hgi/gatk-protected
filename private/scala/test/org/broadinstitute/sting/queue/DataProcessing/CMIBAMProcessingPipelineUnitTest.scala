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
package org.broadinstitute.sting.queue.DataProcessing

import org.broadinstitute.sting.queue.pipeline.{PipelineTestSpec, PipelineTest}
import org.testng.annotations.Test
import org.broadinstitute.sting.BaseTest


class CMIBAMProcessingPipelineUnitTest {

 //   @Test
    def testSimpleBAM {
      val projectName = "test1"
      val testOut = projectName + ".exampleBAM.bam.clean.dedup.recal.bam"
      val spec = new PipelineTestSpec
      spec.name = "CMIBAMProcessingPipeline"
      spec.args = Array(
        " -S private/scala/qscript/org/broadinstitute/sting/queue/qscripts/DataProcessing/CMIBAMProcessingPipeline.scala",
        " -R " + BaseTest.publicTestDir + "exampleFASTA.fasta",
        " -m testdata/CMITestData/metadata.test.txt",
        " -D " + BaseTest.publicTestDir + "exampleDBSNP.vcf",
        " -L testdata/CMITestData/exampleFASTA.interval_list",
        " -bwa testdata/CMITestData/bwa",
        " -call "
        ).mkString
      spec.fileMD5s += testOut -> "45d97df6d291695b92668e8a55c54cd0"
      PipelineTest.executeTest(spec)
    }




}
