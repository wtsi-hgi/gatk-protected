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

package org.broadinstitute.sting.queue.qscripts.performance

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.function.QFunction

class BCFvsVCFPerformance extends QScript {
  @Argument(shortName = "BUNDLE", doc = "Directory holding all of our data files", required = false)
  val BUNDLE_DIR: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current")

  @Argument(shortName = "dataDir", doc = "Directory to hold temporary VCF/BCF files", required = false)
  val DATA_DIR: File = new File("./results")

  @Argument(shortName = "testVCF", doc = "VCF files for testing", required = false)
  val TEST_VCFs: List[File] = List(new File("/home/unix/depristo/dev/localData/ESP.chr20.vcf"))

  @Argument(shortName = "iterations", doc = "Number of iterations we should execute", required = false)
  val iterations: Int = 1

  @Argument(shortName = "myIntervals", doc = "Intervals", required = false)
  val myIntervals: String = null

  def withBundle(filename: String): File = {
    new File(BUNDLE_DIR.getAbsolutePath + "/" + filename)
  }

  def referenceFile = withBundle("human_g1k_v37.fasta")

  var TEST_RUNNERS: List[(File) => QFunction] = List()
  def addRunner(runner: (File) => QFunction) {
    TEST_RUNNERS +:= runner
  }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO"
    this.reference_sequence = referenceFile
    this.memoryLimit = 2
    if ( myIntervals != null ) {
      this.intervalsString = Seq(myIntervals)
    }
  }

  def script = {
    // performance of converting VCFs and BCFs amoung each other
    for ( rawVCF <- TEST_VCFs ) {
      val testVCF = rawVCFToTestVCF(rawVCF, true)
      val testBCF = vcfToBCF(testVCF)

      for (iteration <- 1 until iterations + 1) {
        for ( inputType <- List("vcf", "bcf")) {
          // test conversion from genotypes files to sites file
          val input = if ( inputType == "vcf" ) testVCF else testBCF
          val sv = new SelectVariants() with UNIVERSAL_GATK_ARGS
          sv.V = input
          sv.analysisName = "genotypesToSites"
          sv.out = swapExt(DATA_DIR, testVCF, ".sites." + inputType, outExt)
          sv.configureJobReport(Map( "iteration" -> iteration, "inputType" -> inputType))
          add(sv)

          // test conversion from XCF -> XCF
          for ( forceGenotypesDecode <- List(true, false)) {
            for ( fastGenotypes <- List(true, false)) {
              for ( outputType <- List("vcf", "bcf")) {
                val sv = new SelectVariants() with UNIVERSAL_GATK_ARGS
                sv.V = input
                sv.analysisName = "fileTypeConversion"
                val outExt = ".it_%d.decode_%b.fastGT_%b.%s".format(iteration, forceGenotypesDecode, fastGenotypes, outputType)
                sv.out = swapExt(DATA_DIR, testVCF, "." + inputType, outExt)
                sv.useFastGenotypes = fastGenotypes
                sv.forceGenotypesDecode = forceGenotypesDecode
                sv.configureJobReport(Map(
                  "iteration" -> iteration,
                  "forceGenotypesDecode" -> forceGenotypesDecode,
                  "fastGenotypes" -> fastGenotypes,
                  "outputType" -> outputType,
                  "inputType" -> inputType))
                add(sv)
              }
            }
          }
        }
      }



      // generic tests of VCF vs. BCF inputs
      for ( includeGenotypes <- List(true, false)) {
        // create the corresponding BCF
        val testVCF = rawVCFToTestVCF(rawVCF, includeGenotypes)
        val testBCF = vcfToBCF(testVCF)

        // now run each test file
        for (iteration <- 1 until iterations + 1) {
          for ( makeRunner <- TEST_RUNNERS ) {
            setupRunner(makeRunner(testVCF), "VCF", iteration, includeGenotypes)
            setupRunner(makeRunner(testBCF), "BCF", iteration, includeGenotypes)
          }
        }
      }
    }
  }

  def setupRunner(job: QFunction, name: String, iteration: Int, includesGenotypes: Boolean) {
    job.configureJobReport(Map(
      "includesGenotypes" -> includesGenotypes,
      "iteration" -> iteration,
      "inputType" -> name))
    add(job)
  }

  def rawVCFToTestVCF(rawVCF: File, includeGenotypes: Boolean): File = {
    val sv = new SelectVariants() with UNIVERSAL_GATK_ARGS
    sv.V = rawVCF
    if ( ! includeGenotypes )
      sv.sites_only = true
    val name = if ( includeGenotypes ) ".test.gt.vcf" else ".test.nogt.vcf"
    sv.out = swapExt(DATA_DIR, rawVCF, ".vcf", name)
    add(sv)
    sv.out
  }

  def vcfToBCF(testVCF: File): File = {
    val sv = new SelectVariants() with UNIVERSAL_GATK_ARGS
    sv.V = testVCF
    sv.out = swapExt(DATA_DIR, testVCF, ".vcf", ".bcf")
    add(sv)
    sv.out
  }

  /**
   * Create a test that just counts records in the input file
   * @param inputFile
   */
  def makeCountVariants(inputFile: File): QFunction = {
    val cr = new CountRODs() with UNIVERSAL_GATK_ARGS
    cr.analysisName = "countRODs"
    cr.rod :+= inputFile
    cr
  }

  /**
   * Create a test that computes per-sample metrics with VariantEval
   * @param inputFile
   */
  def makeVariantSummaryTest(inputFile: File): QFunction = {
    val ev = new VariantEval() with UNIVERSAL_GATK_ARGS
    ev.analysisName = "VariantEval"
    ev.eval :+= inputFile
    ev.noEV = true
    ev.EV = Seq("VariantSummary")
    ev.noST = true
    ev
  }

  addRunner(makeCountVariants)
  addRunner(makeVariantSummaryTest)
}

