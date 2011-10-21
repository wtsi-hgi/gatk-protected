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

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

class RodPerformanceGoals extends QScript {
  @Argument(shortName = "BUNDLE", doc = "Directory holding all of our data files", required=false)
  val BUNDLE_DIR: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current")

  @Argument(shortName = "dataDir", doc = "Directory holding all of our data files", required=false)
  val DATA_DIR: File = new File("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/rodPerformanceGoals");

  @Argument(shortName = "short", doc = "If provided, we will use the shorter version of the tests", required=false)
  val SHORT: Boolean = false;

  @Argument(shortName = "largeFile", doc = "If provided, we will use the large VCF file to test the GATK vs. Tribble", required=false)
  val LARGE_FILE: Boolean = false;

  @Argument(shortName = "iterations", doc = "Number of iterations we should execute", required=false)
  val iterations: Int = 3;

  @Argument(shortName = "test", doc = "If provided, we will use only the basic test", required=false)
  val TEST: Boolean = false;

  @Argument(shortName = "multithreaded", doc = "If provided, we will include multi-threaded tests but this should only run on hosts with > 8 cores", required=false)
  val multithreaded: Boolean = false;

  @Argument(shortName = "sc", doc = "X", required=false)
  val SC: Int = 1;

  def withBundle(filename: String): File = {
    new File(BUNDLE_DIR.getAbsolutePath + "/b37/" + filename)
  }

  def referenceFile = withBundle("human_g1k_v37.fasta")
  def dbsnp = withBundle("dbsnp_132.b37.vcf")
  def bam = withBundle("NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam")

  def OMNI_GENOTYPES = withBundle("1000G_omni2.5.b37.vcf")
  def OMNI_SITES = withBundle("1000G_omni2.5.b37.sites.vcf")
  def LARGE_1000G_FILE = new File("/humgen/1kg/processing/official_release/phase1/projectConsensus/ALL.wgs.projectConsensus_v2b.20101123.snps.sites.vcf");

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = referenceFile;
    this.memoryLimit = 4
    if ( SHORT )
      this.intervalsString = List("20:1-10000")
  }

  def script = {
    for ( iteration <- 1 until iterations + 1 ) {
      sitesVsGenotypesTest(iteration)
      if ( ! TEST ) {
        countCovariatesTest(iteration)
        bigCombineVariantsTest(iteration)
	if ( LARGE_FILE ) {
	   lowLevelTribbleVsGATK(iteration, LARGE_1000G_FILE)
	} else {
           lowLevelTribbleVsGATK(iteration, OMNI_SITES)
	}
      }
    }
  }

  /**
   * In effect two tests.
   *
   * What's the performance cost of running CountCovariates with dbsnp VCF,
   * compared to running without dbsnp at all?
   *
   * Second, in both cases, what's the performance with -nt 1 vs. -nt 8?  We can
   * then compare scaling w/ w/o dbsnp and w/ and w/o parallelism
   */
  def countCovariatesTest(iteration:Int) {
    for ( usedbsnp <- List(true, false))
      for ( nt <- (if (multithreaded) List(1, 8) else List(1)) ) {
        val cc = new CountCovariates() with UNIVERSAL_GATK_ARGS
        cc.configureJobReport(Map("nt" -> nt, "dbsnp" -> usedbsnp, "iteration" -> iteration))
        cc.analysisName = "CountCovariates"
        cc.input_file :+= bam
        cc.recal_file = "/dev/null"
        cc.nt = nt
        if ( usedbsnp ) {
          cc.knownSites :+= dbsnp
        } else {
          cc.run_without_dbsnp_potentially_ruining_quality = true;
        }

        if ( cc.intervalsString == Nil )
          cc.intervalsString = List("20:1-10,000,000")

        add(cc)
      }
  }

  /**
   * What's the relative performance of running CountRods, which doesn't do anything with the
   * tribble feature, when running on the OMNI VCF with sites only vs. with genotypes
   */
  def sitesVsGenotypesTest(iteration:Int) {
    for ( vcf <- List(OMNI_SITES, OMNI_GENOTYPES) ) {
      val cr = new CountRODs() with UNIVERSAL_GATK_ARGS
      cr.configureJobReport(Map("includesGenotypes" -> (vcf == OMNI_GENOTYPES), "iteration" -> iteration))
      cr.analysisName = "SitesVsGenotypes"
      cr.scatterCount = SC
      cr.rod :+= vcf
      if ( cr.intervalsString == Nil )
        cr.intervalsString = List("1")

      add(cr)
    }
  }

  /**
   * Compares the performance of GATK CombineVariants to cat-ing VCFs and grep -v #
   */
  def bigCombineVariantsTest(iteration:Int) {
    def chunkFile(i: Int): File = new File(DATA_DIR.getAbsolutePath + "/chunks/chunk_" + i + ".vcf")
    val vcfs: List[File] = List.range(1, 50).map(chunkFile(_))
    // cat
    val cgc = new CatGrepCombineVCFs(vcfs)
    cgc.analysisName = "BigCombine"
    cgc.configureJobReport(Map("mode" -> "CatGrep", "iteration" -> iteration))
    add( cgc )

    // combine variants
    val cv = new CombineVariants() with UNIVERSAL_GATK_ARGS
    cv.configureJobReport(Map("mode" -> "CombineVariants", "iteration" -> iteration))
    cv.analysisName = "BigCombine"
    cv.variant = vcfs
    cv.assumeIdenticalSamples = true
    cv.out = "/dev/null"
    add(cv)
  }

  class CatGrepCombineVCFs(vcfs: List[File]) extends CommandLineFunction {
    @Output(doc="foo") var out: File = new File("/dev/null")
    val files = vcfs.map(_.getAbsolutePath).reduceLeft(_ + " " + _)
    def commandLine = "cat %s | grep -v \"#\" > %s".format(files, out)
  }

  def lowLevelTribbleVsGATK(iteration:Int, vcf:File) {
    val justTribbleMode = org.broadinstitute.sting.gatk.walkers.performance.ProfileRodSystem.ProfileType.JUST_TRIBBLE_DECODE
    val justGATKMode = org.broadinstitute.sting.gatk.walkers.performance.ProfileRodSystem.ProfileType.JUST_GATK

    def makeTest(name: String) = {
      val prs = new ProfileRodSystem() with UNIVERSAL_GATK_ARGS
      prs.intervalsString = null
      prs.analysisName = "TribbleVsGATK"
      prs.configureJobReport(Map("mode" -> name, "iteration" -> iteration))
      prs.out = "profile.rod." + name + ".txt"
      prs
    }

    val justGATK = makeTest("GATK")
    justGATK.vcf = vcf
    justGATK.mode = justGATKMode

    val justTribble = makeTest("Tribble")
    justTribble.vcf = vcf
    justTribble.mode = justTribbleMode

    val justGATKStream = makeTest("GATK-STREAM")
    justGATKStream.vcf = TaggedFile(vcf,"storage=STREAM")
    justGATKStream.mode = justGATKMode

    add(justGATK)
    add(justTribble)
    add(justGATKStream)
  }
}

