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
import org.broadinstitute.sting.queue.util.QJobReport

class RodPerformanceGoals extends QScript {
  @Argument(shortName = "BUNDLE", doc = "Directory holding all of our data files", required=false)
  val BUNDLE_DIR: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current")

  @Argument(shortName = "dataDir", doc = "Directory holding all of our data files", required=false)
  val DATA_DIR: File = new File("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/rodPerformanceGoals");

  @Argument(shortName = "short", doc = "If provided, we will use the shorter tests", required=false)
  val SHORT: Boolean = false;


  def withBundle(filename: String): File = {
    new File(BUNDLE_DIR.getAbsolutePath + "/b37/" + filename)
  }

  def referenceFile = withBundle("human_g1k_v37.fasta")
  def dbsnp = withBundle("dbsnp_132.b37.vcf")
  def bam = withBundle("NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam")

  def OMNI_GENOTYPES = withBundle("1000G_omni2.5.b37.vcf")
  def OMNI_SITES = withBundle("1000G_omni2.5.b37.sites.vcf")

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = referenceFile;
    this.memoryLimit = 2
    if ( SHORT )
      this.intervalsString = List("20:1-10000")
  }

  def script = {
    //countCovariatesTest()
    sitesVsGenotypesTest()
    //bigCombineVariantsTest()
    //lowLevelTribbleVsGATK()
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
  def countCovariatesTest() {
    for ( usedbsnp <- List(true, false))
      for ( nt <- List(1, 8) ) {
        val cc = new CountCovariates() with UNIVERSAL_GATK_ARGS with QJobReport
        cc.setJobLogging("CountCovariates", Map("nt" -> nt, "dbsnp" -> usedbsnp))
        cc.analysisName = "CountCovariatesTest"
        cc.input_file :+= bam
        cc.recal_file = "/dev/null"
        cc.nt = nt
        if ( usedbsnp ) {
          cc.knownSites :+= dbsnp
        } else {
          cc.run_without_dbsnp_potentially_ruining_quality = true;
        }

        add(cc)
      }
  }

  /**
   * What's the relative performance of running CountRods, which doesn't do anything with the
   * tribble feature, when running on the OMNI VCF with sites only vs. with genotypes
   */
  def sitesVsGenotypesTest() {
    for ( vcf <- List(OMNI_SITES, OMNI_GENOTYPES) ) {
      val cr = new CountRODs() with UNIVERSAL_GATK_ARGS with QJobReport
      cr.setJobLogging("SitesVsGenotypes", Map("includesGenotypes" -> (vcf == OMNI_GENOTYPES)))
      cr.analysisName = "SitesVsGenotypesTest"
      cr.scatterCount = 2
      cr.rod :+= vcf
      add(cr)
    }
  }

  /**
   * Compares the performance of GATK CombineVariants to cat-ing VCFs and grep -v #
   */
  def bigCombineVariantsTest() {
    def chunkFile(i: Int): File = new File(DATA_DIR.getAbsolutePath + "/chunks/chunk_" + i + ".vcf")
    val vcfs: List[File] = List.range(1, 50).map(chunkFile(_))
    // cat
    val cgc = new CatGrepCombineVCFs(vcfs) with QJobReport
    cgc.analysisName = "BigCombineVariantsTest"
    cgc.setJobLogging("BigCombine", Map("mode" -> "CatGrep"))
    add( cgc )

    // combine variants
    val cv = new CombineVariants() with UNIVERSAL_GATK_ARGS with QJobReport
    cv.setJobLogging("BigCombine", Map("mode" -> "CombineVariants"))
    cv.analysisName = "BigCombineVariantsTest"
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

  def lowLevelTribbleVsGATK() {
    val justTribble = org.broadinstitute.sting.gatk.walkers.performance.ProfileRodSystem.ProfileType.JUST_TRIBBLE_DECODE
    for ( mode <- List(org.broadinstitute.sting.gatk.walkers.performance.ProfileRodSystem.ProfileType.JUST_TRIBBLE_DECODE,
                       org.broadinstitute.sting.gatk.walkers.performance.ProfileRodSystem.ProfileType.JUST_GATK)) {
      val prs = new ProfileRodSystem() with UNIVERSAL_GATK_ARGS with QJobReport
      prs.intervalsString = null
      prs.analysisName = "lowLevelTribbleVsGATK"
      prs.setJobLogging("LogLevelTribbleVsGATK", Map("mode" -> (if (mode == justTribble) "Tribble" else "GATK")))
      prs.vcf = OMNI_SITES
      prs.mode = mode
      prs.out = "profile.rod." + mode + ".txt"
      add(prs)
    }
  }
}

