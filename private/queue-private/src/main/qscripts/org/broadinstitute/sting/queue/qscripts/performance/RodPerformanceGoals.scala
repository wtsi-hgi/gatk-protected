/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
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

