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

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import java.lang.Math
import org.broadinstitute.sting.utils.PathUtils
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel
import org.broadinstitute.sting.utils.baq.BAQ
import org.broadinstitute.sting.commandline.ClassType

class GATKPerformanceOverTime extends QScript {
  @Argument(shortName = "results", doc="results", required=false)
  val resultsDir: File = new File("runResults")

  @Argument(shortName = "resources", doc="resources", required=true)
  val resourcesDir: String = ""

  @Argument(shortName = "myJarFile", doc="Path to the current GATK jar file", required=true)
  val myJarFile: File = null

  @Argument(shortName = "iterations", doc="it", required=false)
  val iterations: Int = 2

  @Argument(shortName = "smallData", doc="it", required=false)
  val smallData: Boolean = false

  @Argument(shortName = "justDeepWGS", doc="it", required=false)
  val justDeepWGS: Boolean = false

  @Argument(shortName = "skipBAQ", doc="it", required=false)
  val skipBAQ: Boolean = false

  @Argument(shortName = "assessment", doc="Which assessments should we run?", required=false)
  val assessmentsArg: Set[String] = Assessment.values map(_.toString)
  var assessments: Set[Assessment.Assessment] = _

  @Argument(shortName = "ntTest", doc="For each value provided we will use -nt VALUE in the multi-threaded tests", required=false)
  @ClassType(classOf[Int])
  val ntTests: List[Int] = List(1, 2, 3, 4, 8, 12, 16)

  @Argument(shortName = "steps", doc="steps", required=false)
  val steps: Int = 10

  @Argument(shortName = "maxNSamples", doc="maxNSamples", required=false)
  val maxNSamples: Int = 1000000

  @Argument(shortName = "manyVersions", doc="manyVersions", required=false)
  val manyVersions: Boolean = false

  val singleTestsPerIteration = 3

  val RECAL_BAM_FILENAME = "CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20GAV.8.bam"
  val RECAL_GATKREPORT_FILENAME = "CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20GAV.8.grp"
  val dbSNP_FILENAME = "dbsnp_132.b37.vcf"
  val BIG_VCF_WITH_GENOTYPES = "ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
  val BIG_VCF_WITH_GENOTYPES_16_COMPATIBLE = new File("/humgen/gsa-hpprojects/GATK/bundle/1.5/b37/1000G_omni2.5.b37.vcf")
  val b37_FILENAME = "human_g1k_v37.fasta"

  val manyDeepExomeSamples = "combined_5022.bam"
  val manyDeepExomeIntervals = "1:1104385-1684501"

  def makeResource(x: String): File = new File("%s/%s".format(resourcesDir, x))
  def makeChunk(x: Int): File = makeResource("chunk_%d.vcf".format(x))
  def COMBINE_FILES: List[File] = Range(1,10).map(makeChunk).toList

  class AssessmentParameters(val name: String,
                             val bamList: File,
                             val fullIntervals: String,
                             val shortIntervals: String,
                             val nSamples: Int,
                             val dcov: Int,
                             val baq: Boolean) {
    def addIntervals(gatkCmd : CommandLineGATK, useFull: Boolean): CommandLineGATK = {
      val intervals = if (useFull && ! smallData) fullIntervals else shortIntervals
      val maybeFile = makeResource(intervals)
      if ( maybeFile.exists() )
        gatkCmd.intervals :+= maybeFile
      else
        gatkCmd.intervalsString :+= intervals
      gatkCmd
    }
  }

  // TODO -- count the number of lines in the bam.list file
  val WGSAssessment = new AssessmentParameters("WGS.multiSample.4x", "wgs.bam.list.local.list", "wgs.bam.list.select.intervals", "20:10000000-11000000", 1103, 50, true)
  val WGSDeepAssessment = new AssessmentParameters("WGS.singleSample.60x", "wgs.deep.bam.list.local.list", "wgs.deep.bam.list.select.intervals", "1", 1, 250, true)
  val WGSDeepAssessmentReduced = new AssessmentParameters("WGS.singleSample.60x.reduced", "CEUTrio.HiSeq.WGS.b37.NA12878.clean.dedup.recal.reduced.bam", "20:10,000,000-10,500,000", "20:10,000,000-10,500,000", 1, 250, true)
  val WExAssessment = new AssessmentParameters("WEx.multiSample.150x", "wex.bam.list.local.list", "wex.bam.list.select.intervals", "wex.bam.list.small.intervals", 140, 500, true)

  val GATK_RELEASE_DIR = new File("/humgen/gsa-hpprojects/GATK/bin/")
  var GATKs: Map[String, File] =
    Map(
      "v2.cur" -> myJarFile, // TODO -- how do I get this value?
      "v2.5" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "2.5"),
      "v2.4" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "2.4"),
      "v2.3" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "2.3"),
      "v2.0" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "2.0"),
      "v1.6" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "1.6")
    )

  object Assessment extends Enumeration {
    type Assessment = Value
    val UG, UG_NT, CL, CL_NT, CV, CV_NT, VE, VE_NT, SV, BQSR_NT, PRINT_READS_NT, MANY_SAMPLES_NT, HC_VS_UG = Value
  }

  val NCT_ASSESSMENTS = List(Assessment.UG_NT, Assessment.CL_NT, Assessment.BQSR_NT, Assessment.PRINT_READS_NT, Assessment.MANY_SAMPLES_NT)
  def supportsNCT(assessment: Assessment.Assessment) = NCT_ASSESSMENTS contains assessment

  val NO_NT_ASSESSMENTS = List(Assessment.PRINT_READS_NT, Assessment.MANY_SAMPLES_NT)
  def supportsNT(assessment: Assessment.Assessment) = ! (NO_NT_ASSESSMENTS contains assessment)

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO"
    this.reference_sequence = makeResource(b37_FILENAME)
    this.memoryLimit = 4
  }

  def script() {
    assessments = assessmentsArg.map(Assessment.withName(_))

    if ( manyVersions ) {
      GATKs = GATKs ++ Map(
        "v2.3" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "2.3"),
        "v2.2" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "2.2"),
        "v2.1" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "2.1"),
        "v1.5" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "1.5"),
        "v1.4" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "1.4"),
        "v1.3" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "1.3"),
        "v1.2" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "1.2"),
        "v1.1" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "1.1"),
        "v1.0" -> PathUtils.findMostRecentGATKVersion(GATK_RELEASE_DIR, "1.0")
      )
    }

    if ( ! resultsDir.exists ) resultsDir.mkdirs()

    // iterate over GATK's and data sets
    for ( iteration <- 0 until iterations ) {
      for ( (gatkName, gatkJar) <- GATKs ) {

        enqueueCommandsForEachAssessment(iteration, gatkName, gatkJar)
        enqueueCommandsForHCvsUG(iteration, gatkName, gatkJar)
        enqueueMultiThreadedCommands(iteration, gatkName, gatkJar)
        enqueueSingleTestCommands(iteration, gatkName, gatkJar)
      }
    }
  }

  def enqueueMultiThreadedCommands(iteration: Int, gatkName: String, gatkJar: File) {
    // GATK v2 specific tests
    if ( assessments.contains(Assessment.CV_NT) ) {
      if ( gatkName.contains("v2") ) {
        for ( outputBCF <- List(true, false) ) {
          def makeCV(): CommandLineGATK = {
            val outputName = if ( outputBCF ) "bcf" else "vcf"
            val CV = new CombineVariants with UNIVERSAL_GATK_ARGS
            CV.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName, "assessment" -> outputName))
            CV.jarFile = gatkJar
            CV.variant = List(makeResource(BIG_VCF_WITH_GENOTYPES))
            CV.out = new File("/dev/null")
            CV.bcf = outputBCF
            CV
          }
          addMultiThreadedTest(gatkName, Assessment.CV_NT, makeCV)
        }
      }
    }

    if ( assessments.contains(Assessment.BQSR_NT) ) {
      def makeBQSR(): BaseRecalibrator = {
        val BQSR = new MyBaseRecalibrator(gatkName.contains("v2"))
        BQSR.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName, "assessment" -> "20GAV.8.bam"))
        BQSR.jarFile = gatkJar
        BQSR
      }
      addMultiThreadedTest(gatkName, Assessment.BQSR_NT, makeBQSR, 8) // max nt until BQSR is performant
    }

    if ( assessments.contains(Assessment.MANY_SAMPLES_NT) ) {
      def makeUG(): UnifiedGenotyper = {
        val ug = new Call(makeResource(manyDeepExomeSamples), 1, "manyDeepExomes", false) with UNIVERSAL_GATK_ARGS
        ug.intervalsString :+= manyDeepExomeIntervals
        ug.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName, "assessment" -> "manyDeepExomes"))
        ug.jarFile = gatkJar
        ug.glm = GenotypeLikelihoodsCalculationModel.Model.SNP
        ug.memoryLimit = 16
        ug
      }
      addMultiThreadedTest(gatkName, Assessment.MANY_SAMPLES_NT, makeUG, scaleMem = false)
    }

    if ( assessments.contains(Assessment.PRINT_READS_NT) ) {
      for ( assessment <- List("BQSR", "BAQ") ) {
        def makePrintReads(): PrintReads = {
          val PR = new PrintReads with UNIVERSAL_GATK_ARGS
          PR.intervalsString = List("1", "2", "3", "4", "5")
          PR.input_file :+= makeResource(RECAL_BAM_FILENAME)
          PR.out = new File("/dev/null")
          PR.memoryLimit = 4
          if ( assessment.equals("BQSR") )
            PR.baq = BAQ.CalculationMode.RECALCULATE
          else
            PR.BQSR = makeResource(RECAL_GATKREPORT_FILENAME)
          PR.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName, "assessment" -> assessment))
          PR.jarFile = gatkJar
          PR
        }
        addMultiThreadedTest(gatkName, Assessment.PRINT_READS_NT, makePrintReads, scaleMem = false)
      }
    }
  }

  def enqueueSingleTestCommands(iteration: Int, gatkName: String, gatkJar: File) {
    for ( subiteration <- 0 until singleTestsPerIteration ) {
      trait VersionOverrides extends CommandLineGATK {
        this.jarFile = gatkJar
        this.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName))
      }

      val CV = new CombineVariants with UNIVERSAL_GATK_ARGS with VersionOverrides
      CV.variant = COMBINE_FILES
      CV.intervalsString = List("1", "2", "3", "4", "5")
      CV.out = new File("/dev/null")
      if ( assessments.contains(Assessment.CV) )
        addGATKCommand(CV)

      val SV = new SelectVariants with UNIVERSAL_GATK_ARGS with VersionOverrides
      SV.variant = BIG_VCF_WITH_GENOTYPES_16_COMPATIBLE
      SV.sample_name = List("HG00096") // IMPORTANT THAT THIS SAMPLE BE IN CHUNK ONE
      SV.out = new File("/dev/null")
      if ( assessments.contains(Assessment.SV) )
        addGATKCommand(SV)

      def makeVE(): CommandLineGATK = {
        val VE = new VariantEval with UNIVERSAL_GATK_ARGS with VersionOverrides
        VE.eval :+= BIG_VCF_WITH_GENOTYPES_16_COMPATIBLE
        VE.out = new File("/dev/null")
        VE.comp :+= new TaggedFile(makeResource(dbSNP_FILENAME), "dbSNP")
        VE.addJobReportBinding("assessment", BIG_VCF_WITH_GENOTYPES_16_COMPATIBLE.getName)
        VE
      }

      if ( assessments.contains(Assessment.VE) ) {
        addGATKCommand(makeVE())
      }

      if ( assessments.contains(Assessment.VE_NT) )
        addMultiThreadedTest(gatkName, Assessment.VE_NT, makeVE)
    }
  }

  def enqueueCommandsForEachAssessment(iteration: Int, gatkName: String, gatkJar: File) {
    val dataSets = if ( justDeepWGS ) List(WGSDeepAssessment) else List(WGSAssessment, WGSDeepAssessment, WExAssessment)

    for ( assess <- dataSets ) {
      for (nSamples <- divideSamples(assess.nSamples) ) {
        val sublist = new SliceList(assess.name, nSamples, makeResource(assess.bamList))
        if ( iteration == 0 ) add(sublist) // todo - remove condition when Queue bug is fixed
        val name: String = "%s/assess.%s_gatk.%s_iter.%d".format(resultsDir, assess.name, gatkName, iteration)

        trait VersionOverrides extends CommandLineGATK {
          this.jarFile = gatkJar
          this.dcov = assess.dcov

          this.configureJobReport(Map(
            "iteration" -> iteration,
            "gatk" -> gatkName,
            "nSamples" -> nSamples,
            "assessment" -> assess.name))
        }

        // SNP calling
        if ( assessments.contains(Assessment.UG) )
          addGATKCommand(assess.addIntervals(new Call(sublist.list, nSamples, name, assess.baq) with VersionOverrides, false))
        if ( assessments.contains(Assessment.UG_NT) && nSamples == assess.nSamples )
          addMultiThreadedTest(gatkName, Assessment.UG_NT, () => assess.addIntervals(new Call(sublist.list, nSamples, name, assess.baq) with VersionOverrides, false))

        // CountLoci
        if ( assessments.contains(Assessment.CL) )
          addGATKCommand(assess.addIntervals(new MyCountLoci(sublist.list, nSamples, name) with VersionOverrides, true))
        if ( assessments.contains(Assessment.CL_NT) && nSamples == assess.nSamples )
          addMultiThreadedTest(gatkName, Assessment.CL_NT, () => assess.addIntervals(new MyCountLoci(sublist.list, nSamples, name) with VersionOverrides, true))
      }
    }
  }

  def enqueueCommandsForHCvsUG(iteration: Int, gatkName: String, gatkJar: File) {
    val G1K4x = new AssessmentParameters("WGS.G1K4x.1100samples", "20120522.chrom20.20_10000000_10100000.bam", "20:10,000,000-10,100,000", "20:10,000,000-10,010,000", 1, 50, true)
    val dataSets = if ( justDeepWGS ) List(WGSDeepAssessment) else List(G1K4x, WGSDeepAssessment, WGSDeepAssessmentReduced, WExAssessment)

    val intervalsForAssessments = Map(
      G1K4x -> "20:10,000,000-10,010,000",
      WGSDeepAssessment -> "1:10,000,000-11,000,000",
      WGSDeepAssessmentReduced -> WGSDeepAssessmentReduced.fullIntervals,
      WExAssessment -> makeResource("wex.bam.list.tiny.intervals")
    )

    if ( assessments.contains(Assessment.HC_VS_UG) ) {
      for ( assess <- dataSets ) {
        for ( tool <- List("HaplotypeCaller", "UnifiedGenotyper") ) {
          val bams = {
            if ( assess.bamList.getAbsolutePath.endsWith(".list") ) {
              val sublist = new SliceList(assess.name, assess.nSamples, makeResource(assess.bamList))
              if ( iteration == 0 ) add(sublist)
              sublist.list
            } else {
              makeResource(assess.bamList)
            }
          }
          val name: String = "%s/assess.%s_gatk.%s_iter.hc_vs_ug.%d".format(resultsDir, assess.name, gatkName, iteration)

          trait VersionOverrides extends CommandLineGATK {
            this.jarFile = gatkJar
            this.intervalsString :+= intervalsForAssessments.get(assess).get
            this.memoryLimit = 8

            this.configureJobReport(Map(
              "iteration" -> iteration,
              "gatk" -> gatkName,
              "tool" -> tool,
              "assessment" -> assess.name))
            this.analysisName = "HCvsUG"
          }

          if ( tool == "UnifiedGenotyper" )
            addGATKCommand(new Call(bams, assess.nSamples, name, assess.baq) with VersionOverrides)
          else if ( gatkName.contains("v2") ) // must be HC, only do it if version is 2+
            addGATKCommand(new HCCall(bams, assess.nSamples, name) with VersionOverrides)
        }
      }
    }
  }

  def addMultiThreadedTest(gatkName: String,
                           assessment: Assessment.Assessment,
                           makeCommand: () => CommandLineGATK,
                           maxNT : Int = 1000,
                           scaleMem : Boolean = true) {
    if ( ntTests.size >= 1 ) {
      for ( nt <- ntTests ) {
        if ( nt <= maxNT ) {
          for ( useNT <- List(true, false) ) {
            if ( ( useNT && supportsNT(assessment)) ||
              (! useNT && supportsNCT(assessment) && gatkName.contains("v2.cur") )) {
              // TODO -- fix v2.cur testing
              val cmd = makeCommand()
              if ( useNT )
                cmd.nt = nt
              else
                cmd.nct = nt
              if ( scaleMem )
                cmd.memoryLimit = cmd.memoryLimit * (if ( nt >= 8 ) (if (nt>=16) 4 else 2) else 1)
              cmd.addJobReportBinding("nt", nt)
              cmd.addJobReportBinding("ntType", if ( useNT ) "nt" else "nct")
              cmd.analysisName = cmd.analysisName + ".nt"
              addGATKCommand(cmd)
            }
          }
        }
      }
    }
  }

  def divideSamples(nTotalSamples: Int): List[Int] = {
    val maxLog10: Double = Math.log10(Math.min(maxNSamples, nTotalSamples))
    val stepSize: Double = maxLog10 / steps
    val ten: Double = 10.0
    def deLog(x: Int): Int = Math.round(Math.pow(ten, stepSize * x)).toInt
    dedupe(Range(0, steps+1).map(deLog).toList)
  }

  class Call(@Input(doc="foo") bamList: File, n: Int, name: String, useBaq: Boolean) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.input_file :+= bamList
    this.stand_call_conf = 10.0
    this.glm = GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.baq = if ( ! skipBAQ && useBaq ) org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE else org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
    @Output(doc="foo") var outVCF: File = new File("/dev/null")
    this.o = outVCF
  }

  class HCCall(@Input(doc="foo") bamList: File, n: Int, name: String) extends HaplotypeCaller with UNIVERSAL_GATK_ARGS {
    this.input_file :+= bamList
    @Output(doc="foo") var outVCF: File = new File("/dev/null")
    this.o = outVCF
  }

  class MyCountLoci(@Input(doc="foo") bamList: File, n: Int, name: String) extends CountLoci with UNIVERSAL_GATK_ARGS {
    this.input_file :+= bamList
    @Output(doc="foo") var outFile: File = new File("/dev/null")
    this.o = outFile
  }

  class SliceList(prefix: String, n: Int, @Input bamList: File) extends CommandLineFunction {
    this.analysisName = "SliceList"
    @Output(doc="foo") var list: File = new File("%s/%s.bams.%d.list".format(resultsDir.getPath, prefix, n))
    def commandLine = "head -n %d %s | awk '{print \"%s/\" $1}' > %s".format(n, bamList, resourcesDir, list)
  }

  def addGATKCommand(gatkCmd: CommandLineGATK) {
    add(gatkCmd)
  }

  def dedupe(elements:List[Int]):List[Int] = {
    if (elements.isEmpty)
      elements
    else
      elements.head :: dedupe(for (x <- elements.tail if x != elements.head) yield x)
  }

  /**
   * This is a total abuse of the Queue system.  Override the command line function and remove the
   * -o /dev/null argument which isn't present in BQSR v1.  Otherwise this system magically
   * @param v2
   */
  class MyBaseRecalibrator(val v2: Boolean) extends BaseRecalibrator with UNIVERSAL_GATK_ARGS {
    this.intervalsString = List("1", "2", "3", "4", "5")
    this.knownSites :+= makeResource(dbSNP_FILENAME)
    // must explicitly list the covariates so that BQSR v1 works
    this.input_file :+= makeResource(RECAL_BAM_FILENAME)
    this.out = new File("/dev/null")
    this.memoryLimit = 12

    if ( ! v2 ) {
      this.analysis_type = "CountCovariates" // BQSR v1 name is CountCovariates
    }

    // terrible terrible hack.  Explicitly remove the -o output which isn't present in v1
    override def commandLine(): String = {
      val covariates = "-cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate"
      if ( ! v2 )
        super.commandLine.replace("'-o' '/dev/null'", covariates)
      else
        super.commandLine
    }
  }
}
