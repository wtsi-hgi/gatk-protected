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

package org.broadinstitute.sting.queue.qscripts.calling

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import scala.io.Source._

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 3/17/11
 * Time: 11:29 AM
 * To change this template use File | Settings | File Templates.
 */


class Downsampling extends QScript {

  @Input(doc="path to GenomeAnalysisTK.jar", shortName="gatk", required=true)
  var GATKjar: File = _

  @Input(doc="input BAM file - or list of BAM files", shortName="i", required=true)
  var input: File = _

  @Input(doc="target intervals", shortName="t", required=true)
  var targetIntervals: File = _

  @Input(doc="bootstrap number", shortName="b", required=false)
  var bootstrap: Int = 1

  @Input(doc="downsampling step", shortName="ds", required=true)
  var downsamplingStep: Double = _

  @Input(doc="downsampling floor", shortName="df", required=false)
  var downsamplingFloor: Double = 0.0

  @Input(doc="downsampling ceiling", shortName="dc", required=false)
  var downsamplingCeiling: Double = 1.0

  @Input(doc="Reference fasta file", shortName="R", required=false)
  var reference: File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")

  @Input(doc="HapMap file", shortName="H", required=false)
  var hapmap: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf")

  @Input(doc="Omni file", shortName="O", required=false)
  var omni: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf")

  @Input(doc="dbSNP file", shortName="D", required=false)
  var dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_132_b37.leftAligned.vcf")

  @Input(doc="project name", shortName="p", required=false)
  var base: String = "prj"

  def countLines(file: File):Int = {
    var count: Int = 0
    for (l <- fromFile(file).getLines) {
      count = count + 1
    }
    return count
  }

  val queueLogDir: String = ".qlog/"
  val outFile: String = "cov.out"
  val fullCoverageVCF = new File("/humgen/gsa-hpprojects/dev/carneiro/downsampling/analysis/fullcov/fullcov.F1.filtered.vcf")
  val trancheTarget = 99.0

  def script = {
    val nIntervals = math.min(200, countLines(targetIntervals))

    var f: Double = downsamplingCeiling
    var i: Int = 1
    while (f>=downsamplingFloor) {
      var b: Int = bootstrap
      while(b > 0) {
        val file = swapExt(outFile, ".out", ".F" + i + "." + b + ".out")
        add(cov(f, file))
        b = b - 1
      }
      val snp_out = new File(base + ".F" + i + ".raw.vcf")
      val filter_out = new File(base + ".F" + i + ".filtered.vcf")
      val eval_out = new File(base + ".F" + i + ".eval")

      add( snps(f, snp_out, nIntervals),
           filter(snp_out, filter_out),
           eval(filter_out, eval_out))

      f = f - downsamplingStep
      i = i + 1
    }
  }

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.intervals :+= targetIntervals
    this.jarFile = GATKjar
    this.reference_sequence = reference
    this.memoryLimit = 4
  }

  case class cov (fraction: Double, outFile: File) extends Percent20xCoverage with CommandLineGATKArgs {
    this.input_file :+= input
    this.out = outFile
    this.ndrs = true
    this.downsample_to_fraction = fraction
    this.jobName = queueLogDir + outFile + ".cov"
  }

  case class snps (fraction: Double, outFile: File, nIntervals: Int) extends UnifiedGenotyper with CommandLineGATKArgs {
    this.memoryLimit = 6
    this.downsample_to_coverage = 600
    this.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
    this.input_file :+= input
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.downsample_to_fraction = fraction
    this.scatterCount = nIntervals
    this.out = outFile
    this.analysisName = outFile + "_snps"
    this.jobName = queueLogDir + outFile
  }

  case class filter (inFile: File, outFile: File) extends VariantFiltration with CommandLineGATKArgs {
    this.filterName ++= List("SNPSBFilter","SNPQDFilter","SNPHRunFilter")
    this.filterExpression ++= List("\"SB>=0.10\"","\"QD<5.0\"","\"HRun>=4\"")
    this.clusterWindowSize = 10
    this.clusterSize = 3
    this.variantVCF = inFile
    this.out = outFile
    this.analysisName = outFile + "_filter"
    this.jobName = queueLogDir + outFile
  }

  // 3.) Variant Quality Score Recalibration - Generate Recalibration table
  case class VQSR(inFile: File, tranchesFiles: File, outFile: File) extends VariantRecalibrator with CommandLineGATKArgs {
    this.rodBind :+= RodBind("input", "VCF", inFile)
    this.rodBind :+= RodBind("hapmap", "VCF", hapmap, "known=false,training=true,truth=true,prior=15.0")
    this.rodBind :+= RodBind("omni", "VCF", omni, "known=false,training=true,truth=true,prior=12.0")
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP, "known=true,training=false,truth=false,prior=10.0")
    this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "HRun")
    this.tranches_file = tranchesFile
    this.recal_file = outFile
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
    this.analysisName = inFile + "_VQSR"
    this.jobName =  queueLogDir + outFile
  }

  // 4.) Apply the recalibration table to the appropriate tranches
  case class applyVQSR (inFile: File, tranchesFiles: File, outFile: File) extends ApplyRecalibration with CommandLineGATKArgs {
    this.rodBind :+= RodBind("input", "VCF", inFile)
    this.tranches_file = tranchesFile
    this.recal_file = inFile
    this.ts_filter_level = Some(trancheTarget)
    this.out = outFile
    this.analysisName = outFile + "_AVQSR"
    this.jobName =  queueLogDir + outFile
  }

  case class eval (inFile: File, outFile: File) extends VariantEval with CommandLineGATKArgs {
    this.noST = true
    this.noEV = true
    this.evalModule ++= List("TiTvVariantEvaluator", "CountVariants", "ValidationReport")
    this.stratificationModule ++= List("EvalRod", "CompRod", "Novelty")
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.rodBind :+= RodBind("eval", "VCF", inFile)
    this.rodBind :+= RodBind("comp", "VCF", fullCoverageVCF)
    this.out = outFile
    this.analysisName = outFile + "_VariantEval"
    this.jobName = queueLogDir + outFile
  }
}
