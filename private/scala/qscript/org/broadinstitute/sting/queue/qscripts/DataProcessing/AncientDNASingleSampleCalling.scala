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

package org.broadinstitute.sting.queue.qscripts.DataProcessing

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.commandline.{ClassType, Hidden}
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine
import org.broadinstitute.sting.gatk.downsampling.DownsampleType
import org.broadinstitute.sting.queue.extensions.gatk.TaggedFile
import org.broadinstitute.sting.queue.extensions.gatk.UnifiedGenotyper
import org.broadinstitute.sting.queue.extensions.gatk.VariantFiltration
import org.broadinstitute.sting.queue.extensions.gatk.CommandLineGATK
import org.broadinstitute.sting.queue.extensions.gatk.VariantEval
import java.io.File

/**
 * Created with IntelliJ IDEA.
 * User: delangel
 * Date: 3/11/13
 * Time: 1:45 PM
 * To change this template use File | Settings | File Templates.
 */
class AncientDNASingleSampleCalling extends QScript{
  /** ***************************************************************************
    * Required Parameters
    * ***************************************************************************/

  @Input(doc = "input bam", fullName = "inputBAM", shortName = "I", required = false)
  var inputBAM:File = _

  /** ******************************************************************************
    * Additional Parameters that the pipeline should have pre-defined in the image
    * ******************************************************************************/

  @Argument(doc="Reference fasta file", fullName="reference", shortName="R", required=false)
  var reference: File = new File("/groups/reich/reference-genomes/human_hg19/human_g1k_v37/human_g1k_v37.fasta")

  @Argument(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=false)
  var dbSNP: Seq[File] = Seq(new File("/groups/reich/sw/gatk/GenomeAnalysisTK-2.3-9/bundle/dbsnp_137.b37.vcf"))

  @Argument(doc = "job queue for LSF", fullName = "queue", shortName = "queue", required = false)
  var queue: String = "short"

  @Argument(doc = "tmp dir", fullName = "tmpDir", shortName = "tmpDir", required = false)
  var tmpDir: String = "/scratch/gd73/tmp/"

  @Argument(doc = "output dir", fullName = "outDir", shortName = "outDir", required = false)
  var outDir: String = "/scratch/gd73/"

  @Argument(doc = "output vcf", fullName = "outputVCF", shortName = "o", required = true)
  var outputVCF:File = _

  @Argument(doc = "splitByContig", fullName = "splitByContig", shortName = "splitByContig", required = false)
  var splitByContig: Boolean = false


  @Argument(doc = "Interval file with targets used in exome capture (used for QC metrics)", fullName = "targets", shortName = "targets", required = false)
  var targets: File = _

  @Argument(doc = "call indels as well as SNPs", fullName = "callIndels", shortName = "indels", required = false)
  var callIndels: Boolean = false

  @Argument(doc = "Do hard filtering on raw calls", fullName = "doFiltering", shortName = "doFiltering", required = false)
  var doFiltering: Boolean = false

  @Argument(doc = "Do just one chr for testing", fullName = "doOneChr", shortName = "doOneChr", required = false)
  var doOneChr: Boolean = false

  @Argument(doc = "Only variants in this mask file are to be used. Only used if filtering enabled", fullName = "includeMask", shortName = "includeMask", required = false)
  var includeMask: File = _

  @Argument(doc = "Minimum MQ", fullName = "minMQ", shortName = "minMQ", required = false)
  var minMQ: Double = 20.0

  @Argument(doc = "Minimum QD", fullName = "minQD", shortName = "minQD", required = false)
  var minQD: Double = 1.0

  @Argument(doc = "Minimum DP", fullName = "minDP", shortName = "minDP", required = false)
  var minDP: Int = 5

  @Argument(doc = "Maximum DP", fullName = "maxDP", shortName = "maxDP", required = false)
  var maxDP: Int = 500

  @Argument(doc = "Minimum QUAL", fullName = "minQual", shortName = "minQual", required = false)
  var minQual: Double = 10.0

  @Argument(doc = "Maximum ReadPosRankSum", fullName = "minReadPosRankSum", shortName = "minRP", required = false)
  var minRP: Double = -2.0
  @Argument(doc = "Maximum HaplotypeScore", fullName = "maxHS", shortName = "maxHS", required = false)
  var maxHS: Double = 20.0


  @Argument(doc = "Default memory limit per job", fullName = "mem_limit", shortName = "mem", required = false)
  var memLimit: Int = 2

  @Argument(doc = "How many ways to scatter/gather", fullName = "scatter_gather", shortName = "sg", required = false)
  var nContigs: Int = 0

  @Argument(doc = "The path to the binary of tabix ", fullName = "tabix", shortName = "tabix", required = false)
  var tabixPath: File = new File("/groups/reich/sw/tabix-0.2.6/tabix")

  @Argument(doc = "The path to the binary of bgzip ", fullName = "bgzip", shortName = "bgzip", required = false)
  var bgzipPath: File = new File("/groups/reich/sw/tabix-0.2.6/bgzip")

  @Argument(doc="Only filter, no calling", shortName="onlyFilter", required=false)
  var onlyFilter: Boolean = false

  @Argument(doc="Compress outputs", shortName="compress", required=false)
  var compress: Boolean = false

  @ClassType(classOf[Double])
  @Argument(doc="Input AF priors", shortName="inputPriors", required=false)
  var inputPriors:Seq[Double] = Nil

  /** **************************************************************************
    * Main script
    * ***************************************************************************/
  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = memLimit
    this.isIntermediate = true
    this.jobQueue = queue
  }

  val contigs:List[String] = List("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT")
  val contigShort:List[String] = List("20")
  def script() {

    val ve = new Eval
    if (splitByContig) {

      for (chr <- {if(doOneChr==false)contigs else contigShort}) {
        ve.eval :+= callAndFilter(chr)
      }
    } else {
      ve.eval :+= callAndFilter("")
    }
    if (doOneChr)
      ve.intervalsString = contigShort
    ve.out =  outputVCF+".eval"
    add(ve)

  }

  def callAndFilter(chr:String) = {
    val caller = new call(inputBAM, outputVCF,chr)
    val bg = bgzip(caller.out)
    val indexer = index(bg.out)

    if (!onlyFilter)  {
      add(caller)
      if (compress) {
        add(bg)
        add(indexer)
      }
    }

    if (doFiltering) {
      val filt = filter({if (compress)bg.out else caller.out},indexer.out, chr)
//      val bgf = bgzip(filt.out)
      val idf = index(filt.out)
      add(filt)
      add(idf)
      filt.out

    } else {
      if (compress)
        bg.out
      else
        caller.out
    }

  }
   // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK  {
    this.memoryLimit = memLimit
 //   this.isIntermediate = true
    this.jobQueue = queue
    this.reference_sequence = reference
  }

  case class Eval() extends VariantEval with CommandLineGATKArgs {

//    this.eval :+= evalVCF
    this.dbsnp = dbSNP(0)
    this.doNotUseAllStandardModules = true
    this.evalModule = List("TiTvVariantEvaluator", "CountVariants", "CompOverlap")
    this.doNotUseAllStandardStratifications = true
    this.stratificationModule = Seq("Novelty", "Filter","AlleleCount")
    this.num_threads = 4
    this.memoryLimit = 8
  }

  case class filter(inVCF: File, inIdx: File, chr:String) extends VariantFiltration with CommandLineGATKArgs {
    @Input(doc = "input file") var inp = inVCF
    @Input(doc = "input file") var inpidx = inIdx

    this.isIntermediate = false
    this.variant = inVCF
    if (inVCF.endsWith(".vcf"))
      this.out = swapExt(outDir, inVCF, ".vcf", ".filtered.vcf")
    else if (inVCF.endsWith(".bcf"))
      this.out = swapExt(outDir,inVCF, ".bcf", ".filtered.bcf")
    else if (inVCF.endsWith(".vcf.gz"))
      this.out = swapExt(outDir,inVCF, ".vcf.gz", ".filtered.vcf.gz")

    @Output(doc = "output  file") var outp = this.out

    this.filterExpression = Seq("QUAL<"+minQual)
    this.filterName = Seq("LowQual")
    this.logging_level = "ERROR"
    this.filterExpression :+= "QD<"+minQD
    this.filterName :+= "LowQD"
    this.filterExpression :+=("MQ<"+minMQ)
    this.filterName :+=("LowMQ")
    this.filterExpression :+=("DP<"+minDP)
    this.filterName :+=("LowDP")
    this.filterExpression :+=("DP>"+maxDP)
    this.filterName :+=("HighDP")
    this.filterExpression :+= ("ReadPosRankSum<"+minRP)
    this.filterName :+= "LowReadPosRankSum"
    this.filterExpression :+= ("HaplotypeScore>"+maxHS)
    this.filterName :+= "HaplotypeScore"

    if (includeMask != null) {
      if (includeMask.toUpperCase.endsWith("BED")|| includeMask.toUpperCase.endsWith("BED.GZ")) {
        this.mask = new TaggedFile( includeMask, "BED" )
      } else this.mask = includeMask
      this.maskName = "InMask"
    }
    if (!chr.isEmpty)
      this.intervalsString :+= chr
  }
  case class call(inBAM: File, outVCF: File, chr:String) extends UnifiedGenotyper with CommandLineGATKArgs {
    this.input_file :+= inBAM
    this.isIntermediate = (doFiltering) // if filtering == true, caller output is intermediate
    this.dbsnp = dbSNP(0)

    this.downsample_to_coverage = 600
    if (callIndels)
      this.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    else
      this.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP

     this.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES

    this.stand_call_conf = Some(5.0)
    this.stand_emit_conf = Some(5.0) // will override with VF later

    if (!inputPriors.isEmpty) {
      this.inputPrior = inputPriors

    }

    if (splitByContig) {
      if (outVCF.endsWith(".vcf"))
        this.out = swapExt(outDir, outVCF, ".vcf", ".chr"+chr+".vcf")
      else if (outVCF.endsWith(".bcf"))
        this.out = swapExt(outDir,outVCF, ".bcf", ".chr"+chr+".bcf")
      else if (outVCF.endsWith(".vcf.gz"))
        this.out = swapExt(outDir,outVCF, ".vcf.gz", ".chr"+chr+".vcf.gz")


      if (chr.equals("MT"))
        this.dt = DownsampleType.NONE
      this.intervalsString :+= chr
      this.scatterCount = nContigs

    } else {
      this.out = outVCF
      this.scatterCount = nContigs
      if (targets != null)
        this.intervals :+= targets

    }
    // add useful annotations which are not default:
    this.A :+= "GCContent"
    this.A :+= "BaseCounts"
  }

  case class index(inVCF: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "input file") var inp = inVCF
    @Output(doc = "output  file") var out = swapExt(outDir,inVCF,"gz","gz.tbi")

    def commandLine = tabixPath + "  -f -p vcf -s 1 -b 2 -e 2 -c \\#  " + " " + inp

    this.memoryLimit = 2
    this.analysisName = out + ".tabix"
    this.jobName = out + ".tabix"

  }

  case class bgzip(inVCF: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "input file") var inp = inVCF
    @Output(doc = "output  file") var out = swapExt(outDir,inVCF,"vcf","vcf.gz")

    def commandLine = bgzipPath + " "+ inp

    this.memoryLimit = 2
    this.analysisName = out + ".bgzip"
    this.jobName = out + ".bgzip"

  }

}
