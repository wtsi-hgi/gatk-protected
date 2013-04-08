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
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.gatk.{VariantFiltration, CommandLineGATK, UnifiedGenotyper}
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine
import org.broadinstitute.sting.gatk.downsampling.DownsampleType

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

  @Input(doc = "a table with all the necessary information to process the data", fullName = "metadata", shortName = "m", required = false, exclusiveOf = "individual")
  var metaData: File = _

  /** ******************************************************************************
    * Additional Parameters that the pipeline should have pre-defined in the image
    * ******************************************************************************/

  @Argument(doc="Reference fasta file", fullName="reference", shortName="R", required=false)
  var reference: File = new File("/groups/reich/reference-genomes/human_hg19/human_g1k_v37/human_g1k_v37.fasta")

  @Argument(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=false)
  var dbSNP: Seq[File] = Seq(new File("/groups/reich/sw/gatk/GenomeAnalysisTK-2.3-9/bundle/dbsnp_137.b37.vcf"))

  @Argument(doc = "job queue for LSF", fullName = "queue", shortName = "queue", required = false)
  var queue: String = "short"

  @Argument(doc = "job queue for LSF", fullName = "project", shortName = "project", required = false)
  var project: String = "default"

  @Argument(doc = "tmp dir", fullName = "tmpDir", shortName = "tmpDir", required = false)
  var tmpDir: String = "/scratch/gd73/tmp/"

  @Argument(doc = "output dir", fullName = "outDir", shortName = "outDir", required = false)
  var outDir: String = "/scratch/gd73/"

  @Argument(doc = "job queue for LSF", fullName = "runname", shortName = "runname", required = false)
  var runname: String = "default"

  @Argument(doc = "input bam", fullName = "inputBAM", shortName = "I", required = false)
  var inputBAM:File = _

  @Argument(doc = "output vcf", fullName = "outputVCF", shortName = "o", required = false)
  var outputVCF:File = _

  @Argument(doc = "splitByContig", fullName = "splitByContig", shortName = "splitByContig", required = false)
  var splitByContig: Boolean = false

  @Argument(doc = "noPrior", fullName = "noPrior", shortName = "noPrior", required = false)
  var noPriors: Boolean = false


  @Input(doc = "Interval file with targets used in exome capture (used for QC metrics)", fullName = "targets", shortName = "targets", required = false)
  var targets: File = _

  @Argument(doc = "call indels as well as SNPs", fullName = "callIndels", shortName = "indels", required = false)
  var callIndels: Boolean = false

  @Argument(doc = "Default memory limit per job", fullName = "mem_limit", shortName = "mem", required = false)
  var memLimit: Int = 2

  @Argument(doc = "How many ways to scatter/gather", fullName = "scatter_gather", shortName = "sg", required = false)
  var nContigs: Int = 0


  @Argument(doc="Use Genotype Given Alleles mode for calling", shortName="doGGA", required=false)
  var doGGA: Boolean = false


  /** **************************************************************************
    * Main script
    * ***************************************************************************/

  val contigs:List[String] = List("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT")
  def script() {

    if (splitByContig) {

      for (chr <- contigs) {
        add(call(inputBAM, outputVCF,chr))
      }
    } else {
      add(call(inputBAM, outputVCF,""))
    }


  }

   // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK  {
    this.memoryLimit = memLimit
    this.isIntermediate = true
    this.jobQueue = queue
    this.reference_sequence = reference
  }

  case class filter(inVCF: File, outVCF: File) extends VariantFiltration with CommandLineGATKArgs {

  }
  case class call(inBAM: File, outVCF: File, chr:String) extends UnifiedGenotyper with CommandLineGATKArgs {
    this.input_file :+= inBAM
    this.isIntermediate = false
//    this.analysisName = outVCF + ".singleSampleCalling"
//    this.jobName = outVCF + ".singleSampleCalling"
    this.dbsnp = dbSNP(0)

    this.downsample_to_coverage = 600
    if (callIndels)
      this.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    else
      this.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP

     this.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES

    if (noPriors)
      this.noPrior = true
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
  }


}
