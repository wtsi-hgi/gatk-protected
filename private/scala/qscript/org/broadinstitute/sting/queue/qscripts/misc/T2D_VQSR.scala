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

package org.broadinstitute.sting.queue.qscripts.misc

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.text.XReadLines
import scala.collection.JavaConversions._
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.library.ipf.vcf.VCFExtractIntervals
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantQualityScore
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection
import org.broadinstitute.sting.gatk.walkers.genotyper.{GenotypeLikelihoodsCalculationModel, UnifiedGenotyperEngine}
import org.broadinstitute.sting.utils.baq.BAQ

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/18/12
 * Time: 9:36 PM
 * To change this template use File | Settings | File Templates.
 */

class T2D_VQSR extends QScript {

  @Input(shortName="BI",fullName="BroadCalls",required=true,doc="Broad unfiltered calls")
  var broadCalls : File = _

  @Input(shortName="U", fullName="Union",required=true,doc="U-Michigan Filtered SVM Union calls")
  var unionCalls : File = _

  @Input(shortName="I",fullName="Bams",required=true,doc="File listing all bam files for squaring off")
  var bamList : File = _

  @Input(shortName="NC",fullName="NumBamChunks",required=false,doc="Number of samples in each bam chunk")
  var nBamChunk : Int = 150

  val callDir : File = new File("/broad/shptmp/chartl/t2d/calling/")
  val bamDir : File = new File("/broad/shptmp/chartl/t2d/calling/bams/")
  val vqsrDir : File = new File("/broad/shptmp/chartl/t2d/vqsr/")
  val vqsrBase : File = new File(vqsrDir,"vqsr.base")

  val ref : File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  val hapMap : File = new TaggedFile("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf",
    "hapmap,vcf,known=false,training=true,truth=true,prior=15.0")
  val dbSNP : File = new TaggedFile("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.vcf",
    "dbsnp,vcf,known=true,training=false,truth=false,prior=8.0")
  // todo -- perhaps restrict ourselves only to EUR-polymorphic omni sites?
  val omni : File = new TaggedFile("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1856_samples.b37.vcf",
  "omni,vcf,known=false,training=true,truth=false,prior=12.0")

  val MAX_GAUSS = List(8,10,12,15,17,6,11,9)
  val DIRICH = List(0.0001,0.001,0.01)

  trait STAND_ARGS extends CommandLineGATK {
    this.reference_sequence = ref
    this.memoryLimit = Some(4)
  }

  def VQSR(inVCF: File, inList : File,  maxGauss: Int, dir: Double) : VariantRecalibrator = {
    val vqsr = new VariantRecalibrator with STAND_ARGS
    vqsr.an ++= List("QD","HaplotypeScore","MQRankSum","ReadPosRankSum","FS","MQ","InbreedingCoeff","DP")
    vqsr.mode = VariantRecalibratorArgumentCollection.Mode.SNP
    vqsr.input :+= inVCF
    vqsr.resource ++= List(hapMap,dbSNP,omni)
    vqsr.recalFile = swapExt(vqsrDir,vqsrBase, ".base", ".mg%d.dir%f.recal.txt".format(maxGauss,dir) )
    vqsr.tranchesFile = swapExt(vqsrDir,vqsrBase,".base",".mg%d.dir%f.tranche.txt".format(maxGauss,dir) )
    vqsr.intervals :+= inList
    vqsr.maxGaussians = Some(maxGauss)
    vqsr.dirichlet = Some(dir)
    vqsr.tranche ++= (8500 to 10000).map(u => (u.toDouble/100).toString)
    vqsr
  }

  def script = {

    // step 0: generate BI interval list
    val biCaInt = new VCFExtractIntervals(broadCalls,new File(callDir,"BI_calls.intervals.list"),true)
    add(biCaInt)

    // step 1: determine MI-unique calls by removing broad calls
    var remBroad = new SelectVariants with STAND_ARGS
    remBroad.variant = unionCalls
    remBroad.excludeIntervals :+= biCaInt.listOut
    remBroad.out = new File(callDir,"MI_Unique_calls.vcf")
    add(remBroad)

    // make sure to unfilter them
    var filt = new VariantFiltration with STAND_ARGS
    filt.variant = remBroad.out
    filt.out = new File(callDir,"MI_Unique_calls.noFilt.vcf")
    filt.invalidatePreviousFilters = true
    add(filt)

    val miUqInt = new VCFExtractIntervals(filt.out,new File(callDir,"MI_Unique_calls.intervals.list"),true)
    add(miUqInt)

    // step 2: merge the bam files into sub-chunks at MI-unique sites
    val printReads : List[PrintReads] = asScalaIterator(new XReadLines(bamList)).grouped(nBamChunk).zipWithIndex.map( (u : (Seq[String],Int)) => {
      val pr = new PrintReads with STAND_ARGS
      pr.input_file ++= u._1
      pr.out = new File(bamDir,"bam_chunk_%d.bam".format(u._2))
      pr.intervals :+= miUqInt.listOut
      pr.memoryLimit = Some(2)
      pr
    }).toList
    addAll(printReads)

    // step 3: recall using alleles only at MI unique sites
    var callMI = new UnifiedGenotyper with STAND_ARGS
    callMI.alleles = filt.out
    callMI.input_file ++= printReads.map( u => u.out )
    callMI.intervals :+= miUqInt.listOut
    callMI.genotyping_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
    callMI.out = new File(callDir,"BI_Calls_MI_Unique_Sites.vcf")
    callMI.stand_call_conf = Some(4.0)
    callMI.stand_emit_conf = Some(4.0)
    callMI.baq = BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    callMI.scatterCount = 50

    add(callMI)

    // step 4: combine initial broad calls with the unique ones
    var combineBroad = new CombineVariants with STAND_ARGS
    combineBroad.variant :+= broadCalls
    combineBroad.variant :+= callMI.out
    combineBroad.intervals :+= biCaInt.listOut
    combineBroad.intervals :+= miUqInt.listOut
    combineBroad.out = new File(callDir,"BI_Union_calls.vcf")
    combineBroad.scatterCount = 50

    add(combineBroad)

    // step 5: make an interval list
    var combInt = new VCFExtractIntervals(combineBroad.out,new File(vqsrDir,"Combined_Calls.intervals.list"),true)
    add(combInt)

    // step 6: run the VQSR for various selections of parameters
    for ( mg <- MAX_GAUSS ) {
      for ( dir <- DIRICH ) {
        add(VQSR(combineBroad.out,combInt.listOut,mg,dir))
      }
    }
  }
}