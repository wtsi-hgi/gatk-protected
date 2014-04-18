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

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import scala.io.Source
import org.broadinstitute.sting.queue.function.RetryMemoryLimit

/**
 * Created with IntelliJ IDEA.
 * User: delangel
 * Date: 7/2/13
 * Time: 4:03 PM
 * To change this template use File | Settings | File Templates.
 */
class Phase3Staging extends QScript{
  @Argument(doc="the chromosome to process", shortName="minChrToProcess", required=false)
  var minChrToProcess: Int = 1
  @Argument(doc="the chromosome to process", shortName="maxChrToProcess", required=false)
  var maxChrToProcess: Int = 23
  @Argument(doc="queue", shortName="queue", required=false)
  var queue: String = "gsa"
  @Argument(doc="nct", shortName="nct", required=false)
  var nct: Int = 4
  @Argument(doc="chunkSize", shortName="chunkSize", required=false)
  var chunkSize: Int = 20000000
  @Argument(doc="bamsPerBlock", shortName="bamsPerBlock", required=false)
  var bamsPerBlock: Int = 5
  @Argument(doc="skipUG", shortName="skipUG", required=false)
  var skipUG: Boolean = false
  @Argument(doc="skipStaging", shortName="skipStaging", required=false)
  var skipStaging: Boolean = false

  @Argument(doc="pbase", shortName="pbase", required=false)
  var pbase: String = "/humgen/1kg/processing/production_wgs_final/"

  val analysisPanels = List("ASN","EUR","SAS","AFR","AMR")
  val reference = "/humgen/1kg/reference/human_g1k_v37_decoy.fasta"
  val chromosomeLength = List(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560)

  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = reference
    this.memoryLimit = Some(2)
    this.jobQueue = queue

  }

  case class PR1(listPiece:List[File], outFile:File, chunkString:String) extends PrintReads with CommandLineGATKArgs {
    @Output(doc = "chunk output", required = true) var outp: File = outFile
    @Input(doc = "sam files") var bams = listPiece

    this.nct = nct
    this.I = bams
    this.out = outp
    this.intervalsString :+= chunkString
    this.isIntermediate = true
    this.jobQueue = "hour"
  }

  case class PR2(listPiece:List[File], outFile:File, previousOut:File, chunkString:String) extends PrintReads with CommandLineGATKArgs {
    @Output(doc = "chunk output", required = true) var outp: File = outFile
    @Input(doc = "sam files") var bams = listPiece
    @Input(doc = "sam files") var dum = previousOut

    this.nct = nct
    this.I = bams
    this.out = outp
    this.intervalsString :+= chunkString
    this.isIntermediate = true
    this.jobQueue = "hour"

  }
  def blockStage(fileBunch: List[File], chrStr:String, start:Int, lim:Int, jobNumber:Int, ap:String) = {
    @Input(doc = "sam files") var bams = fileBunch

    val iter:Iterator[List[File]] = fileBunch grouped bamsPerBlock
    var idx:Int = 1
    val outDir = pbase + "chr" + chrStr + "/"  + "/" +  ap+"/"
    val stage2 =  new PrintReads with CommandLineGATKArgs
    val chunkString = "%s:%d-%d".format(chrStr, start, lim)
    var previousOut:File = fileBunch(1) //dummy
    var isFirst:Boolean = true
    while (iter.hasNext) {
      val listPiece:List[File] = iter.next()
      val outFile = new File(outDir+ap+"."+"chr"+chrStr+".chunk"+jobNumber.toString+".block"+idx.toString+".bam")

      // hacky way to create fake dependencies to avoid running the blocks in parallel.
      // Just a way to limit the concurrent IO processes into /broad/1kg.. ugh
      if (!skipStaging) {
        if (isFirst) {
          add(new PR1(listPiece, outFile,chunkString)  )
        }
        else {
          add(new PR2(listPiece, outFile, previousOut,chunkString))
        }
      }

      previousOut = outFile
      idx += 1
      stage2.I :+= outFile
      isFirst = false
    }
    // second part of staging: merge all samples into per-chunk analysis panel BAM
    stage2.nct = nct
    stage2.out = new File(outDir+ap+"."+"chr"+chrStr+".chunk"+jobNumber.toString+".bam")
    stage2.intervalsString :+= chunkString
    stage2.jobQueue = "week"
    if (!skipStaging) {
      add(stage2)
    }

    @Output(doc = "output bam file") var outbams = stage2.out

    stage2.out
  }


  def stageThisChunk(chrStr:String,start:Int, lim:Int, jobNumber:Int, chr:Int) = {

    /*
    1) For each of the continental groups, for current chunk:
    a) print sample BAM from thumper into staging area, for just this chunk
    b) Once chunk is in staging area, aggregate all samples in continental group
     */
//Console.out.println("chrStr:%s start:%d lim:%d jobN:%d".format(chrStr,start,lim,jobNumber) )

    var chunkBams = List[File]()
    for (cgroup <- analysisPanels) {
      val bamList = "/humgen/1kg/processing/production_wgs_final/phase3."+cgroup + ".20130502.low_coverage.alignment.bam.list"
      // read file and convert to list of files
      var fileBunch = List[File]()
      for (line: String <- Source.fromFile(bamList).getLines()) {
        fileBunch :+= line
      }

      chunkBams :+= blockStage(fileBunch, chrStr, start, lim, jobNumber,cgroup)

    }
    chunkBams
  }
  def script = {
    for(chr:Int <- minChrToProcess to maxChrToProcess) {
      var chrStr: String = chr.toString
      if(chr == 23) { chrStr = "X" }

      val lastBase: Int = chromosomeLength(chr-1)
      var start: Int = 1
      var stop: Int = start - 1 + chunkSize

      val commonPrefix = "/humgen/1kg/processing/production_wgs_final/chr"+chrStr+"/calls/ALL.chr"+chrStr

      val hcCombine = new CombineVariants with CommandLineGATKArgs
      hcCombine.assumeIdenticalSamples = true
      hcCombine.o = new File(commonPrefix + ".raw.HC.vcf")

      val ugCombine = new CombineVariants with CommandLineGATKArgs
      if (!skipUG) {
        ugCombine.assumeIdenticalSamples = true
        ugCombine.o = new File(commonPrefix + ".raw.UG.vcf")
      }

      var jobNumber: Int = 1
      do {
        val lim:Int = {
          if ( stop < lastBase )
            stop
          else
            lastBase
        }

        val chunkBams = stageThisChunk(chrStr,start, lim, jobNumber, chr)
        val outPrefix = commonPrefix + ".chunk"+jobNumber.toString
        val chunkstr = "%s:%d-%d".format(chrStr, start, lim)
        // call chunk
        val hcout = outPrefix + ".HC.vcf"
        add(new HC(chunkBams,hcout,chunkstr) )
        hcCombine.V :+= hcout
        if (!skipUG) {
          val ugout = outPrefix + ".UG.vcf"
          add(new UG(chunkBams,ugout,chunkstr) )
          ugCombine.V :+= ugout
        }
        // process loop
        start += chunkSize
        stop += chunkSize
        jobNumber += 1
      } while (start <= lastBase)
      if (!skipUG) {
        add(ugCombine)
      }
      add(hcCombine)


    }

  }

  case class HC( inputBamList: List[File], outFile: File, chunkStr:String ) extends HaplotypeCaller with RetryMemoryLimit with CommandLineGATKArgs {
    this.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
    this.intervalsString :+= chunkStr
    this.out = outFile
    this.input_file = inputBamList
    this.memoryLimit = 6
    this.scatterCount = 400
    this.stand_call_conf = 12
    this.stand_emit_conf = 12
    this.dcov = 200
    this.numPruningSamples = 3
    this.maxPathsPerSample = 5
    this.max_alternate_alleles = 3
//    this.isIntermediate = true
    this.jarFile = new File("/humgen/1kg/processing/production_wgs_final/dist/GenomeAnalysisTK.jar")
 //   @Output(doc = "output bam file") var outbams = outFile
  }

  case class UG( inputBamList: List[File], outFile: File,chunkStr:String ) extends UnifiedGenotyper with RetryMemoryLimit with CommandLineGATKArgs {
    this.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
    this.intervalsString :+= chunkStr
    this.out = outFile
    this.input_file = inputBamList
    this.memoryLimit = 3
    this.scatterCount = 10
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.stand_call_conf = 10
    this.stand_emit_conf = 10
    this.dcov = 200
    this.max_alternate_alleles = 3
//    this.isIntermediate = true
//@Output(doc = "output bam file") var outbams = outFile

  }

}