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

/**
 * Created with IntelliJ IDEA.
 * User: delangel
 * Date: 7/11/12
 * Time: 11:47 AM
 * To change this template use File | Settings | File Templates.
 */
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import collection.immutable.HashMap

class LargeScaleValidationPerPoolAnalysis extends QScript {
  qscript =>

  @Argument(doc="output path", shortName="outputDir", required=false)
  var outputDir: String =  "/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/"


  @Argument(doc="reference BAM file", shortName="refSample", required=false)
  var referenceSample: String = "NA12878"

  @Argument(doc="base output filename", shortName="runnName", required=false)
  var runName: String = "/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/finalPaperData/calls/finalRunGGANew"

  @Argument(doc="scatterCount", shortName="scatterCount", required=false)
  var variantCallerScatterCount: Int = 1

  @Argument(doc="ploidy", shortName="ploidy", required=false)
  var ploidy: Int = 24

  @Argument(doc="samplesPerPool", shortName="maxPool", required=false)
  var maxPool: Int = 93
  private val tmpDir: File = new File("/broad/hptmp/delangel/tmp/")
  private val reference: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  val pbase = "/humgen/gsa-hpprojects/dev/validationExperiments/largeScaleValidation/"

  val bamList: File = pbase+"dataAnalysis/bams_and_samples.list"

  val lof: File = new File(pbase+"dataAnalysis/assessment/ALL.wgs.phase1_release_v3.20101123.exomeSNPs.genotypes.atHC_LOFSites.vcf")
  val lofIntervals = pbase+"inputSets/LOF.DanielMacArthur_20120910.HighConfidence.interval_list"
  val allBaitIntervals = pbase+"baitDesign/ALL.wgs.allCombinedValidationSites.intervals_1bp.interval_list"

  val lofFile: File = new File(pbase+"outputVCFs/LOF.DanielMacArthur_20130212.fixed.vcf")
  val oneKGReleaseMinusBadPools: File = new File(pbase+"finalPaperData/ALL.wgs.phase1_release_v3.20101123.snps_indels_svs.genotypes.InGoodSamples.vcf")
  val exomeChip: File = new File(pbase+"finalPaperData/1000G.exomechip.20121009.subsetToGoodPoolSamples.genotypes.vcf")
  val axiomChip: File = new File(pbase+"finalPaperData/ALL.wex.axiom.20120206.snps_and_indels.subsetToGoodPoolSamples.genotypes.vcf")
  val omniChip: File = new File(pbase+"finalPaperData/Omni25_genotypes_2141_samples.b37.samplesInGoodPools.vcf")

  val sampleMapFile: File = new File(pbase+"finalPaperData/samples_and_pools.txt")
  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = qscript.reference
    this.memoryLimit = 2
    this.jobTempDir = qscript.tmpDir
    this.jobQueue = "hour"
  }

  case class SelectVariantsAtPool(pool:Int, inputFile: File, intervalss: File, sampleList: Seq[String]) extends SelectVariants with CommandLineGATKArgs {
    var poolStr:String = "%d".format(pool)
    if (poolStr.length == 1){
      poolStr = "0"+poolStr
    }

    this.variant = inputFile
    this.sn = sampleList
    this.ef = false
    this.ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES = true
    this.o = swapExt(qscript.outputDir + "/allPoolData/poolData"+ poolStr +"/", inputFile, ".vcf",".atPool"+poolStr+".vcf")
    this.intervals :+= intervalss
  }

  case class AnnotateWithAC(inputFile: File, acFile: File, intervalss: File, acName: String) extends VariantAnnotator with CommandLineGATKArgs {

    this.variant = inputFile

    this.o = swapExt(inputFile, ".vcf",".acAnnot.vcf")
    this.intervals :+= intervalss
    this.resource :+= acFile
    this.E :+= "resource."+acName
  }
  class VToti(inputFile: File, intervalss:File) extends VTot(inputFile) {
    this.intervals :+= intervalss
  }

  class VTot(inputFile: File) extends VariantsToTable with CommandLineGATKArgs {
    this.variant :+= inputFile
    this.allowMissingData = true
    this.raw = true
    this.o = swapExt(inputFile.getParentFile, inputFile, ".vcf",".table")

    this.F = Seq("CHROM","POS","REF","ALT","AC","DP","QUAL","FILTER","MLEAC","REFDEPTH","EVENTLENGTH","TYPE","resource.set","resource.LOF","NCALLED")
 //   this.GF = Seq("DP","MLPSAC")

  }
  def script = {
  // select omni genotypes for each pool

    // read bam file and assign to each pool number, taking care of skipped bad pools

    var poolFileArray = new Array[String](94)
    var poolInd:Int = 0

    var sampleMap = new HashMap[Int,Seq[String]]()

    for (line <- scala.io.Source.fromFile(sampleMapFile).getLines()) {
      // each line has format (sampleName pool#)
      val tokens:Array[String] = line.split("\t")
      val sampleName = tokens(0)
      val poolInd = Integer.parseInt(tokens(1))
      if (!(sampleMap contains poolInd))
        sampleMap += (poolInd -> Seq[String]())

      var list = sampleMap(poolInd)
      list :+= sampleName
      sampleMap += (poolInd -> list)
    }

    for (line <- scala.io.Source.fromFile(qscript.bamList).getLines()) {
      val tokens:Array[String] = line.split(" ")
      var bamFile = tokens(0)
      val sampleName = tokens(1)
      if (!sampleName.contentEquals(qscript.referenceSample) ) {
        poolInd =  Integer.parseInt(tokens(1).replace("1KG_Pool_",""))
//        println("%s %s".format(poolInd,bamFile)+line.toString)
        poolFileArray(poolInd) = bamFile

      }
    }
    // for each pool, call at omni/LOF sites
    // annotate with omni/LOF AC
    // variants to table
    // first get global non-pool tables
    add(new VToti(exomeChip,allBaitIntervals))
    add(new VToti(axiomChip,allBaitIntervals))
    add(new VToti(oneKGReleaseMinusBadPools,allBaitIntervals))
    add(new VToti(omniChip,allBaitIntervals))

    var doneA:Boolean = false
    for (pool <- 1 to qscript.maxPool ) {
      var poolStr:String = "%d".format(pool)
      if (poolStr.length == 1){
        poolStr = "0"+poolStr
      }

      if (poolFileArray(pool) != null) {

        val sampleNames:Seq[String] = sampleMap(pool)
        // OMNI, 1000G, axiom, exome chip
        val selExome = new SelectVariantsAtPool(pool,exomeChip, allBaitIntervals, sampleNames)
        add(selExome)
        val selExomeToTable = new VTot(selExome.o )
        add(selExomeToTable)


        // do the same for axiom
        val selAxiom = new SelectVariantsAtPool(pool,axiomChip, allBaitIntervals, sampleNames)
        add(selAxiom)
        val selAxiomToTable = new VTot(selAxiom.o )
        add(selAxiomToTable)

        // 1000 G
        val sel1000G = new SelectVariantsAtPool(pool,oneKGReleaseMinusBadPools, allBaitIntervals, sampleNames)
        add(sel1000G)
        val sel1000GToTable = new VTot(sel1000G.o )
        add(sel1000GToTable)

        // omni
        val selOmni = new SelectVariantsAtPool(pool,omniChip, allBaitIntervals, sampleNames)
        add(selOmni)
        val selOmniToTable = new VTot(selOmni.o )
        add(selOmniToTable)

        for (dataSet <- List("LOFSNP","LOFINDEL","exomeChip","omniPoly","millsPoly","omniMono","afSNPs","afIndels","unifSNPs","unifIndels")) {
          // Step 2: select genotypes from input file, for all datasets in which we have genotype data:
          // Step 1: select corresponding pool results from data
          var inputPoolFile = new File(runName+"."+dataSet + ".withRef.filtered.annotated.vcf")
          if (dataSet.contains("LOF"))
            inputPoolFile = new File(runName+"."+dataSet + ".withRef.filtered.annotated.LOF.vcf")
          var selPoolVar = new SelectVariantsAtPool(pool,inputPoolFile, allBaitIntervals, List("1KG Pool " + poolStr))
          add(selPoolVar)
          val selPoolToTable = new VTot(selPoolVar.o )
          add(selPoolToTable)

          if (!doneA) {
            // export non-pool data to table
            add(new VToti(inputPoolFile,allBaitIntervals))
          }
        }
        doneA = true

      }
    }

///       Same for LOF from 1KG
 /*
        var selVarLOF = new SelectVariants with CommandLineGATKArgs
        selVarLOF.variant = new File(qscript.lof)
        selVarLOF.sf :+= new File("/humgen/gsa-hpprojects/dev/delangel/LargeScaleValidation/poolData%d/samples.txt".format(pool))
        selVarLOF.ef = true
        selVarLOF.ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES = true
        selVarLOF.o = new File(qscript.outputDir + "/poolData"+ poolStr +"/ALL.wgs.LOFHCInPool"+ poolStr + ".genotypes.vcf")
        add(selVarLOF)

        var poolCallLOF = new SelectVariants with CommandLineGATKArgs
        poolCallLOF.V = new File("/humgen/gsa-hpprojects/dev/largeScaleValidation/dataAnalysis/results/pools."+qscript.runName+".LOF.filtered.vcf")
        poolCallLOF.sn :+= "1KG Pool " + poolStr
        poolCallLOF.o = qscript.outputDir + "/poolData"+ poolStr +"/ALL.wgs.PoolCaller_at_LOFHCInPool"+ poolStr + "."+qscript.runName+".snp.genotypes.vcf"
        add(poolCallLOF)

        var annotLOF = new VariantAnnotator with CommandLineGATKArgs
        annotLOF.variant = new File(poolCallLOF.out)
        annotLOF.out = qscript.outputDir + "/poolData"+ poolStr +"/ALL.wgs.PoolCaller_at_LOFHCInPool"+ poolStr +"."+qscript.runName+".annotated.snp.genotypes.vcf"
        annotLOF.resource :+= selVarLOF.o
        annotLOF.E :+= "resource.AC"
        annotLOF.intervalsString :+= qscript.lofIntervals
        add(annotLOF)

        var vtotLOF = new VariantsToTable with CommandLineGATKArgs
        vtotLOF.variant :+= new File (annotLOF.out)
        vtotLOF.allowMissingData = true
        vtotLOF.out = qscript.outputDir + "/poolData"+ poolStr +"/ALL.wgs.PoolCaller_at_LOFHCInPool"+ poolStr +"."+qscript.runName+".annotated.snp.genotypes.table"
        vtotLOF.F = Seq("CHROM","POS","AC","resource.AC","DP","QUAL","FILTER","MLEAC","REFDEPTH","EVENTLENGTH","TYPE")
        vtotLOF.GF = Seq("DP")
        add(vtotLOF)

    for (pool <- 1 to qscript.maxPool ) {
      if (!poolFileArray(pool-1).isEmpty) {
        val poolStr:String = "%d".format(pool)
        var samplesFile = new File(qscript.outputDir + "/poolData"+ poolStr +"/samples.txt")

        for (sample <- scala.io.Source.fromFile(samplesFile).getLines) {
          var selVar = new SelectVariants with CommandLineGATKArgs
    //      selVar.sn :+= sample
          selVar.env = true
          selVar.select :+= "AC==1"
          selVar.variant = new File(qscript.outputDir + "/poolData"+ poolStr +"/ALL.wgs.OmniInPool"+ poolStr +".correctedAll.genotypes.vcf")
          selVar.o = new File(qscript.outputDir + "/poolData"+ poolStr +"/ALL.wgs.OmniInPool"+ poolStr +".Singletons" +
            "."+qscript.runName+".genotypes.vcf")
          add(selVar)

          var selVar2 = new SelectVariants with CommandLineGATKArgs
          selVar2.sn :+= sample
          selVar2.env = true
          selVar2.variant = selVar.o
          selVar2.o = new File(qscript.outputDir + "/poolData"+ poolStr +"/ALL.wgs.OmniInPool"+ poolStr +".Singletons_sample_"+sample +
            "."+qscript.runName+".genotypes.vcf")
          add(selVar2)

          var poolCall = new UnifiedGenotyper with CommandLineGATKArgs
          poolCall.scatterCount = qscript.variantCallerScatterCount // the smallest interval list has 63 intervals, one for each Mb on chr20
          poolCall.input_file :+= new File(poolFileArray(pool-1))
          poolCall.input_file :+= qscript.refBam
          poolCall.ploidy = qscript.ploidy
          poolCall.out = qscript.outputDir + "/poolData"+ poolStr +"/ALL.wgs.PoolCaller_at_OmniInPool"+ poolStr +".Singletons_sample_"+sample +
          "."+qscript.runName+".snp.genotypes.vcf"
          poolCall.refsample = "gsa878"
          poolCall.referenceCalls = new File("/humgen/gsa-hpprojects/NA12878Collection/callsets/snps/NA12878.HiSeq.WGS.b37.recalibrated.99_5_cut_for_heng.vcf")
          poolCall.gt_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
          poolCall.out_mode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES
          poolCall.alleles = selVar2.o
          poolCall.ignoreLane = true
          poolCall.maxAltAlleles = Some(1)  // memory usage will overflow without this
          poolCall.dt = DownsampleType.NONE
          add(poolCall)

          var vtot = new VariantsToTable with CommandLineGATKArgs
          vtot.variant :+= new File (poolCall.out)
          vtot.allowMissingData = true
          vtot.out = qscript.outputDir + "/poolData"+ poolStr +"/ALL.wgs.PoolCaller_at_OmniInPool"+ poolStr +".Singletons_sample_"+sample +
            "."+qscript.runName+".annotated.snp.genotypes.table"
          vtot.F = Seq("CHROM","POS","AC","DP","QUAL","FILTER","MLEAC","REFDEPTH")
          vtot.GF = Seq("AD","GQ")
          add(vtot)
        }

      }
    }
    */
  }
}