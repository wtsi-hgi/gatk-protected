/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.queue.qscripts.qc

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
import org.broadinstitute.gatk.utils.commandline.ClassType
import org.broadinstitute.gatk.utils.commandline

class CreateSyntheticContaminationScript extends QScript {
  qscript =>

  //@Input(shortName = "vcf", required = true, doc = "ContaminationVCF") var ContaminationVCF: File = _
  @Input(shortName = "b1", required = true, doc = "First BAM file") var bam1: File = _
  @Input(shortName = "b2", required = true, doc = "Second BAM file") var bam2: File = _
  //@Input var ContigFile: File = _
  @Input(shortName = "r", required = false, doc = "Reference sequence") var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")

  @ClassType(classOf[Double])
  @Argument(shortName = "f", doc = "list of fractions to contaminate with") var fractionList: List[Double] = _
  @Argument(required=false) var includeSwapped:Boolean=false

  @Argument(shortName = "pp", doc = "Path to the directory where all the picard Jars are stored", required=false) var picardPath: File = new File("/seq/software/picard/current/bin/picard.jar") 
  @Argument(required=false) var ContaminationJarPath:File=_

  @Argument(required=false) var bamSuffix:String="bam"


  trait CommonArguments extends CommandLineGATK {
    this.reference_sequence = referenceFile
    //this.intervalsString :+= ContaminationVCF
  }


  trait DocArguments extends DepthOfCoverage with CommonArguments {
    this.minBaseQuality = 17
    this.minMappingQuality = 30
    //this.omitLocusTable = true
    //this.omitBaseOutput = true
    //this.omitIntervals = true
    //this.omitDepthOutputAtEachBase = true
    //this.omitIntervalStatistics = true

  }

  class MyDepthOfCoverage extends DepthOfCoverage with DocArguments{
    @Output lazy val out_summary:File=swapExt(this.out,"",".sample_summary")
  }

  def script {

    if(ContaminationJarPath==null)ContaminationJarPath=picardPath

    assert(fractionList.reduce(math.max(_:Double,_:Double)) <= 1.0)
    assert(fractionList.reduce(math.min(_:Double,_:Double)) >= 0.0)


    for(swap:Boolean<-Set(false,includeSwapped)){
      doWork(swap)
    }

  }

  def doWork(swapBams:Boolean){

    val (innerbam1,innerbam2)=if(swapBams)(bam1,bam2) else (bam2,bam1)

    //printreads (using the contamination VCF, to save space)

    val pr1 = new PrintReads with CommonArguments
    pr1.input_file :+= innerbam1
    pr1.out = swapExt(innerbam1, bamSuffix, "small.bam")
    add(pr1)

    val pr2 = new PrintReads with CommonArguments
    pr2.input_file :+= innerbam2
    pr2.out = swapExt(innerbam2, bamSuffix, "small.bam")
    add(pr2)

    val doc1=new MyDepthOfCoverage with DocArguments
    doc1.input_file:+=pr1.out
    doc1.out=swapExt(pr1.out,".bam",".coverage")
    add(doc1)

    val doc2=new MyDepthOfCoverage with DocArguments
    doc2.input_file:+=pr2.out
    doc2.out=swapExt(pr2.out,".bam",".coverage")
    add(doc2)


    for (frac <- fractionList) {
      val oneMinusfrac=1-frac
      //downsample (using downsamplesam)

      val gnf=new GetNormalizedFractions()
      gnf.depthFile1=swapExt(doc1.out,"",".sample_summary")
      gnf.depthFile2=swapExt(doc2.out,"",".sample_summary")
      gnf.fraction=frac
      gnf.normalizedFractionAlpha=swapExt(doc1.out,"",s".$frac.fraction")
      gnf.normalizedFractionOneMinusAlpha=swapExt(doc2.out,"",s".$oneMinusfrac.fraction")
      add(gnf)

      val dss1 = new DownSampleSam()
      dss1.bam = pr1.out
      dss1.fraction=gnf.normalizedFractionAlpha
      dss1.out = swapExt(dss1.bam, "bam", frac.toString + ".bam")
      add(dss1)

      val dss2 = new DownSampleSam()
      dss2.fraction=gnf.normalizedFractionOneMinusAlpha
      dss2.bam = pr2.out
      dss2.out = swapExt(dss2.bam, "bam", oneMinusfrac.toString + ".bam")
      add(dss2)


      //extract headers
      val gh1 = new GetHeader()
      gh1.bam = dss1.bam
      gh1.out = swapExt(gh1.bam, ".bam", ".header")
      add(gh1)

      val gh2 = new GetHeader()
      gh2.bam = dss2.bam
      gh2.out = swapExt(gh2.bam, ".bam", ".header")
      add(gh2)

      //extract sample names

      val gs1 = new GetSample()
      gs1.bamHeader = gh1.out
      gs1.out = swapExt(gs1.bamHeader, "header", "sample")
      add(gs1)

      val gs2 = new GetSample()
      gs2.bamHeader = gh2.out
      gs2.out = swapExt(gs2.bamHeader, "header", "sample")
      add(gs2)

      //prepare header

      val rsih2 = new ReplaceSampleInHeader()
      rsih2.bamHeader = gs2.bamHeader
      rsih2.originalSample = gs2.out
      rsih2.newSample = gs1.out
      rsih2.out = swapExt(gs2.bamHeader, "header", "new.header")
      add(rsih2)

      //replace header  (replaceSAMHEader)

      val rh2 = new ReplaceHeader()
      rh2.bam = dss2.out
      rh2.header = rsih2.out
      rh2.out = swapExt(dss2.out, "bam", "swapped.sample.bam")
      add(rh2)

      //merge bams (mergeSamFiles)
      val ms1 = new MergeSamples()
      ms1.bam1 = dss1.out
      ms1.bam2 = rh2.out
      ms1.out = swapExt(dss1.out, "bam", "contaminated.with." + dss2.out)
      ms1.isIntermediate=false
      add(ms1)

     
     val bi1 = new BuildBamIndex()
     bi1.bam = ms1.out
     bi1.isIntermediate=false
     add(bi1)
	//estimate contamination

      //val ec1=new EstimateContamination()
      //ec1.bam=ms1.out
      //ec1.ContaminationJarPath=ContaminationJarPath
      //ec1.metrics=swapExt(ms1.out,"bam","metrics")
      //ec1.plotMetrics=swapExt(ms1.out,"bam","plot_metrics")
      //ec1.CentromereFile=ContigFile
      //add(ec1)
   }


    class GetNormalizedFractions extends CommandLineFunction{
      @Input var depthFile1:File=_
      @Input var depthFile2:File=_
      @Argument var fraction:Double=_
      @Output var normalizedFractionAlpha:File=_
      @Output var normalizedFractionOneMinusAlpha:File=_


      lazy private val rscript:String=s"""depthFile1="$depthFile1"
                                   |depthFile2="$depthFile2"
                                   |
                                   |alpha=$fraction
                                   |
                                   |resultFile1="$normalizedFractionAlpha"
                                   |resultFile2="$normalizedFractionOneMinusAlpha"
                                   |
                                   |Files=c(depthFile1,depthFile2)
                                   |Alpha=c(alpha,1-alpha)
                                   |Results=c(resultFile1,resultFile2)
                                   |
                                   |Depth=unlist(Map(function(x)read.table(x,skip=2)[["V3"]],c(Files)))
                                   |
                                   |A=(Alpha*Depth[2:1])/(Alpha[2:1]*Depth)
                                   |
                                   |Answer=rep(1,2)
                                   |Answer[which(A>=0&A<=1,)[1]]=A[which(A>=0&A<=1,)[1]]
                                   |
                                   |Map(function(x) write.table(Answer[x],Results[x],col.names=F,row.names=F),1:2)
                                   |""".stripMargin


      override def commandLine=required("echo", rscript) +
        required("|",escape=false) +
        required("R")+
        required("--vanilla")+
        required("--slave")

    }




    class DownSampleSam extends JavaCommandLineFunction {
      @Input var bam: File = _
      @Input var fraction: File = _

      @Output var out: File = _
      @Argument(required = false) var ip: Int = 50

      this.jarFile = new File(picardPath)


      override def commandLine = super.commandLine +
	required("DownsampleSam") +
        required("I=", bam, spaceSeparated = false) +
        required("O=", out,spaceSeparated = false) +
        required(s"P=$$(cat $fraction)",escape=false)

    }

    class GetHeader extends CommandLineFunction {
      @Input var bam: File = _
      @Output var out: File = _

      override def commandLine = required("samtools") + required("view") + required("-H", bam) + required(">",escape = false) + required(out)
    }

    class GetSample extends CommandLineFunction {
      @Input var bamHeader: File = _
      @Output var out: File = _

      override def commandLine = required("sed") +
        required("-n", "/SM:/{s/.*\\(SM:.*\\)\\t.*/\\1/p ;q};") +
        required(bamHeader) +
        required(">", escape = false) +
        required(out)

    }

    class ReplaceSampleInHeader extends CommandLineFunction {
      @Input var bamHeader: File = _
      @Input(doc = "file containing the sample_Alias1") var originalSample: File = _
      @Input(doc = "file containing the sample_Alias2") var newSample: File = _
      @Output var out: File = _

      override def commandLine = required("sed") +
        required("\'s/\'\"$(cat " + originalSample + ")\"\'/\'\"$(cat " + newSample + ")\"\'/\'", escape = false) +
        required(bamHeader)+
        required(">", escape = false) +
        required(out)
    }

    class ReplaceHeader extends JavaCommandLineFunction {
      @Input var bam: File = _
      @Input var header: File = _
      @Output var out: File = _

      this.jarFile = new File(picardPath)

      override def commandLine = super.commandLine +
        required("ReplaceSamHeader") +
	required("I=", bam,spaceSeparated = false) +
        required("O=", out,spaceSeparated = false) +
        required("HEADER=", header, spaceSeparated = false)
    }

    class MergeSamples extends JavaCommandLineFunction {
      @Input var bam1: File = _
      @Input var bam2: File = _
      @Output var out: File = _

      this.jarFile = new File(picardPath)

      override def commandLine = super.commandLine +
        required("MergeSamFiles") +
	required("I=", bam1,spaceSeparated = false) +
        required("I=", bam2,spaceSeparated = false) +
        required("O=", out,spaceSeparated = false)
    }

    class BuildIndex extends JavaCommandLineFunction {
     @Input var bam: File = _

     this.jarFile = new File(picardPath)

     override def commandLine = super.commandLine +
        required("BuildBamIndex") + 
        required("I=", bam, spaceSeparated = false)
    }

    class EstimateContamination extends JavaCommandLineFunction {
      @Input var bam: File = _
      @Output var metrics: File = _
      @Output var plotMetrics: File = _
      @Input var CentromereFile:File = _
      @Argument var ContaminationJarPath:File = _


      override def freezeFieldValues(){
        this.jarFile = new File(ContaminationJarPath, "CalculateContamination.jar")
        super.freezeFieldValues()
      }

      override def commandLine = {
        super.commandLine +
          //required("CONTAMINATION_VCF=", ContaminationVCF,spaceSeparated = false) +
          required("I=", bam,spaceSeparated = false) +
          required("CENTROMERE_FILE=", CentromereFile, spaceSeparated = false) +
          required("O=", metrics,spaceSeparated = false) +
          required("PLOT_METRICS_OUTPUT_FILE=", plotMetrics, spaceSeparated = false) +
          required("GRID_RESOLUTION=",0.001,spaceSeparated = false)  +
          required("CONFIDENCE_INTERVAL_LAMBDA=0.01")
      }

    }

  }
}



