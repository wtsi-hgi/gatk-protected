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

package org.broadinstitute.sting.queue.qscripts.misc.chartl

import org.broadinstitute.sting.queue.extensions.gatk.ByTranscriptEvaluator
import org.broadinstitute.sting.queue.library.ipf.vcf.VCFExtractIntervals
import org.broadinstitute.sting.queue.QScript

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 6/28/12
 * Time: 2:24 PM
 * To change this template use File | Settings | File Templates.
 */

class VariantEffectPipeline extends QScript {

  @Input(doc="The raw (un-annotated) VCF file",required=true,shortName="vcf",fullName="vcfToAnalyse")
  var rawInputVCF : File = _

  @Argument(doc="Project name",required=true,shortName="pname",fullName="projectName")
  var projectName : String = _

  @Argument(doc="Predict temporary dir",required=false,fullName="veptmpdir")
  var veptmpdir : String = "/broad/hptmp/chartl/vep/"

  override def script() {
    var ext : String = null
    if ( rawInputVCF.getAbsolutePath.endsWith(".vcf")) {
      ext = ".vcf"
    } else if ( rawInputVCF.getAbsolutePath.endsWith(".vcf.gz") ){
      ext = ".vcf.gz"
    } else {
      ext = ".bcf"
    }

    // compute the interval list
    var extract : VCFExtractIntervals = new VCFExtractIntervals(rawInputVCF)
    add(extract)
    // generate the predictions
    var predictRefseq : PVE = new PVE()
    predictRefseq.memoryLimit = Some(4)
    predictRefseq.variants = rawInputVCF
    predictRefseq.intervals = extract.listOut
    predictRefseq.tempdir = veptmpdir
    predictRefseq.annotVCF = swapExt(rawInputVCF,ext,".predicted.otherDB.vcf")
    predictRefseq.vCommand = "source /humgen/gsa-hphome1/chartl/projects/variantEffect/resources/runPredictorRefseqOnline.src | grep -v \"0+0k\" | sort -nk2,2 | perl /humgen/gsa-hphome1/chartl/sting/public/perl/sortByRef.pl - /humgen/1kg/reference/human_g1k_v37.fasta.fai | python /humgen/gsa-hphome1/chartl/projects/variantEffect/resources/headerToTop.py"
    //add(predictRefseq)

    var predictStandard :PVE = new PVE()
    predictStandard.memoryLimit=Some(4)
    predictStandard.variants = rawInputVCF
    predictStandard.intervals = extract.listOut
    predictStandard.tempdir = veptmpdir
    predictStandard.annotVCF = swapExt(rawInputVCF,ext,".predicted.standard.vcf")
    predictStandard.vCommand = "source /humgen/gsa-hphome1/chartl/projects/variantEffect/resources/runPredictor.src | grep -v \"0+0k\" | sort -nk2,2 | perl /humgen/gsa-hphome1/chartl/sting/public/perl/sortByRef.pl - /humgen/1kg/reference/human_g1k_v37.fasta.fai | python /humgen/gsa-hphome1/chartl/projects/variantEffect/resources/headerToTop.py"
    add(predictStandard)

    val stdEvals : List[ByTranscriptEvaluator] = generateEvals(predictStandard.annotVCF)
    val refEvals : List[ByTranscriptEvaluator] = generateEvals(predictRefseq.annotVCF)

    val stdParsed : List[ParseInfo] = stdEvals.map( u => {
      var pi = new ParseInfo
      pi.bteval = u.out
      pi.perGeneAggroFile = new File(u.out.getAbsolutePath+".genes.txt")
      pi.geneComposite = new File(u.out.getAbsolutePath+".worst.txt")
      pi
    })

    //addAll(stdParsed)

    val refParse : List[ParseInfo] = refEvals.map(u => {
      var pi = new ParseInfo
      pi.bteval = u.out
      pi.perGeneAggroFile = new File(u.out.getAbsolutePath+".genes.txt")
      pi.geneComposite = new File(u.out.getAbsolutePath+".worst.txt")
      pi
    })

    //addAll(refParse)

  }

  def generateEvals(inVCF : File ) : List[ByTranscriptEvaluator] = {
  // evaluate the VCF
    trait StdArgs extends ByTranscriptEvaluator {
      this.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
      this.eval = inVCF
      this.memoryLimit = Some(4)
    }

    var all_all : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    all_all.out = swapExt(inVCF,".vcf",".allTranscripts.allFrequency.eval")
    var all_01 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    all_01.out = swapExt(inVCF,".vcf",".allTranscripts.minAAF_0_01.eval")
    var all_05 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    all_05.out = swapExt(inVCF,".vcf",".allTranscripts.minAAF_0_05.eval")
    //add(all_all,all_01,all_05)
    var ccds_all : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    ccds_all.out = swapExt(inVCF,".vcf",".ccds.allFrequency.eval")
    var ccds_01 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    ccds_01.out = swapExt(inVCF,".vcf",".ccds.minAAF_0_01.eval")
    var ccds_05 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    ccds_05.out = swapExt(inVCF,".vcf",".ccds.minAAF_0_05.eval")
    //add(ccds_all,ccds_01,ccds_05)
    var refseq_all : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    refseq_all.out = swapExt(inVCF,".vcf",".refseq.allFrequency.eval")
    var refseq_01 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    refseq_01.out = swapExt(inVCF,".vcf",".refseq.minAAF_0_01.eval")
    var refseq_05 : ByTranscriptEvaluator = new ByTranscriptEvaluator with StdArgs
    refseq_05.out = swapExt(inVCF,".vcf",".refseq.minAAF_0_05.eval")
    //add(refseq_all,refseq_01,refseq_05)

    return List(all_all,all_01,all_05,ccds_all,ccds_01,ccds_05,refseq_all,refseq_01,refseq_05)
  }

  class PVE extends CommandLineFunction {
    @Input(doc="Queue jar file")
    var queueJar : File = new File("/humgen/gsa-hphome1/chartl/sting/dist/Queue.jar")

    @Input(doc="The input VCF file. Preferably sites-only.")
    var variants : File = _

    @Input(doc="Intervals file for the VCF")
    var intervals : File = _

    @Argument(doc="The number of variants to scatter into")
    var numVariantsPerScatter : Int = 10000

    @Argument(doc="temp dir")
    var tempdir : String = "./"

    @Argument(doc="javaTempDir")
    var javatmp = "/broad/hptmp/chartl/tmp/"

    @Output(doc="annotated vcf")
    var annotVCF : File = _

    @Argument(doc="The script for running Variant Effect Predictor. Must read from stdin and write to stdout.")
    var vCommand : String = _

    def commandLine : String = {
      val scriptPath : String = "/humgen/gsa-hphome1/chartl/sting/private/scala/qscript/org/broadinstitute/sting/queue/qscripts/misc/chartl/PredictVariantEffects.scala"
      "java -Xmx4g -Djava.io.tmpdir=%s -jar %s -S %s -V %s --numVariantsPerScatter %d -L %s --tempdir %s --vScript \"%s\" -out %s -run -bsub -jobQueue week".format(javatmp.getAbsolutePath,
        queueJar.getAbsolutePath,scriptPath,variants.getAbsolutePath,numVariantsPerScatter,
        intervals.getAbsolutePath,tempdir, vCommand, annotVCF.getAbsolutePath)
    }
  }

  class ParseInfo extends CommandLineFunction {
    @Input(doc="The eval file")
    var bteval : File = _

    @Output(doc="The per-gene aggregation")
    var perGeneAggroFile : File = _

    @Output(doc="The worst composite file")
    var geneComposite : File = _

    def commandLine : String = {
      "python /humgen/gsa-hphome1/chartl/projects/variantEffect/resources/runRefseq/merged %s %s %s".format(bteval.getAbsolutePath,
      perGeneAggroFile.getAbsolutePath,geneComposite.getAbsolutePath)
    }
  }


}

