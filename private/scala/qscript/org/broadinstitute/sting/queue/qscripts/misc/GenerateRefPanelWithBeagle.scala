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
import scala.collection.JavaConversions._
import org.broadinstitute.sting.utils.text.XReadLines
import org.broadinstitute.sting.commandline.Argument
import java.io.PrintStream
import org.broadinstitute.sting.queue.extensions.gatk._
import net.sf.picard.reference.FastaSequenceIndex
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource
import collection.mutable.HashSet

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/20/12
 * Time: 11:36 AM
 * To change this template use File | Settings | File Templates.
 */

class GenerateRefPanelWithBeagle extends QScript {
  @Input(shortName="LP",fullName="lowPass",doc="VCF from low-pass calls and genotypes, containing likelihoods",required=true)
  var lowPassCalls : File = _
  @Input(shortName="EX",fullName="exome",doc="VCF from exome calls and genotypes, containing likelihoods",required=true)
  var exomeCalls : File = _
  @Input(shortName="CHIP",fullName="genotypeChip",doc="VCF from genotype chip genotypes, not necessarily containing likelihoods",required=true)
  var chipCalls : File = _
  @Input(shortName="CI",fullName="chunkIntervals",doc="Interval list containing the chunking strategy for running this analysis",required=true)
  var chunks : File = _
  @Input(shortName="R",fullName="referenceFile",doc="The reference fasta. Should be obvious.",required=false)
  var ref : File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  @Input(shortName="B",fullName="beagleJar",doc="Path to the beagle jar",required=false)
  var beagleJar = new File("/humgen/gsa-hpprojects/software/beagle/beagle.jar")
  @Input(shortName="PC",fullName="phasedComparison",doc="A file phased with beagle for comparison of phasing and genotypes. Usually lowpassed phased on its own.",required=true)
  var phasedComp : File = _
  @Input(shortName="X",fullName="excludeIdx",doc="Exclude these job indeces (useful for beagle failing when no variants present)",required=false)
  var exclude : File = _

  val BEAGLE_MEM_IN_GB : Int = 4
  val TMPDIR : String = System.getProperty("java.io.tmpdir");


  class RefineGenotypesWithBeagle(@Input beagleInput: File) extends CommandLineFunction {
    this.memoryLimit = BEAGLE_MEM_IN_GB

    // Note: These get set
    @Output val beaglePhasedFile: File = new File(beagleInput.getAbsolutePath +".phased.gz")
    @Output val beagleLikelihoods: File = new File(beagleInput.getAbsolutePath +".gprobs.gz")
    @Output val beagleRSquared: File = new File(beagleInput.getAbsolutePath +".r2")

    def commandLine = "java -Djava.io.tmpdir=%s -Xmx%dg -jar %s out=foo like=%s niterations=50 omitprefix=true".format(
      TMPDIR, BEAGLE_MEM_IN_GB, beagleJar,beagleInput.getAbsolutePath)
  }

  class Gunzip(@Input zippedFile : File ) extends CommandLineFunction {
    @Output val outputFile : File = new File(zippedFile.getAbsolutePath.replace(".gz",""))

    def commandLine = "gunzip -c %s > %s".format(zippedFile.getAbsolutePath,outputFile.getAbsolutePath)
  }

  def script = {

    // identify what to exclude
    var exIdx : HashSet[Int] = new HashSet[Int]
    if ( exclude != null ) {
      exIdx.addAll(asJavaCollection(asScalaIterable(new XReadLines(exclude)).map(u => u.toInt)))
    }
    // for each chunk
    var idx : Int = 0
    var panelChunks : List[File] = Nil
    asScalaIterator(new XReadLines(chunks)).foreach( chunk => {
      if ( ! exIdx.contains(idx) ) {
        // calculate the flanks (150kb)
        val chr = chunk.split(":")(0)
        val start = Integer.parseInt(chunk.split(":")(1).split("-")(0))
        val stop = Integer.parseInt(chunk.split(":")(1).split("-")(1))
        var flanks : List[String] = Nil
        var leftHalfFlank : String = null
        if ( start > 150000 ) {
          flanks :+= "%s:%d-%d".format(chr,start-150000,start)
          leftHalfFlank = "%s:%d-%d".format(chr,start-75000,start)
        }
        var rightHalfFlank : String = null
        val dif : Int = (new ReferenceDataSource(ref)).getReference.getSequenceDictionary.getSequence(chr).getSequenceLength - stop
        if ( dif > 0 ) {
          val ept = scala.math.min(dif,150000)
          flanks :+= "%s:%d-%d".format(chr,stop,stop+ept)
          rightHalfFlank = "%s:%d-%d".format(chr,stop,stop+ept/2)
        }
        // 1: multiply together the likelihoods
        val mLik = new MultiplyLikelihoods
        mLik.reference_sequence = ref
        mLik.Variants :+= lowPassCalls
        mLik.Variants :+= exomeCalls
        mLik.Variants :+= chipCalls
        mLik.intervalsString :+= chunk
        mLik.intervalsString ++= flanks
        mLik.out = new File(TMPDIR,"Mult_Likelihoods.chunk%d.vcf".format(idx))
        mLik.memoryLimit = 2
        add(mLik)

        // 2: create the input for beagle
        val beagOut = new ProduceBeagleInput
        beagOut.reference_sequence = ref
        beagOut.variant = mLik.out
        beagOut.out = new File(TMPDIR,"Mult_Likelihoods.chunk%d.beagle".format(idx))
        beagOut.intervalsString :+= chunk
        beagOut.intervalsString ++= flanks
        beagOut.memoryLimit = 2
        add(beagOut)

        // 3: refine the genotypes (and phase)
        val beagRefine = new RefineGenotypesWithBeagle(beagOut.out)
        add(beagRefine)

        // 3a: gunzip this stuff
        val phase = new Gunzip(beagRefine.beaglePhasedFile)
        add(phase)
        val like = new Gunzip(beagRefine.beagleLikelihoods)
        add(like)

        // 4: convert beagle output back to VCF
        val beag2vcf = new BeagleOutputToVCF
        beag2vcf.reference_sequence = ref
        beag2vcf.memoryLimit = Some(4)
        beag2vcf.out = new File(TMPDIR,"Ref_Panel.chunk%d.vcf".format(idx))
        beag2vcf.keep_monomorphic = true
        beag2vcf.beaglePhased = new TaggedFile(phase.outputFile,"ph,BEAGLE")
        beag2vcf.beagleProbs = new TaggedFile(like.outputFile,"pr,BEAGLE")
        beag2vcf.beagleR2 = new TaggedFile(beagRefine.beagleRSquared,"r2,BEAGLE")
        beag2vcf.intervalsString :+= chunk
        if ( leftHalfFlank != null ) {
          beag2vcf.intervalsString :+= leftHalfFlank
        }
        if ( rightHalfFlank != null ) {
          beag2vcf.intervalsString :+= rightHalfFlank
        }
        beag2vcf.variant = mLik.out
        panelChunks :+= beag2vcf.out
        beag2vcf.memoryLimit = 2
        add(beag2vcf)
      }
      idx = 1 + idx;
    })

    val gather = new MyVCFGather
    gather.vcfs ++= panelChunks
    gather.outVCF = new File("Reference_Panel_Final.vcf")
    add(gather)

    // now we want to compare this with a pre-phased file to ensure no truly gross errors
    val gtEval = new VariantEval
    gtEval.reference_sequence = ref
    gtEval.eval :+= new TaggedFile(gather.outVCF,"Consensus,VCF")
    gtEval.comp :+= new TaggedFile(phasedComp,"PhasedComp,VCF")
    gtEval.EV :+= "GenotypeConcordance"
    gtEval.EV :+= "GenotypePhasingEvaluator"
    gtEval.out = new File("Reference_Panel_Eval.gatk")
    gtEval.intervals :+= chunks
    gtEval.memoryLimit = 4
    add(gtEval)
  }

    class MyVCFGather extends InProcessFunction {
    @Input(doc="VCFs to be merged") var vcfs: List[File] = Nil
    @Output(doc="The final VCF to write to") var outVCF : File = _

    def run : Unit = {
      var stream : PrintStream = new PrintStream(outVCF)
      var first : XReadLines = new XReadLines(vcfs(0))
      var line : String = first.next()
      while ( line != null && line.startsWith("#") ) {
        stream.printf("%s%n",line)
        if ( ! first.hasNext ) {
          line = null
        } else {
          line = first.next()
        }
      }
      first.close()
      vcfs.map(
        u => asScalaIterator[String](new XReadLines(u)) ).foreach(
        x => x.filter( z => ! z.startsWith("#") ).foreach(
         s => stream.printf("%s%n",s)))
      stream.close()
    }
  }

}