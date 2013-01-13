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

package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils.{GenotypeMergeType, VariantMergeType}
import org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript

class Onekg_table1 extends QScript {
  @Argument(doc="stage")
  var stage: String = _

  @Argument(doc="gatkJarFile")
  var gatkJarFile: File = _

  @Argument(shortName = "R", doc="gatkJarFile")
  var referenceFile: File = _

trait UNIVERSAL_GATK_ARGS extends CommandLineGATK { logging_level = "INFO"; jarFile = gatkJarFile; reference_sequence = referenceFile } // -L 1

class Target(project: String, snpVCF: String, indelVCF: String, calledGenome: Double, targetGenome: Double, pop: String, pilot : String, bam: String = null) {
    def reportFile: String = List(pop, pilot, "report").mkString(".")
    def extraArgs = { 
      val basic = "--project %s --snps %s --calledGenome %f --totalGenome %f --pop %s".format(project, snpVCF, calledGenome, targetGenome, pop)
      basic + (if ( indelVCF == null ) "" else " --indels " + indelVCF)
    }
    
    def getPilot = pilot
    def getProject = project
    def getPop = pop
    def getSNPVCF = snpVCF
    def getIndelVCF = indelVCF
    def hasIndelVCF = indelVCF != null 
    def getBAM = bam
    def hasBAM = bam != null
    def getDOC = List(getPilot, getPop, getProject, "doc").mkString(".")
    def getDOCSummaryFile = "doc/" + getDOC + ".sample_summary"
    def hasDOC = hasBAM
    private def getEval(t: String) = List(getPilot, getPop, getProject, t, "eval").mkString(".")
    def getSNPEval = getEval("snps")
    def getIndelEval = getEval("indels")
}

val RELEASE = "/humgen/1kg/DCC/ftp/release/2010_07/"

var targets: List[Target] = List()

val p1Targets = List(("CEU", 2.43e9), ("YRI", 2.39e9), ("CHBJPT", 2.41e9))

for ( (pop: String,called) <- p1Targets ) 
  targets ::= new Target("SRP000031", pop + ".pilot1.vcf", "v1/dindel-v2/"+pop+".low_coverage.2010_06.indel.genotypes.vcf", called, 2.85e9, pop, "pilot1")

// pilot 2
val p2Targets = List(("CEU", 2.264e9), ("YRI", 2.214e9))
for ( (pop: String, called) <- p2Targets )
  targets ::= new Target("SRP000032", RELEASE + "trio/snps/" + pop + ".trio.2010_03.genotypes.vcf.gz", "v1/dindel-v2/"+pop+".trio.2010_06.indel.genotypes.vcf", called, 2.85e9, pop, "pilot2")

// pilot 3
for (pop <- List("CEU", "CHB", "CHD", "JPT", "LWK", "TSI", "YRI")) {
  val indels = if ( pop != "LWK" ) "exon/indel/"+pop+".exon.2010_06.genotypes.vcf.gz" else null
  targets ::= new Target("SRP000033", "exon/snps/" + pop + ".exon.2010_03.genotypes.vcf.gz", indels, 1.43e6, 1.43e6, pop, "pilot3", "/humgen/gsa-hpprojects/1kg/1kg_pilot3/useTheseBamsForAnalysis/pilot3.%s.cleaned.bam".format(pop))
}

// merged files
targets ::= new Target("SRP000031", "pilot1.snps.merged.vcf", "pilot1.indels.merged.vcf", 2.42e9, 2.85e9, "all", "pilot1.merged")
targets ::= new Target("SRP000032", "pilot2.snps.merged.vcf", "pilot2.indels.merged.vcf", 2.565e9, 2.85e9, "all", "pilot2.merged")
targets ::= new Target("SRP000033", "pilot3.snps.merged.vcf", "pilot3.indels.merged.vcf", 1.43e7, 1.43e7, "all", "pilot3.merged")
targets ::= new Target("SRP00003.", "1kg.snps.merged.vcf", "1kg.indels.merged.vcf", 2.42e7, 2.85e9, "all", "1kg.merged")

val INTERVALS = Map(
    "pilot1" -> null,
    "pilot2" -> null,
    "pilot3" -> "/humgen/gsa-hpprojects/1kg/1kg_pilot3/documents/CenterSpecificTargetLists/results/p3overlap.targets.b36.interval_list"
    )

def script = stage match {
    case "ALL" =>
        // initial pilot1 merge -- autosomes + x
        for ( (pop: String,called) <- p1Targets ) {
            val auto = RELEASE + "low_coverage/snps/"+ pop +".low_coverage.2010_07.genotypes.vcf.gz"
            // todo -- remove fixed when Laura gives us the official calls
            val x = RELEASE + "low_coverage/snps/"+ pop +".low_coverage.2010_07.xchr.fixed.genotypes.vcf"
            val combineSNPs = new Combine(List(auto, x), pop + ".pilot1.vcf")
            add(combineSNPs)
        }

        // create pilot wide merges
        val pilots = List("pilot2", "pilot1", "pilot3") // order of perference in merging
        for ( pilot <- pilots ) {
            val pilotTargets = targets filter (_.getPilot == pilot)
            val combineSNPs = new Combine(pilotTargets.map(_.getSNPVCF), pilot + ".snps.merged.vcf")
            add(combineSNPs)
    
            if ( pilotTargets(0).getIndelVCF != null ) {
                val combineIndels = new Combine(pilotTargets.map(_.getIndelVCF).filter((x: String) => x != null), pilot + ".indels.merged.vcf")
                add(combineIndels)
            }
        }

        // create project wide merges
        val snps = "1kg.snps.merged.vcf"
        val indels = "1kg.indels.merged.vcf"

        //add(new Combine(pilots.map(_ + ".snps.merged.vcf"), snps))
        add(new Combine(pilots.map(_ + ".indels.merged.vcf"), indels))
	

    case "EVAL" =>
        // VariantEval of the SNPs
        for (target <- targets) {
          add(new VariantEval(target.getSNPVCF, target.getSNPEval))
          //add(new StatPop(target))
        }

    case "DOC" => 
        for (target <- targets) {
          if ( target.hasBAM ) 
            add(new DepthOfCoverage(target.getBAM, target.getDOC, INTERVALS(target.getPilot)))
        }
    case "MASK" => 
        for ( pop <- List("CEU", "YRI", "CHBJPT") ) 
            add(new MaskStats(pop))

    case _ => throw new Exception("Unknown stage" + stage)
}

// Using scala anonymous classes
class VariantEval(vcfIn: String, evalOut: String, vcfType: String = "VCF") extends org.broadinstitute.sting.queue.extensions.gatk.VariantEval with UNIVERSAL_GATK_ARGS {
    val vcfFile = new File(vcfIn)
    this.rodBind :+= RodBind("eval", vcfType, vcfFile)
    this.out = new File(evalOut)
    this.DBSNP = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b36.rod")
    this.reportType = VE2TemplateType.Grep
    this.noStandard = true;
    this.evalModule :+= "CompOverlap"
    this.memoryLimit = 3

    override def dotString = "VariantEval: " + vcfFile.getName
}

class StatPop(target: Target) extends CommandLineFunction {
    @Input(doc="foo") var snpVCF = new File(target.getSNPVCF)
    @Input(doc="foo") var snpEval = new File(target.getSNPEval)
    @Input(doc="foo", required=false) var indelVCF: File = if (target.hasIndelVCF) new File(target.getIndelVCF) else { null }
    @Output(doc="foo") var reportFile: File = new File(target.reportFile)
    override def dotString = "1kgStats: " + reportFile
    def commandLine = "python ~/dev/GenomeAnalysisTK/trunk/python/1kgStatsForCalls.py -v -a pilot_data.alignment.index -s pilot_data.sequence.index -r /broad/1KG/DCC/ftp/ -o " + target.reportFile + " " + target.extraArgs + (if (target.hasDOC) " -c " + target.getDOCSummaryFile else "") + " --snpsEval " + target.getSNPEval + (if (target.hasIndelVCF) " --indels " + target.getIndelVCF else "")
}

class Combine(vcfsInArg: List[String], vcfOutPath: String) extends org.broadinstitute.sting.queue.extensions.gatk.CombineVariants with UNIVERSAL_GATK_ARGS {
  val vcfs = vcfsInArg.map((x: String) => new File(x))
  val vcfFile = new File(vcfOutPath)
  this.variantmergeoption = VariantMergeType.UNION
  this.genotypemergeoption = GenotypeMergeType.PRIORITIZE
  this.out = vcfFile
  this.rodBind ++= vcfs.map( input => RodBind(input.getName,"VCF",input) )
  this.rod_priority_list = vcfs.map( _.getName ).mkString(",")
  override def dotString = "CombineVariants: " + vcfs.map(_.getName).mkString(",") + " => " + vcfFile.getName
}

class MaskStats(pop: String) extends CommandLineFunction {
    @Output(doc="foo") var outFile: File = new File(pop + ".stats")
    def commandLine = "python ~/dev/GenomeAnalysisTK/trunk/python/maskStats.py masks/" + pop + ".mask.fa.gz -x MT -x Y -o " + outFile
}

class DepthOfCoverage(bam: String, docOutPath: String, interval: String) extends org.broadinstitute.sting.queue.extensions.gatk.DepthOfCoverage with UNIVERSAL_GATK_ARGS {
  val bamFile = new File(bam)
  this.omitIntervalStatistics = true
  this.omitDepthOutputAtEachBase = true
  this.minBaseQuality = 0
  this.minMappingQuality = 0
  this.out = new File(docOutPath)
  this.input_file :+= bamFile
  if (interval != null) {
    this.intervalsString :+= interval
    this.excludeIntervalsString ++= List("MT", "Y")
  }

  override def dotString = "DOC: " + bamFile.getName
}
}
