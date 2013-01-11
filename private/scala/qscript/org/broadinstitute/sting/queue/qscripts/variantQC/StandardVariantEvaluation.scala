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

package org.broadinstitute.sting.queue.qscripts.variantQC

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk.RodBind
import org.broadinstitute.sting.queue.extensions.gatk._

class StandardVariantEvaluation extends QScript {
  // todo -- update to released version when things stabilize
  @Argument(doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("/home/radon01/depristo/dev/GenomeAnalysisTKFromLaptop/trunk/dist/GenomeAnalysisTK.jar")

  @Argument(shortName = "R", doc="B37 reference sequence: defaults to broad standard location", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "intervals", doc="intervals to evaluate.  Only supports evaluation on chromosome 20 now, as most evaluation data is there", required=false)
  val TARGET_INTERVAL: String = "20"

  @Argument(shortName = "includeUnion", doc="If provided, we'll create a union of the evaluation data sets for evaluation", required=false)
  val CREATE_UNION: Boolean = false

  @Argument(shortName = "dataDir", doc="Path to the standard evaluation data files", required=false)
  val DATA_DIR = "/humgen/gsa-hpprojects/GATK/data/Comparisons/StandardForEvaluation/b37/"

  @Argument(shortName = "evalStandard1000GCalls", doc="If provided, we'll include some standard 1000G data for evaluation", required=false)
  val EVAL_STANDARD_1000G_CALLS: Boolean = false

  val COMPS_DIR = DATA_DIR + "/comps/"
  val EVALS_DIR = DATA_DIR + "/evals/"

  @Argument(shortName = "moreSNPsToEval", doc="Path to additional SNP call sets for evaluation", required=false)
  val moreSNPsToEval: List[File] = Nil

  @Argument(shortName = "moreIndelsToEval", doc="Path to additional Indel call sets for evaluation", required=false)
  val moreIndelsToEval: List[File] = Nil

  val VARIANT_TYPES: List[String] = List("indels", "snps")
  val VARIANT_TYPE_VT: Map[String, List[org.broad.tribble.util.variantcontext.VariantContext.Type]] = Map(
    "indels" -> List(org.broad.tribble.util.variantcontext.VariantContext.Type.INDEL, org.broad.tribble.util.variantcontext.VariantContext.Type.MIXED, org.broad.tribble.util.variantcontext.VariantContext.Type.NO_VARIATION),
    "snps" -> List(org.broad.tribble.util.variantcontext.VariantContext.Type.SNP, org.broad.tribble.util.variantcontext.VariantContext.Type.NO_VARIATION)
  )

  val SITES_DIR: String = "sitesFiles"

  // path to b37 DBSNP
  @Argument(shortName = "dbsnp", doc="Path to DBSNP **VCF** for evaluation", required=false)
  val MY_DBSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf")
  //val MY_DBSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf");

  class Comp(val name: String, val evalType: String, val filename: String, val MakeHomVar: Boolean = false) {
    val originalFile = new File(COMPS_DIR + filename)
    val file: File = if ( MakeHomVar ) swapExt(originalFile, ".vcf",".homvar.vcf") else originalFile
    val sitesFile = new File(SITES_DIR + "/" + swapExt(file, ".vcf", ".sites.vcf").getName)
  }

  class Eval(val name: String, val evalType: String, val filename: String, val overrideFile: File = null ) {
    val file: File = if ( overrideFile != null ) overrideFile else new File(EVALS_DIR + "/" + filename)
  }

  var COMPS: List[Comp] = Nil
  def addComp(comp: Comp) { COMPS = comp :: COMPS }

  var EVALS: List[Eval] = Nil
  def addEval(eval: Eval) { EVALS = eval :: EVALS }
  def addEvalFromCMD(file: File, t: String) { addEval(new Eval(file.getName, t, null, file)) }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    this.intervalsString = List(TARGET_INTERVAL);
    this.reference_sequence = referenceFile;
    this.memoryLimit = 2
  }

  def initializeStandardDataFiles() = {
    //
    // Standard evaluation files for indels
    //
    addComp(new Comp("NA12878.homvar.GATK", "indels", "Indels.NA12878_WGS.filtered_Q50.0_QD5.0_SB-1.0_HR18.vcf", true))
    addComp(new Comp("CG.38samples", "indels", "CG.Indels.leftAligned.b37.vcf"))
    addComp(new Comp("NA12878.homvar.CG", "indels", "NA12878.CG.b37.indels.vcf", true))
    addComp(new Comp("g1k.pilot1.validation", "indels", "pilot1_indel_validation_2009.b37.vcf"))
    addComp(new Comp("NA12878.hand_curated", "indels", "NA12878.validated.curated.polymorphic.indels.vcf"))
    addComp(new Comp("NA12878.Mullikin", "indels", "NA12878.DIPline.NQScm.expanded.chr20.b37.minReads_2_or_gt2bp.vcf"))


    //
    // INDEL call sets
    //
    if ( EVAL_STANDARD_1000G_CALLS ) {
      addEval(new Eval("dindel", "indels", "20110208.chr20.dindel2.EUR.sites.vcf"))
      addEval(new Eval("si", "indels", "20101123.chr20.si.v2.EUR.sites.vcf"))
      addEval(new Eval("gatk", "indels", "EUR.phase1.chr20.broad.filtered.indels.sites.vcf"))
    }

    //
    // Standard evaluation files for SNPs
    //
    addComp(new Comp("NA12878.homvar.GATK", "snps", "NA12878.HiSeq19.cut.vcf", true))
    addComp(new Comp("CG.38samples", "snps", "CG.38samples.b37.vcf"))
    addComp(new Comp("NA12878.homvar.CG", "snps", "NA12878.CG.b37.snps.vcf", true))
    addComp(new Comp("HapMap3.3", "snps", "hapmap3.3.sites_r27_nr.b37_fwd.vcf"))
    addComp(new Comp("OMNI.2.5M", "snps", "omni2.5.1212samples.b37.sites.chr20.monoAreAC0.vcf"))
    addComp(new Comp("g1k.pilot1.validation", "snps", "1000G.snp.validation.b37.vcf"))

    //
    // SNP call sets
    //
    if ( EVAL_STANDARD_1000G_CALLS ) {
      addEval(new Eval("1000G.gatk.eurPlus.phase1", "snps", "EUR+.phase1.chr20.broad.recal.vrcut1p0.sites.vcf"))
      addEval(new Eval("1000G.high_specificity.phase1", "snps", "ALL.phase1.chr20.projectConsensus.highSpecificity.snps.genotypes.sites.vcf"))
    }
  }

  def script = {
    val sitesDir = new File(SITES_DIR)
    if ( ! sitesDir.exists ) sitesDir.mkdirs()

    initializeStandardDataFiles();

    // add additional files for evaluation, if necessary
    moreSNPsToEval.foreach(addEvalFromCMD(_, "snps"))
    moreIndelsToEval.foreach(addEvalFromCMD(_, "indels"))

    //
    // create hom-var versions of key files
    //
    for ( comp <- COMPS )
      if ( comp.MakeHomVar )
        add(new SelectHomVars(comp.originalFile, comp.file))

    for ( comp <- COMPS )
        add(new JustSites(comp.file, comp.sitesFile))

    //
    // Loop over evaluation types
    //
    for ( evalType <- VARIANT_TYPES ) {
      var evalsOfType = EVALS.filter(_.evalType == evalType)
      val compsOfType = COMPS.filter(_.evalType == evalType)

      if ( evalsOfType.size > 0 ) {

        // if desired and possible, create a union.X.vcf file
        if ( CREATE_UNION && evalsOfType.size > 1 ) {
          val union: File = new File("union.%s.vcf".format(evalType))
          add(new MyCombine(evalsOfType.map(_.file), union));
          evalsOfType = new Eval("union", evalType, null, union) :: evalsOfType
        }

        // our root VE
        val VE = new MyEval()
        VE.VT = VARIANT_TYPE_VT(evalType)
        VE.o = new File(evalType + ".eval")

        // add evals
        for ( calls <- evalsOfType )
          VE.rodBind :+= RodBind("eval_" + calls.name, "VCF", calls.file)

        // add comps
        //VE.rodBind :+= RodBind("dbsnp", "VCF", MY_DBSNP)
        for ( comp <- compsOfType )
          VE.rodBind :+= RodBind("comp_" + comp.name, "VCF", comp.sitesFile)

        add(VE)
      }
    }
  }

  /**
   * Select homozygous non-reference sites from a single deep data set
   */
  class SelectHomVars(@Input(doc="foo") vcf: File, @Output(doc="foo") out: File) extends SelectVariants with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("variant", "VCF", vcf)
    this.o = out
    this.select ++= List("\"AC == 2\"")
  }

  /**
   * A simple union
   */
  class MyCombine(@Input(doc="foo") vcfs: List[File], @Output(doc="foo") out: File) extends CombineVariants with UNIVERSAL_GATK_ARGS {
    for ( vcf <- vcfs )
      this.rodBind :+= RodBind(vcf.getName, "VCF", vcf)
    this.o = out
  }

  /**
   * A command line (cut) that removes all genotyping information from a file
   */
  class JustSites(@Input(doc="foo") in: File, @Output(doc="foo") out: File) extends CommandLineFunction {
    def commandLine = "cut -f 1-8 %s > %s".format(in, out)
  }

  /**
   * Base class for VariantEval used here
   */
  class MyEval() extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.noST = true
    this.evalModule :+= "ValidationReport"
  }
}

