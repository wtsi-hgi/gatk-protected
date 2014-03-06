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

package org.broadinstitute.sting.queue.qscripts.analysis


import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import java.io.FileWriter
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.utils.variantcontext.VariantContext
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantSummary

class G1KPhaseISummaryTable extends QScript {
  qscript =>

  @Argument(shortName = "L", fullName = "intervals", doc="intervals", required=false)
  val myIntervals: List[String] = null;

  @Argument(shortName = "nt", fullName = "nt", doc="Number of threads to use", required=false)
  val NumThreads: Int = 4;

  @Argument(shortName = "allPops", fullName = "allPops", doc="Run all populations, not just ALL", required=false)
  val allPops: Boolean = false;

//  val b37_decoy = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  val bundle = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/")
  val b37 = new File(bundle.getPath + "/human_g1k_v37.fasta")
  val dbSNP_b37 = new File(bundle.getPath + "/dbsnp_132.b37.vcf")

  val dbSNP_b37_129 = new File(bundle.getPath + "/dbsnp_132.b37.excluding_sites_after_129.vcf")
  val dbSNP_b37_135_minus_1000g = new File("resources/dbSNP_135.no1000GProduction.vcf")

  // CNV information
  val knownCNVsFile = new File("resources/known_deletions.bed")
  // an inclusive bed file that contains all human SVs from dbVAR classified as 'germline SVs'.
  val knownCNVsInclusive = new File("resources/dbvar.human.all.sets.GRCh37.ucsc.bed")
  // is a high precision bed file that contains all human SVs from dbVAR
  // that are classified as 'germline SVs' and are annotated with the
  // "Method"-tag "Sequencing" or "Sequence alignment" (http://www.ncbi.nlm.nih.gov/dbvar/studies/)
  // -- i.e., a bed file of SVs with basically nucleotide resolution breakpoint information.
  val knownCNVsPrecise = new File("resources/dbvar.human.sequencing.sets.GRCh37.ucsc.bed")

  val populations = List("EUR", "ASN", "AFR", "AMR", "ALL")

  val callsets = Range(1,23).map("/humgen/1kg/DCC/ftp/release/20110521/ALL.chr%d.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz".format(_))
  val X_callset = "/humgen/1kg/DCC/ftp/release/20110521/ALL.chrX.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"

  val CCDS_BED = new File("resources/ucsc.ccds.bed")
  //val CAPTURE_BED = new File("resources/20110225.exome.consensus.annotation.bed")
  val GENCODE_BED = new File("resources/gencode7.coding.bed")
  // we've converged on using GENCODE
  //val INTERVALS = Map("CCDS" -> CCDS_BED, "GENCODE" -> GENCODE_BED) // "CAPTURE" -> CAPTURE_BED,
  val INTERVALS = Map("GENCODE" -> GENCODE_BED) // "CAPTURE" -> CAPTURE_BED,

  def script = {
    for ( population <- if ( allPops ) populations else List("ALL") ) {
      for ( (cnvName, cnvFile) <- Map("inclusive" -> knownCNVsInclusive, "precise" -> knownCNVsPrecise) ) {
        for ( (geneSetName, geneIntervals) <- INTERVALS ) {
          add(new evalVariants(population, geneSetName, geneIntervals, cnvName, cnvFile, "autosome"))
          val evX = new evalVariants(population, geneSetName, geneIntervals, cnvName, cnvFile, "X")
          evX.eval = List(X_callset)
          add(evX)
        }
      }
    }
  }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    logging_level = "INFO";
    memoryLimit = 8;
    reference_sequence = b37
    intervalsString = myIntervals
  }

  // 5.) Variant Evaluation Base(OPTIONAL)
  class evalVariants(pop: String, geneSetName: String, geneIntervals: File, cnvName: String, cnvFile: File, callsetName: String) extends VariantEval with UNIVERSAL_GATK_ARGS {
    for ( callset <- callsets )
      this.eval :+= new File(callset)
    this.mergeEvals = true
    //this.comp :+= new TaggedFile(dbSNP_b37_129, "dbSNP_129")
    this.comp :+= new TaggedFile(dbSNP_b37_135_minus_1000g, "dbSNP_135_minus_1000g")
    //this.comp :+= new TaggedFile(dbSNP_b37, "dbSNP_132")
    this.sample = List("%s.samples.list".format(pop))
    this.out = new File("%s.samples.%s_calls.genes_%s.cnvs_%s.eval".format(pop, callsetName, geneSetName, cnvName))
    this.noEV = true
    this.EV = List("VariantSummary")
    this.noST = true
    this.keepAC0 = true
    this.stratIntervals = geneIntervals
    this.ST = List("IntervalStratification")
    this.nt = NumThreads
    this.knownCNVs = cnvFile
  }
}
