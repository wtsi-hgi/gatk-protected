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

package org.broadinstitute.sting.queue.qscripts.techdev

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel
import org.broadinstitute.sting.queue.util.QScriptUtils

/**
 * Created with IntelliJ IDEA.
 * User: carneiro
 * Date: 8/19/13
 * Time: 1:50 PM
 */

/**
 * this script downsamples an exome BAM several times and makes a coverage distribution
 * analysis (of bases that pass filters) as well as haplotype caller calls with a NA12878
 * Knowledge Base assessment.
 *
 * The R script "dwn_exome_analysis.r" produces plots of the loss of sensitivity/specificity
 * relative to coverage.
 *
 * This script was used for the "downsampling the exome" presentation
 */
class DownsampledExome extends QScript {
  qscript =>

  val b37 = new File ("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")
  val illumina  = new File ("/seq/references/HybSelOligos/whole_exome_illumina_coding_v1/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list")
  val agilent   = new File ("/seq/references/HybSelOligos/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list")


  @Argument(doc = "Input BAM file (only 1 BAM or list of BAMs). A list of bams will be expanded and each bam (in order) will match the interval list in order.",
    shortName = "I", required = false)
  var bams = Seq(new File ("/seq/picard_aggregation/DEV-2930/Pond-233795,_Spin_Column_Single_Index_Midi_Plate/current/Pond-233795,_Spin_Column_Single_Index_Midi_Plate.bam"),
    new File("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WEx.b37.NA12878.bam"))

  @Argument(doc = "Interval lists to use (one at a time) in order, matching the bams list. If less intervals than bams are provided the last interval in the list will be used to all remaining bams.",
    shortName = "L", required = false)
  var exome_intervals = Seq(illumina, agilent)

  @Argument(doc = "Expand the list of bams -- process each bam in the list individually",
    shortName = "eb", required = false)
  var expand_bams = false

  @Argument(doc = "Reference",
    shortName = "R", required = false)
  var reference = b37

  @Argument(doc = "List of samples to add to multi-sample calling",
    shortName = "M", required = false)
  var multiSample: File = null

  @Argument(doc = "number of iterations to downsample",
    shortName = "N", required = false)
  var niterations: Int = 50

  @Argument(doc = "do not run the HaplotypeCaller",
    shortName = "noHC", required = false)
  var no_hc = false

  @Argument(doc = "do not run the UnifiedGenotyper",
    shortName = "noUG", required = false)
  var no_ug = false

  @Argument(doc = "Scatter gather number of jobs",
    shortName = "sg", required = false)
  var scatter: Int = 25

  trait gatk_stuff extends CommandLineGATK {
    this.reference_sequence = reference
    this.memoryLimit = 4
  }

  def script() {

    // expand bam lists if any provided
    if (qscript.expand_bams) {
      var expanded_bams: Seq[File] = Seq()
      for (bam_list <- bams) {
        expanded_bams ++= QScriptUtils.createSeqFromFile(bam_list)
      }
      bams = expanded_bams
    }

    println("Processing " + bams.size + " bam files.")

    // process all the BAMs
    for (index <- 0 until bams.size) {
      val bam = bams(index)

      // make sure intervals and BAMs are matched in order. If interval list is smaller, repeat the last interval
      // (see description in the argument doc above)
      val interval = if (index < exome_intervals.size) {
        Seq(exome_intervals(index))
      } else {
        Seq(exome_intervals(exome_intervals.size-1))
      }

      val samples = QScriptUtils.getSamplesFromBAM(bam)
      val useKB = samples.size == 1 && samples("NA12878")

      val step = math.round(100 * 1.0/qscript.niterations).toDouble / 100
      for (ds <- step to 1.0 by step) {

        val fixed_ds = math.round(ds * 100) / 100.0

        val pr = new PrintReads with gatk_stuff
        pr.input_file :+= bam
        pr.dfrac = fixed_ds
        pr.out = swapExt(bam, ".bam", ".%.2f.bam".format(fixed_ds))
        pr.isIntermediate = true

        val bcd = new BaseCoverageDistribution with gatk_stuff
        bcd.input_file :+= pr.out
        bcd.intervals = interval
        bcd.out = swapExt(pr.out, ".bam", ".cov.grp")

        val cb = new CountBases with gatk_stuff
        cb.input_file :+= pr.out
        cb.jobOutputFile = swapExt(pr.out, ".bam", ".cb.out").toString

        add(pr, bcd, cb)

        val vcfs = make_calls(pr.out, interval)
        evaluate_precision(bam, interval, useKB, vcfs, pr.out, fixed_ds)
      }
    }
  }

  /**
   * Makes calls with UG and HC
   *
   * @param bam      the input bam file (downsampled)
   * @param interval interval file to make calls on
   * @return a tuple of vcfs. First the UG vcf, second the HC vcf
   */
  def make_calls(bam: File, interval: Seq[File]) : (File, File) = {
    val ug_out = if (no_ug) {
      null
    } else {
      val callMultiSample = qscript.multiSample != null
      val ug = new UnifiedGenotyper with gatk_stuff
      ug.out = swapExt(bam, ".bam", ".ug.vcf")
      ug.glm = GenotypeLikelihoodsCalculationModel.Model.BOTH
      ug.input_file :+= bam
      ug.intervals = interval
      ug.scatterCount = qscript.scatter

      val output = if (callMultiSample) {
        ug.input_file :+= qscript.multiSample
        val sv = new SelectVariants with gatk_stuff
        sv.V = ug.out
        sv.out = swapExt(ug.out, ".vcf", ".sv.vcf")
        sv.sample_name = Seq("NA12878")
        add(sv)
        sv.out
      } else {
        ug.out
      }
      add(ug)
      output
    }

    val hc_out = if (no_hc) {
      null
    } else {
      val hc = new HaplotypeCaller with gatk_stuff
      hc.out = swapExt(bam, ".bam", ".hc.vcf")
      hc.input_file :+= bam
      hc.intervals = interval
      hc.scatterCount = qscript.scatter
      hc.memoryLimit = 8
      add(hc)
      hc.out
    }
    (ug_out, hc_out)
  }

  def evaluate_precision(bam: File, interval: Seq[File], useKB: Boolean, vcfs: (File, File), pr: File, ds: Double) {
    // unwind the tuple into better names
    val ug = vcfs._1
    val hc = vcfs._2
    if(useKB) {
      val kb = new AssessNA12878 with gatk_stuff
      if (!no_ug)
        kb.V :+= new TaggedFile(ug, "UG")
      if (!no_hc)
        kb.V :+= new TaggedFile(hc, "HC")
      kb.out = swapExt(bam, ".bam", ".%.2f.kb.grp".format(ds))
      kb.BAM = pr
      kb.detailed = true
      kb.badSites = swapExt(bam, ".bam", ".%.2f.badsites.vcf".format(ds))
      kb.ip = 0
      kb.intervals = interval
      add(kb)
    } else {
      if(!no_ug) {
        val gc_ug = new GenotypeConcordance() with gatk_stuff
        gc_ug.comp = new TaggedFile(swapExt(bam, ".bam", ".%.2f.ug.vcf".format(ds)), "UG_full")
        gc_ug.eval = ug
        gc_ug.out = swapExt(bam, ".bam", ".%.2f.gc_ug.grp".format(ds))
        gc_ug.intervals = interval
        add(gc_ug)
      }
      if (!no_hc) {
        val gc_hc = new GenotypeConcordance() with gatk_stuff
        gc_hc.comp = new TaggedFile(swapExt(bam, ".bam", ".%.2f.hc.vcf".format(ds)), "UG_full")
        gc_hc.eval = ug
        gc_hc.out = swapExt(bam, ".bam", ".%.2f.gc_hc.grp".format(ds))
        gc_hc.intervals = interval
        add(gc_hc)
      }
    }
  }

}