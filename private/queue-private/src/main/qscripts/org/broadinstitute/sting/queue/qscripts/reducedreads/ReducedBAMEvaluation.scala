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

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.baq.BAQ

class ReducedBAMEvaluation extends QScript {
  @Argument(shortName = "R",       fullName = "reference", doc="ref", required=false) var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
  @Argument(shortName = "rbam",    fullName = "reduced_bam", doc="If you provide a reduced bam, reducing will be skipped.", required=false) var reducedBam: File = null;
  @Argument(shortName = "fbam",    fullName = "full_bam", doc="If you provide a reduced bam, reducing will be skipped.", required=false) var fullBam: File = null;
  @Argument(shortName = "ri",      fullName = "reduce_interval", doc="Interval to reduce at", required=false) var reduceInterval: String = null;
  @Argument(shortName = "cs",      fullName = "context_size", doc = "", required = false) protected var contextSize: Int = _
  @Argument(shortName = "minmap",  fullName = "minimum_mapping_quality", doc = "", required = false) protected var minMappingQuality: Int = _
  @Argument(shortName = "mintail", fullName = "minimum_tail_qualities", doc = "", required = false) protected var minTailQuality: Byte = _
  @Argument(shortName = "minvar",  fullName = "minimum_alt_proportion_to_trigger_variant", doc = "", required = false) protected var minAltProportionToTriggerVariant: Double = _
  @Argument(shortName = "mindel",  fullName = "minimum_del_proportion_to_trigger_variant", doc = "", required = false) protected var minIndelProportionToTriggerVariant: Double = _
  @Argument(shortName = "minqual", fullName = "minimum_base_quality_to_consider", doc = "", required = false) protected var minBaseQual: Int = _
  @Argument(shortName = "maxqual", fullName = "maximum_consensus_base_qual", doc = "", required = false) protected var maxQualCount: Byte = _
  @Argument(shortName = "ds",      fullName = "downsample_coverage", doc = "", required = false) protected var downsampleCoverage: Int = _

  val dbSNP_b37: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b37.leftAligned.vcf")
  val hapmap_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf"
  val omni_b37 = "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf"



  def script = {
    val sliceBAM =  swapExt(bam,".bam",".printreads.bam")
    val sliceVCF = swapExt(sliceBAM,".bam",".filtered.vcf")
    add(SliceBAM(bam, sliceBAM))
    callAndEvaluateBAM(sliceBAM, sliceVCF)

    for {
      a <- cs
      b <- adav
      c <- mbrc
    } yield {
    val header = ".CS-" + toStringLength(cs, a) + ".mrav-" + toStringLength(adav, b) + ".MBRC-" + toStringLength(mbrc, c)
    val reduceBAM = swapExt(bam, ".bam", header + ".reduced.bam")
    val reduceVCF = swapExt(reduceBAM,".bam",".filtered.vcf")
    val combineVCF = swapExt(reduceVCF, ".bam", ".filtered.combined.vcf")

    // Generate the new BAMs
    add(ReduceBAM(bam, reduceBAM, a, b, c ))

    // Call SNPs, filter and get variant eval report
    callAndEvaluateBAM(reduceBAM, reduceVCF)

    // Combine the callsets
    add(combine(reduceVCF, sliceVCF, combineVCF))

    // Get a report of the combined callsets
    val eval = new Eval(combineVCF) // evaluate the combined VCF
    eval.select = List("'set==\"Intersection\"'", "'set==\"fullBAM\"'", "'set==\"reducedBAM\"'", "'set==\"filterInreducedBAM-fullBAM\"'", "'set==\"reducedBAM-filterInfullBAM\"'")
    eval.selectName = List("Intersection", "fullBAM", "reducedBAM", "filterInreducedBAM-fullBAM", "reducedBAM-filterInfullBAM")
    add(eval)
    }

  }
  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = referenceFile;
    this.memoryLimit = 4
  }


  case class ReduceBAM(bam: File, outVCF: File, a: Int, b: Int, c: Int) extends ReduceReads with UNIVERSAL_GATK_ARGS with CoFoJa {
    this.memoryLimit = 3
    this.input_file = List(bam)
    this.o = outVCF
    this.CS = a   // Best safe value 7
    this.ADAV = b     //Best 30
    this.MBRC = c
    this.baq = BAQ.CalculationMode.OFF
    this.intervalsString = List(reduceInterval);
  }

  case class SliceBAM(bam: File, outVCF: File) extends PrintReads with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.input_file = List(bam)
    this.baq = BAQ.CalculationMode.CALCULATE_AS_NECESSARY
    this.o = outVCF
    this.intervalsString = List(reduceInterval);
  }

  case class combine (reducedBAMVCF: File, fullBAMVCF: File, outVCF: File) extends CombineVariants with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("fullBAM", "VCF", fullBAMVCF)
    this.rodBind :+= RodBind("reducedBAM", "VCF", reducedBAMVCF)
    this.rod_priority_list = "reducedBAM,fullBAM"
    this.filteredrecordsmergetype = org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED
    this.out = outVCF
  }

  def callAndEvaluateBAM(inBAM: File, outVCF: File) {
    val rawVCF   = swapExt(inBAM, ".bam", ".vcf")
    val recalVCF = swapExt(inBAM, ".bam", ".filtered.vcf")
    val recal    = swapExt(inBAM, ".bam", ".recal")
    val tranches = swapExt(inBAM, ".bam", ".tranches")

    add(Call(inBAM, rawVCF),
        HardFilter(rawVCF, recalVCF),
//        VQSR(rawVCF, recal, tranches),
//        applyVQSR(rawVCF, recal, tranches, recalVCF),
        Eval(recalVCF)
//        DiffableTable(rawVCF),   // for convenient diffing
//        DiffableTable(recalVCF)  // for convenient diffing
    )
  }

  case class Eval(@Input vcf: File) extends VariantEval with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind("eval", "VCF", vcf)
    if ( dbSNP.exists() )
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    this.doNotUseAllStandardStratifications = true
    this.doNotUseAllStandardModules = true
    this.evalModule = List("TiTvVariantEvaluator", "CountVariants")
    this.stratificationModule = List("EvalRod", "CompRod", "Novelty", "Filter", "JexlExpression")
    this.out = swapExt(vcf,".vcf",".eval")
    this.intervalsString = List(reduceInterval);
  }

  case class Call(inBAM: File, outVCF: File) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.input_file = List(inBAM)
    this.stand_call_conf = 50.0
    this.stand_emit_conf = 50.0
    this.dcov = DCOV;
    this.baq = BAQ.CalculationMode.OFF
    this.o = outVCF
    this.intervalsString = List(reduceInterval);

    if ( dbSNP.exists() )
      this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)

    if ( minimalVCF )
      this.group = List("none")

    if ( reduceInterval != null ) {
      this.intervalsString = List(reduceInterval)
    }
  }

  case class HardFilter (inVCF: File, outVCF: File) extends VariantFiltration with UNIVERSAL_GATK_ARGS {
    this.variantVCF = inVCF
    this.filterName = List("SNP_SB", "SNP_QD", "SNP_HRun")
    this.filterExpression = List("\"SB>=0.10\"", "\"QD<5.0\"", "\"HRun>=4\"")
    this.clusterWindowSize = 10
    this.clusterSize = 3
    this.out = outVCF
  }

  case class VQSR(inVCF: File, outRecal: File, outTranches: File) extends VariantRecalibrator with UNIVERSAL_GATK_ARGS {
    this.intervalsString = List(reduceInterval)
    this.rodBind :+= RodBind("input", "VCF", inVCF)
    this.rodBind :+= RodBind("hapmap", "VCF", hapmap_b37, "known=false,training=true,truth=true,prior=15.0")
    this.rodBind :+= RodBind("omni", "VCF", omni_b37, "known=false,training=true,truth=false,prior=12.0")
    this.rodBind :+= RodBind("dbsnp", "VCF", dbSNP_b37, "known=true,training=false,truth=false,prior=4.0")
    this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS", "DP")
    this.tranches_file = outTranches
    this.recal_file = outRecal
    this.allPoly = true
    this.tranche ++= List("100.0", "99.9", "99.5", "99.3", "99.0", "98.9", "98.8", "98.5", "98.4", "98.3", "98.2", "98.1", "98.0", "97.9", "97.8", "97.5", "97.0", "95.0", "90.0")
    this.analysisName = outTranches + "_VQSR"
    this.jobName =  outTranches + ".VQSR"
    this.mG = 5
    this.std = 14
    this.percentBad = 0.04
    this.nt = 5
}

  // 4.) Apply the recalibration table to the appropriate tranches
  case class applyVQSR (inVCF: File, inRecal: File, inTranches: File, outVCF: File) extends ApplyRecalibration with UNIVERSAL_GATK_ARGS {
    this.intervalsString = List(reduceInterval)
    this.rodBind :+= RodBind("input", "VCF", inVCF )
    this.tranches_file = inTranches
    this.recal_file = inRecal
    this.ts_filter_level = 99.0
    this.out = outVCF
    this.analysisName = outVCF + "_AVQSR"
    this.jobName =  outVCF + ".applyVQSR"
  }


  case class DiffableTable(@Input vcf: File) extends CommandLineFunction {
    @Output var out: File = swapExt(vcf,".vcf",".table")
    def commandLine = "cut -f 1,2,4,5,7 %s > %s".format(vcf, out)
  }
}

