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

package org.broadinstitute.sting.queue.qscripts.dev

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.QScript
import org.apache.commons.io.FilenameUtils;

class VQSR_parameterSearch extends QScript {
  qscript =>

  @Argument(shortName="gatk", doc="gatk jar file", required=true)
  var gatkJarFile: File = _

  @Argument(shortName="experiment", doc="experiment number", required=true)
  var experiment: String = "0000"

  @Argument(shortName="outputDir", doc="output directory", required=true)
  var outputDir: String = "./"

  @Argument(shortName="skipCalling", doc="If true, skip the calling part of the pipeline and only run VQSR on preset, gold standard VCF files", required=false)
  var skipCalling: Boolean = false

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK { logging_level = "INFO"; jarFile = gatkJarFile; memoryLimit = 2; }

  class Target(val baseName: String, val reference: File, val rodName: String, val bamList: File, val goldStandard_VCF: File, val intervals: String, val titvTarget: Double, val isLowpass: Boolean) {
    def name = qscript.outputDir + baseName
    def clusterFile = new File(name + ".clusters")
    def rawVCF = new File(name + ".raw.vcf")
    def filteredVCF = new File(name + ".filtered.vcf")
    def goldStandardName = qscript.outputDir + "goldStandard/" + baseName
    var goldStandardClusterFile: File = new File("")
    var gaussian: Int = 1
    var shrinkage: Double = 1.0
    var dirichlet: Double = 1.0
    var backoff: Double = 1.0
    var qualCutoff: Int = 1
    var std: Double = 1.0
    var useQD: Int = 1
    var useSB: Int = 1
    var useHS: Int = 1
    var useHRUN: Int = 1
    var useMQRST: Int = 1
    var useBQRST: Int = 1
    var useGC: Int = 1
    var useMQ: Int = 1
    var useSumGL: Int = 1
    var trainOmni: Int = 1
  }

  val hg18 = new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
  val b36 = new File("/humgen/1kg/reference/human_b36_both.fasta")
  val b37 = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  // ToDos:
  // reduce the scope of the datasets so the script is more nimble
  // figure out how to give names to all the Queue-LSF logs (other than Q-1931@node1434-24.out) so that it is easier to find logs for certain steps
  // create gold standard BAQ'd bam files, no reason to always do it on the fly

  // Analysis to add at the end of the script:
  // auto generation of the cluster plots
  // spike in NA12878 to the exomes and to the lowpass, analysis of how much of her variants are being recovered compared to single sample exome or HiSeq calls
  // produce Kiran's Venn plots based on comparison between new VCF and gold standard produced VCF

  // Define the target datasets here
  def lowPass = true
  val HiSeq = new Target("NA12878.HiSeq", hg18, "hg18", // BUGBUG: cut down to chr1
        new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.HiSeq.WGS.bwa.cleaned.recal.bam"),
        new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.vcf"),
        "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.hg18.intervals", 2.07, !lowPass)
  val WEx = new Target("NA12878.WEx", hg18, "hg18",
        new File("/humgen/gsa-hpprojects/NA12878Collection/bams/NA12878.WEx.cleaned.recal.bam"),
        new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.vcf"),
        "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.targets.interval_list", 2.6, !lowPass)
  val LowPassN60 = new Target("lowpass.N60", b36, "b36",                                                              // which reference the data is aligned to
        new File("/humgen/1kg/analysis/bamsForDataProcessingPapers/lowpass_b36/lowpass.chr20.cleaned.matefixed.bam"), // the bam list to call from
        new File("/home/radon01/depristo/work/oneOffProjects/VQSRCutByNRS/lowpass.N60.chr20.filtered.vcf"),           // the gold standard VCF file to run through the VQSR
        "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.b36.intervals", 2.3, lowPass)           // chunked interval list to use with Queue's scatter/gather functionality
  val LowPassAugust = new Target("ALL.august.v4", b37, "b37", // BUGBUG: kill this, it is too large
        new File("/humgen/1kg/processing/allPopulations_chr20_august_release.cleaned.merged.bams/ALL.cleaned.merged.list"),
        new File("/humgen/gsa-hpprojects/dev/data/AugChr20Calls_v4_3state/ALL.august.v4.chr20.filtered.vcf"),
        "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass)
  val LowPassEUR363Nov = new Target("EUR.nov2010", b37, "b37",
        new File("/humgen/1kg/processing/pipeline_test_bams/EUR.363sample.Nov2010.chr20.bam"),
        new File("/humgen/gsa-hpprojects/dev/rpoplin/haplotypeScore/sting_dev_oldQD_hs10/logs/EUR.nov.filtered.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
        "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass)
  val LowPassFIN79Nov = new Target("FIN.nov2010", b37, "b37",
        new File("/humgen/1kg/processing/pipeline_test_bams/FIN.79sample.Nov2010.chr20.bam"),
        new File("/broad/shptmp/rpoplin/pipeline_newHS7/FIN.nov2010.filtered.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
        "/humgen/1kg/processing/pipeline_test_bams/whole_genome_chunked.chr20.hg19.intervals", 2.3, lowPass)
  val TGPWExGdA = new Target("1000G.WEx.GdA", b37, "b37",
        new File("/humgen/1kg/processing/pipeline_test_bams/Barcoded_1000G_WEx_Reduced_Plate_1.cleaned.list"), // BUGBUG: reduce from 60 to 20 people
        new File("/humgen/gsa-scr1/delangel/NewUG/calls/AugustRelease.filtered_Q50_QD5.0_SB0.0.allSamples.SNPs_hg19.WEx_UG_newUG_MQC.vcf"), // ** THIS GOLD STANDARD NEEDS TO BE CORRECTED **
        "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list", 2.6, !lowPass)

  //val targets = List(HiSeq, WEx, LowPassN60, LowPassAugust, LowPassEUR363Nov, LowPassFIN79Nov, TGPWExGdA)
  val targets = List(LowPassEUR363Nov)
  def script = {
      def goldStandard = true

    var gaussianList = List(6)
    var shrinkageList = List(0.0001)
    var dirichletList = List(1000.0)
    var backoffList = List(1.3)
    var qualCutoffList = List(100)
    var stdList = List(4.5)
    var useQDList = List(1)
    var useSBList = List(1)
    var useHSList = List(1)
    var useHRUNList = List(1)
    var useMQRSTList = List(0)
    var useBQRSTList = List(0)
    var useGCList = List(0)
    var useMQList = List(0)
    var useSumGLList = List(0)
    var trainOmniList = List(1)

    if(experiment == "0000") {
      gaussianList = List(6,16)
      trainOmniList = List(0,1)
      useMQRSTList = List(0,1)
    }
    if(experiment == "0001") {
      gaussianList = List(6, 16)
      shrinkageList = List(0.0001, 0.01)
      dirichletList = List(0.001, 1000.0)
      backoffList = List(0.7, 1.0, 1.3)
      useQDList = List(0,1)
      useSBList = List(0,1)
      useHSList = List(0,1)
      useHRUNList = List(0,1)
      useMQRSTList = List(0,1)
      useBQRSTList = List(0,1)
      useSumGLList = List(0,1)
      trainOmniList = List(0,1)
    }
    if(experiment == "0002") {
      gaussianList = List(2, 10, 50)
      stdList = List(2.0, 4.5, 8.5)
      dirichletList = List(0.0001, 0.01)
      backoffList = List(0.5, 0.6, 0.9)
      useQDList = List(1)
      useSBList = List(0,1)
      useHSList = List(0,1)
      useHRUNList = List(0)
      useMQRSTList = List(0,1)
      useBQRSTList = List(0)
      useSumGLList = List(0,1)
      useGCList = List(0,1)
      useMQList = List(0,1)
      trainOmniList = List(0,1)
    }
    if(experiment == "0003") {
      qualCutoffList = List(5, 40, 100, 400)
      shrinkageList = List(0.0001, 0.001, 0.1)
      dirichletList = List(0.0001, 0.001, 0.01)
      useQDList = List(1)
      useSBList = List(0,1)
      useHSList = List(1)
      useHRUNList = List(0)
      useMQRSTList = List(0,1)
      useBQRSTList = List(0,1)
      useGCList = List(0,1)
      useMQList = List(0,1)
      useSumGLList = List(0,1)
      trainOmniList = List(0,1)
    }
    if(experiment == "0004") {
      gaussianList = List(5, 25)
      shrinkageList = List(0.01, 1.0, 100.0)
      dirichletList = List(0.001, 10.0, 1000.0)
      backoffList = List(0.6, 1.0, 1.4)
      useQDList = List(1)
      useSBList = List(1)
      useHSList = List(0,1)
      useHRUNList = List(0,1)
      useMQRSTList = List(0,1)
      useBQRSTList = List(0,1)
      useGCList = List(0,1)
      useMQList = List(0,1)
    }
    if(experiment == "0005") {
      gaussianList = List(4,50,100)
      shrinkageList = List(0.0001, 10.0)
      dirichletList = List(0.0001, 0.001)
      backoffList = List(0.2, 0.3, 0.6)
      stdList = List(0.5, 1.0, 10.0)
      useQDList = List(1)
      useSBList = List(1)
      useHSList = List(1)
      useHRUNList = List(0,1)
      useMQRSTList = List(0,1)
      useBQRSTList = List(0,1)
      useGCList = List(0,1)
      useMQList = List(0)
      trainOmniList = List(0,1)
    }



      for (target <- targets) {



        for(gaussian: Int <- gaussianList) {
          for(shrinkage: Double <- shrinkageList) {
            for(dirichlet: Double <- dirichletList) {
              for(backoff: Double <- backoffList) {
                for(qualCutoff: Int <- qualCutoffList) {
                  for(std: Double <- stdList) {
                  for(useQD: Int <- useQDList ) {
                  for(useSB: Int <- useSBList ) {
                  for(useHS: Int <- useHSList ) {
                  for(useHRUN: Int <- useHRUNList ) {
                  for(useMQRST: Int <- useMQRSTList ) {
                  for(useBQRST: Int <- useBQRSTList ) {
                  for(useGC: Int <- useGCList ) {
                  for(useMQ: Int <- useMQList ) {
                  for(useSumGL: Int <- useSumGLList ) {
                    for(trainOmni: Int <- trainOmniList) {

                      target.gaussian = gaussian
                      target.shrinkage = shrinkage
                      target.dirichlet = dirichlet
                      target.backoff = backoff
                      target.qualCutoff = qualCutoff
                      target.std = std
                      target.useQD = useQD
                      target.useSB = useSB
                      target.useHS = useHS
                      target.useHRUN = useHRUN
                      target.useMQRST = useMQRST
                      target.useBQRST = useBQRST
                      target.useGC = useGC
                      target.useMQ = useMQ
                      target.useSumGL = useSumGL
                      target.trainOmni = trainOmni
                      val clustersName: String = "%s_%d_%.4f_%.4f_%.1f_%d_%.1f_%d%d%d%d%d%d%d%d%d_%d.clusters".format(target.name, target.gaussian, target.shrinkage, target.dirichlet, target.backoff, target.qualCutoff, target.std, target.useQD, target.useSB, target.useHS, target.useHRUN, target.useMQRST, target.useBQRST, target.useGC, target.useMQ, target.useSumGL, target.trainOmni)
                      target.goldStandardClusterFile = new File(clustersName)
                      add(new GenerateVariantClusters(target, goldStandard))
                      add(new VariantRecalibratorTiTv(target, goldStandard))
                      add(new VariantRecalibratorNRS(target, goldStandard))
                    }
                  }
                  }
                  }
                  }
                  }
                  }
                  }
                  }
                  }
                }
              }
            }
          }
        }
      }
      }
  }

  def bai(bam: File) = new File(bam + ".bai")

  val FiltersToIgnore = List("DPFilter", "ABFilter", "ESPStandard", "QualByDepth", "StrandBias", "HomopolymerRun")

  // 3.) VQSR part1 Generate Gaussian clusters based on truth sites
  class GenerateVariantClusters(t: Target, goldStandard: Boolean) extends org.broadinstitute.sting.queue.extensions.gatk.GenerateVariantClusters with UNIVERSAL_GATK_ARGS {
      val name: String = if ( goldStandard ) { t.goldStandardName } else { t.name }
      this.reference_sequence = t.reference
      this.DBSNP = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_" + t.rodName + ".rod")
      this.rodBind :+= RodBind("hapmap", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.2/genotypes_r27_nr." + t.rodName + "_fwd.vcf")
      if(t.trainOmni == 0) {
        this.rodBind :+= RodBind("1kg", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg_pilot1_projectCalls/ALL.low_coverage.2010_07.hg19.vcf")
        this.rodBind :+= RodBind("truth", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.2/genotypes_r27_nr." + t.rodName + "_fwd.vcf")
      } else {
        this.rodBind :+= RodBind("1kg", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/764samples.deduped.b37.annot.vcf")
        this.rodBind :+= RodBind("truth", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/764samples.deduped.b37.annot.vcf")
      }
    this.rodBind :+= RodBind("input", "VCF", if ( goldStandard ) { t.goldStandard_VCF } else { t.filteredVCF } )
      this.clusterFile = if ( goldStandard ) { t.goldStandardClusterFile } else { t.clusterFile }
      //this.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
      if(t.useQD == 1) {
        this.use_annotation ++= List("QD")
      }
      if(t.useSB == 1) {
        this.use_annotation ++= List("SB")
      }
      if(t.useHS == 1) {
        this.use_annotation ++= List("HaplotypeScore1")
      }
      if(t.useHRUN == 1) {
        this.use_annotation ++= List("HRun")
      }
      if(t.useMQRST == 1) {
        this.use_annotation ++= List("MQRankSum")
      }
      if(t.useBQRST == 1) {
        this.use_annotation ++= List("BaseQRankSum")
      }
      if(t.useGC == 1) {
        this.use_annotation ++= List("GC")
      }
      if(t.useMQ == 1) {
        this.use_annotation ++= List("MQ")
      }
      if(t.useSumGL == 1) {
        this.use_annotation ++= List("sumGLbyD+")
      }
      if( t.useQD==0 && t.useSB==0 && t.useHS==0 && t.useHRUN==0 && t.useMQRST==0 && t.useBQRST==0 && t.useGC==0 && t.useMQ==0 && t.useSumGL==0) {
        this.use_annotation ++= List("MQ","QD","DP")
      }
      this.analysisName = name + "_GVC"
      this.intervalsString ++= List(t.intervals)
      this.qual = t.qualCutoff
      this.std = t.std
      this.mG = t.gaussian
      this.ignoreFilter ++= FiltersToIgnore
      this.dirichlet = t.dirichlet
      this.shrinkage = t.shrinkage
  }

  // 4.) VQSR part2 Calculate new LOD for all input SNPs by evaluating the Gaussian clusters
  class VariantRecalibratorBase(t: Target, goldStandard: Boolean) extends org.broadinstitute.sting.queue.extensions.gatk.VariantRecalibrator with UNIVERSAL_GATK_ARGS {
      val name: String = if ( goldStandard ) { t.goldStandardName } else { t.name }
      this.reference_sequence = t.reference
      this.DBSNP = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_" + t.rodName + ".rod")
      this.rodBind :+= RodBind("hapmap", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.2/genotypes_r27_nr." + t.rodName + "_fwd.vcf")
      if(t.trainOmni == 0) {
        this.rodBind :+= RodBind("1kg", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/1kg_pilot1_projectCalls/ALL.low_coverage.2010_07.hg19.vcf")
        this.rodBind :+= RodBind("truth", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.2/genotypes_r27_nr." + t.rodName + "_fwd.vcf")
      } else {
        this.rodBind :+= RodBind("1kg", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/764samples.deduped.b37.annot.vcf")
        this.rodBind :+= RodBind("truth", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/764samples.deduped.b37.annot.vcf")
      }
      this.rodBind :+= RodBind("input", "VCF", if ( goldStandard ) { t.goldStandard_VCF } else { t.filteredVCF } )
      this.clusterFile = if ( goldStandard ) { t.goldStandardClusterFile } else { t.clusterFile }
      this.analysisName = name + "_VR"
      this.intervalsString ++= List(t.intervals)
      this.ignoreFilter ++= FiltersToIgnore
      this.ignoreFilter ++= List("HARD_TO_VALIDATE")
      this.target_titv = t.titvTarget
      this.backOff = t.backoff
  }

  // 4a.) Choose VQSR tranches based on novel ti/tv
  class VariantRecalibratorTiTv(t: Target, goldStandard: Boolean) extends VariantRecalibratorBase(t, goldStandard) {
      this.tranche ++= List("1.0")
      this.out = new File("/dev/null")
      val tranchesName: String = "%s_%d_%.4f_%.4f_%.1f_%d_%.1f_%d%d%d%d%d%d%d%d%d_%d.titv.tranches".format(this.name, t.gaussian, t.shrinkage, t.dirichlet, t.backoff, t.qualCutoff, t.std, t.useQD, t.useSB, t.useHS, t.useHRUN, t.useMQRST, t.useBQRST, t.useGC, t.useMQ, t.useSumGL, t.trainOmni)
      this.tranchesFile = new File(tranchesName)
  }

  // 4b.) Choose VQSR tranches based on sensitivity to truth set
  class VariantRecalibratorNRS(t: Target, goldStandard: Boolean) extends VariantRecalibratorBase(t, goldStandard) {
      this.sm = org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibrator.SelectionMetricType.TRUTH_SENSITIVITY
      if(t.trainOmni == 0 ) {
        this.tranche ++= List("1.0")
      } else {
        this.tranche ++= List("2.5")        
      }
      this.out = new File("/dev/null")
      val tranchesName: String = "%s_%d_%.4f_%.4f_%.1f_%d_%.1f_%d%d%d%d%d%d%d%d%d_%d.ts.tranches".format(this.name, t.gaussian, t.shrinkage, t.dirichlet, t.backoff, t.qualCutoff, t.std, t.useQD, t.useSB, t.useHS, t.useHRUN, t.useMQRST, t.useBQRST, t.useGC, t.useMQ, t.useSumGL, t.trainOmni)
      this.tranchesFile = new File(tranchesName)
  }
}
