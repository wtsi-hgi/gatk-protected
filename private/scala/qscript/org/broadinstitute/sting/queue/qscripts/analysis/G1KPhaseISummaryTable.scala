/*
 * Copyright (c) 2011, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.queue.qscripts.analysis

/*
 * Copyright (c) 2011, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

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
    this.stratIntervals = geneIntervals
    this.ST = List("IntervalStratification")
    this.nt = NumThreads
    this.knownCNVs = cnvFile
  }
}
