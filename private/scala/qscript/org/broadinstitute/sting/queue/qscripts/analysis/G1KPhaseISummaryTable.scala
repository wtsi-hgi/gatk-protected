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

class G1KPhaseISummaryTable extends QScript {
  qscript =>

  @Argument(shortName = "L", fullName = "intervals", doc="intervals", required=false)
  val myIntervals: List[String] = null;

  @Argument(shortName = "nt", fullName = "nt", doc="Number of threads to use", required=false)
  val NumThreads: Int = 1;

//  val b37_decoy = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  val bundle = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/")
  val b37 = new File(bundle.getPath + "/human_g1k_v37.fasta")
  val dbSNP_b37 = new File(bundle.getPath + "/dbsnp_132.b37.vcf")

  val dbSNP_b37_129 = new File(bundle.getPath + "/dbsnp_132.b37.excluding_sites_after_129.vcf")
  val dbSNP_b37_129_with_pilot = new File("resources/dbsnp129_with_pilot.vcf")
  val pilotCalls = new File("resources/pilotSites.vcf")

  val populations = List("EUR", "ASN", "AFR", "AMR", "ALL")

  val callsets = Range(1,22).map("/humgen/1kg/releases/main_project_phaseI/ALL.chr%d.merged_beagle_mach.20101123.snps_indels_svs.genotypes.vcf.gz".format(_))

  val CCDS_BED = new File("resources/ucsc.ccds.bed")
  val CAPTURE_BED = new File("resources/20110225.exome.consensus.annotation.bed")
  val GENE_INTERVALS = Map("CCDS" -> CCDS_BED, "CAPTURE" -> CAPTURE_BED)

  def script = {
    for ( population <- populations ) {
      for ( (geneSetName, geneIntervals) <- GENE_INTERVALS )
        add(new evalVariants(population, geneSetName, geneIntervals))
    }
  }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    logging_level = "INFO";
    memoryLimit = 2;
    reference_sequence = b37
    intervalsString = myIntervals
  }

  // 5.) Variant Evaluation Base(OPTIONAL)
  class evalVariants(pop: String, geneSetName: String, geneIntervals: File) extends VariantEval with UNIVERSAL_GATK_ARGS {
    for ( callset <- callsets )
      this.eval :+= new File(callset)
    this.mergeEvals = true
    this.comp :+= new TaggedFile(dbSNP_b37_129, "dbsnp129")
    this.comp :+= new TaggedFile(pilotCalls, "pilot")
    this.comp :+= new TaggedFile(dbSNP_b37_129_with_pilot, "dbsnp129_and_pilot")
    this.sample = List("%s.samples.list".format(pop))
    this.out = new File("%s.samples.genes_%s.eval".format(pop, geneSetName))
    this.noEV = true
    this.EV = List("G1KPhaseITable")
    this.noST = true
    this.stratIntervals = geneIntervals
    this.ST = List("IntervalStratification")
    this.nt = NumThreads
  }
}
