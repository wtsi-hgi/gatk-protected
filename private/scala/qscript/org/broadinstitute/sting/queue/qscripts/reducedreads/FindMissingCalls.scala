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

package org.broadinstitute.sting.queue.qscripts.reducedreads

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

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.Utils

class FindMissingCalls extends QScript {
  @Argument(shortName = "ref", doc = "Directory holding all of our data files", required=true)
  val referenceFile: File = null

  @Argument(shortName = "bam", doc = "", required=true)
  val BAM: File = null

  @Argument(shortName = "genotypes", doc = "", required=true)
  val GENOTYPES: File = null

  @Argument(shortName = "intervals", doc = "", required=true)
  val INTERVALS: File = null

  @Argument(shortName = "sample", doc = "", required=true)
  val SAMPLE: List[String] = Nil

  @Argument(shortName = "scatterCount", doc = "", required=false)
  val scatterCount: Int = 1

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = referenceFile;
    this.memoryLimit = 4
  }

  def script = {
    val subsites = swapExt(GENOTYPES, ".vcf", "." + SAMPLE.reduceLeft(_ + "." + _) + ".poly.vcf")
    val reducedBAM = swapExt(BAM, ".bam", ".reduced.bam")
    val fullVCF = swapExt(BAM, ".bam", ".vcf")
    val reducedVCF = swapExt(BAM, ".bam", ".reduced.vcf")

    val sv = new SelectVariants() with UNIVERSAL_GATK_ARGS
    sv.V = GENOTYPES
    sv.sample_name = SAMPLE
    sv.out = subsites
    sv.env = true
    sv.intervals :+= INTERVALS
    add(sv)

    // reduce
    val rr = new ReduceReads() with UNIVERSAL_GATK_ARGS
    rr.input_file :+= BAM
    rr.out = reducedBAM
    rr.intervals :+= INTERVALS
    rr.minqual = 20
    rr.scatterCount = scatterCount
    add(rr)

    add(call(BAM, subsites, fullVCF))
    add(call(reducedBAM, subsites, reducedVCF))

    val cv = findMissing(subsites, fullVCF, reducedVCF)
    val missingSitesFile = swapExt(cv.out, ".vcf", ".intervals")
    val missingSites = new MissingSitesList(cv.out, missingSitesFile)
    add(cv)
    add(missingSites)

    // reduce
    for ( thisBam: File <- List(BAM, reducedBAM) ) {
      val pr = new PrintReads() with UNIVERSAL_GATK_ARGS
      pr.input_file :+= thisBam
      pr.out = swapExt(thisBam, ".bam.", ".missing.bam")
      pr.intervals :+= missingSitesFile
      add(pr)
    }
  }

  def call(bam: File, sites: File, vcf: File) = {
    val ug = new UnifiedGenotyper() with UNIVERSAL_GATK_ARGS
    ug.input_file :+= bam
    ug.dbsnp = sites
    ug.out = vcf
    ug.BTI = "dbsnp"
    ug
  }

  def findMissing(allSites: File, fullSites: File, reducedSites: File) = {
    val cv = new CombineVariants() with UNIVERSAL_GATK_ARGS
    cv.V :+= new TaggedFile(allSites, "chip")
    cv.V :+= new TaggedFile(fullSites, "full")
    cv.V :+= new TaggedFile(reducedSites, "reduced")
    cv.out = swapExt(allSites, ".vcf", ".combined.vcf")
    cv.minimalVCF = true
    cv
  }

  class MissingSitesList(@Input combinedVCF: File, @Output outputFile: File) extends CommandLineFunction {
    def commandLine = "cat %s | grep -v Intersection | grep -v -i filter | grep -v \"#\" | grep full | awk '{print $1 \":\" $2}' > %s".format(combinedVCF, outputFile)
  }
}

