/**
 * Quick use of GCContentByInterval walker using an unmerged list of intervals
 * 
 * Takes a list of intervals and calculates the GC Content of each interval
 * separately so the GATK engine doesn't merge the intervals.
 * 
 * @author carneiro
 * @since 1/2/12
 */


import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import io.Source._


class IntervalGCContent extends QScript {
  @Argument(shortName = "i", required=true, doc = "Intervals file")      var intervalsFile: File = _
  @Argument(shortName = "r", required=false, doc = "Reference sequence") var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  @Argument(shortName = "o", required=false, doc = "Output file")        var outputFile: File = new File("interval.out")

  def script {

    var outputList: List[File] = List()
    var index: Int = 1

    for (interval <- fromFile(intervalsFile).getLines()) {
      if (!interval.startsWith("@")) {
        val output = new File("tmp-" + index + ".out")

        val walker = new GCContentByInterval()
        walker.reference_sequence = referenceFile
        walker.intervalsString = List(interval)
        walker.out = output
        walker.memoryLimit = 2
        add(walker)

        outputList :+= output;
        index = index + 1;
        
        if (index % 100 == 0)
          logger.warn("Created " + index + " jobs...")
      }
    }   
    add(cat(outputList, outputFile))
  }

  case class cat (files: List[File], output: File) extends CommandLineFunction {
    var cline: String = "cat "
    for (file <- files) {
      cline += file + " "
    }
    def commandLine = cline + " > " + output
  }

}
