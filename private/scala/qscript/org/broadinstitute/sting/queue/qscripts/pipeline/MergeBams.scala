package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.function.{RetryMemoryLimit, ListWriterFunction}
import org.broadinstitute.sting.queue.extensions.gatk.{CommandLineGATK, PrintReads}

class MergeBams extends QScript {

  @Input(doc="BAM list files. Name must be <projectName>.bam.list", shortName="I")
  var bamList: File = _

  @Input(shortName="R", doc="The reference file for the bam files.")
  var reference: File = _ // _ is scala shorthand for null

  @Input(shortName="L", doc="An optional file with a list of intervals to proccess.",  required=false)
  var intervalsFile: File = _

  @Argument(shortName="fileSize", doc="set the number of the bams in each merged file", required=false)
  var fileSize: Int = 1000

  @Argument(shortName="mergedDir", doc="merge directory", required=false)
  var mergedBamDir: String = "./mergedBams/"

  @Argument(doc="chrs", shortName="C", fullName = "chr", required=false)
  var chrs: Seq[String] = Nil

  def script(){
    val expandIntervals = 50
    val mergeBamMemoryLimit = 16

    trait CommandLineGATKArgs extends CommandLineGATK {
      this.reference_sequence = reference
      this.intervals = Seq(intervalsFile)
      this.interval_padding = expandIntervals
    }

    val bams = io.Source.fromFile(bamList).getLines().toSeq
    val bamGroups = bams.grouped(fileSize).toSeq

    var bamGroupNumber = 0
    var mergeBamLists = Seq.empty[File]

    for (bamGroup <- bamGroups) {
      val mergeBamList = new ListWriterFunction
      mergeBamList.inputFiles = bamGroup
      mergeBamList.listFile = mergedBamDir + "bamLists/%03d.bam.list".format(bamGroupNumber)
      add(mergeBamList)

      mergeBamLists :+= mergeBamList.listFile
      bamGroupNumber += 1
    }

    if (chrs != Nil){
      for (chr <- chrs) {

        val chrDir = "chrs/" + chr + "/"

        trait ChromosomeIntervals extends CommandLineGATKArgs{
          this.intervalsString :+= chr
          this.interval_set_rule = org.broadinstitute.sting.utils.interval.IntervalSetRule.INTERSECTION
        }

        bamGroupNumber = 0
        var chrBams = Seq.empty[File]
        for (bamGroup <- bamGroups) {
          val mergeBam = new PrintReads with ChromosomeIntervals
          mergeBam.input_file = Seq(mergeBamLists(bamGroupNumber))
          mergeBam.out = mergedBamDir + chrDir + "chr%s.%03d.bam".format(chr, bamGroupNumber)
          mergeBam.memoryLimit = mergeBamMemoryLimit
          mergeBam.nct = 4 //todo is it make sense to use it with such an high memoryLimit?
          add(mergeBam)

          chrBams :+= mergeBam.out
          bamGroupNumber += 1
        }


        val chrMergeBamList = new ListWriterFunction
        chrMergeBamList.inputFiles = chrBams
        chrMergeBamList.listFile = mergedBamDir + chrDir + "chr%s.merged.bam.list".format(chr)
        add(chrMergeBamList)
      }
    }

    else{
      bamGroupNumber = 0
      for (bamGroup <- bamGroups) {
        val mergeBam = new PrintReads with CommandLineGATKArgs
        mergeBam.input_file = Seq(mergeBamLists(bamGroupNumber))
        mergeBam.out = mergedBamDir + "%03d.bam".format(bamGroupNumber)
        mergeBam.memoryLimit = mergeBamMemoryLimit
        mergeBam.nct = 4 //todo is it make sense to use it with such an high memoryLimit?
        add(mergeBam)
        bamGroupNumber += 1
      }


    }




  }
}