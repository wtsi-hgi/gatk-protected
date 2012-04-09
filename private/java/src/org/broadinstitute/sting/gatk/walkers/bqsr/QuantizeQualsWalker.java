/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.bqsr;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.recalibration.QualQuantizer;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Quantizes the quality scores in a BAM file
 *
 * <p>
 * x
 *
 * <h2>Input</h2>
 * <p>
 * One or more bam files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A single processed bam file.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T QuantizeQualsWalker \
 *   -o output.bam \
 *   -I input1.bam -Q qual_dist.gatkreport.txt
 * </pre>
 *
 */
public class QuantizeQualsWalker extends ReadWalker<SAMRecord, SAMFileWriter> {
    final protected static Logger logger = Logger.getLogger(QuantizeQualsWalker.class);

    @Output(doc="Write output to this BAM filename instead of STDOUT")
    SAMFileWriter out;

    @Argument(fullName = "qualityHistogram", shortName = "Q", doc="", required = true)
    File qualHistogramFile;

    @Argument(fullName = "quantizationLevels", shortName = "quantizationLevels", doc="The number of quality levels to include in output", required = false)
    int nQualityLevels = 8;

    @Argument(fullName = "minInterestingQual", shortName = "minInterestingQual", doc="Quality scores less than or equal to this value are considered uninteresting, are can be freely merged together", required = false)
    int minInterestingQual = 10;

    @Output(fullName = "report", shortName = "report", doc="Write GATK report of quantization process to this file", required = false)
    PrintStream reportOut = null;

    @Argument(fullName = "maxQualToInclude", shortName = "maxQualToInclude", doc="Only quality scores <= this value are considered for remapping", required = false)
    private int MAX_QUAL_TO_INCLUDE = QualityUtils.MAX_QUAL_SCORE;

    private QualQuantizer quantizer;

    /**
     * The initialize function.
     */
    public void initialize() {
        GATKReport qualHist = new GATKReport(qualHistogramFile);
        GATKReportTable table = qualHist.getTable("QualityScoreDistribution");

        List<Long> nObservationsPerQual = new ArrayList<Long>(MAX_QUAL_TO_INCLUDE+1);
        for ( int q = 0; q <= MAX_QUAL_TO_INCLUDE; q++) {
            long count = Long.valueOf(table.get(q, "Count").toString());
            logger.info(String.format("%3d %10d", q, count));
            nObservationsPerQual.add(count);
        }

        quantizer = new QualQuantizer(nObservationsPerQual, nQualityLevels, minInterestingQual);
        if ( reportOut != null ) {
            quantizer.writeReport(reportOut);
        }
    }

    /**
     * The reads map function.
     *
     * @param ref  the reference bases that correspond to our read, if a reference was provided
     * @param read the read itself, as a SAMRecord
     * @return the read itself
     */
    public SAMRecord map( ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker ) {
        final byte[] quantizedQuals = read.getBaseQualities().clone();

        for ( int i = 0; i < quantizedQuals.length; i++ ) {
            int oq = quantizedQuals[i];
            byte nq = quantizer.getOriginalToQuantizedMap().get(oq);
            quantizedQuals[i] = nq;
        }

//        if ( logger.isDebugEnabled() ) {
//            logger.info("Readname: " + read.getReadName());
//            logger.info("      OQ: " + SAMUtils.phredToFastq(read.getBaseQualities()));
//            logger.info("      NQ: " + SAMUtils.phredToFastq(quantizedQuals));
//        }

        read.setBaseQualities(quantizedQuals); // Overwrite old qualities with new recalibrated qualities
        return read;
    }

    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     *
     * @return SAMFileWriter, set to the BAM output file if the command line option was set, null otherwise
     */
    public SAMFileWriter reduceInit() {
        return out;
    }

    /**
     * given a read and a output location, reduce by emitting the read
     *
     * @param read   the read itself
     * @param output the output source
     * @return the SAMFileWriter, so that the next reduce can emit to the same source
     */
    public SAMFileWriter reduce( SAMRecord read, SAMFileWriter output ) {
        output.addAlignment(read);
        return output;
    }

}
