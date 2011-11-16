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

package org.broadinstitute.sting.gatk.walkers.misc.intronloss;

import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import org.apache.commons.math.stat.descriptive.rank.Max;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportColumn;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Simulates reads from an intron deletion (or not)
 */
public class IntronLossSequenceSimulator extends RefWalker<Pair<Byte,Boolean>,Pair<StringBuffer,StringBuffer>> {

    @Input(shortName = "rs", fullName = "refSeq", doc = "RefGene file", required = true)
    public RodBinding<RefSeqFeature> refSeqRodBinding;

    @Output(shortName = "m1",fullName = "m1OutputFASTQ", doc = "the output FASTQ file for mate 1", required = true)
    public File m1OutputFastQ;

    @Output(shortName = "m2", fullName = "m2OutputFASTQ", doc = "the output FASTQ file for mate 2", required = true)
    public File m2OutputFastQ;

    @Output(shortName = "w", fullName = "wOutputFASTQ", doc = "the output FASTQ for widowed reads", required = true)
    public File widowedOutputFASTQ;

    @Argument(fullName = "ploidy", required = false, doc = "")
    public int ploidy = 2;

    @Argument(fullName = "nVar", required = false, doc = "")
    public int nVar = 1;

    @Argument(fullName = "sequencePrefix", required = false, doc = "")
    public String seqPrefix = "ILSS-SRR";

    @Argument(fullName = "sequenceCoverage", shortName = "cov", required = false, doc = "Coverage that sequence should attain (average)")
	public int AVG_CVG_PER_CHR_POISSON = 15;

    @Input(shortName="H",fullName="insertHistogram",required=false,doc="The insert size histogram per read group, either flat file or GATK report formatted")
    public File readGroupInsertHistogram = null;


    @Input(fullName="raw_qual_histogram",shortName="Q",required=false,doc="File holding the raw quality score histogram")
    public File RAW_QSCORE_HIST = null;

    private static final int READ_SIZE_BP = 76;

    private FastqWriter readWriter;
    private FastqWriter mateWriter;
    private FastqWriter widowedWriter;

    private double[] insertSizeDistribution = null;
    private double[][] qualityScoreHistogram = null;
    private final int MAX_INSERT_SIZE = 1000;

    public void initialize() {
        readWriter = new FastqWriter(m1OutputFastQ);
        mateWriter = new FastqWriter(m2OutputFastQ);
        widowedWriter = new FastqWriter(widowedOutputFASTQ);
        if (readGroupInsertHistogram != null ) {
            try {

                XReadLines xrl = new XReadLines(readGroupInsertHistogram);
                if ( ! xrl.next().startsWith("##:") ) {
                    xrl.close();
                    for ( String entry : new XReadLines(readGroupInsertHistogram) ) {
                        String[] split1 = entry.split("\\t");
                        String id = split1[0];
                        String[] histogram = split1[1].split(";");
                        double[] quals = new double[histogram.length];
                        int idx = 0;
                        for ( String histEntry : histogram ) {
                            quals[idx++] = Double.parseDouble(histEntry);
                        }

                        insertSizeDistribution = quals;
                    }
                } else {
                    xrl.close();
                    GATKReport report = new GATKReport(readGroupInsertHistogram);
                    GATKReportTable reportTable = report.getTable("InsertSizeDistribution");
                    // rows are insert sizes, columns are read groups
                    for (GATKReportColumn reportColumn : reportTable.getColumns() ) {
                        // annoyingly, the column has no knowledge of its own rows
                        int sum = 0;
                        for ( int row = 0; row < reportTable.getNumRows(); row++ ) {
                            sum += Integer.parseInt( (String) reportTable.get(row,reportColumn.getColumnName()));
                        }
                        int remain = sum;
                        double[] rgHist = new double[MAX_INSERT_SIZE];
                        for ( int row = 0; row < rgHist.length; row++) {
                            rgHist[row] = ((double)remain)/sum;
                            if ( reportTable.containsKey(row) ) {
                                remain -= Integer.parseInt( (String) reportTable.get(row,reportColumn.getColumnName()));
                            }
                        }

                        insertSizeDistribution = rgHist;
                    }
                }
            } catch ( Exception e ) {
                throw new ReviewedStingException("Exception: ",e);
            }

            logger.debug(Arrays.deepToString(org.apache.commons.lang.ArrayUtils.toObject(insertSizeDistribution)));
        }

        if ( RAW_QSCORE_HIST != null ) {
            try {
                qualityScoreHistogram = new double[READ_SIZE_BP][41];
                XReadLines xrl = new XReadLines(RAW_QSCORE_HIST);
                for ( int lineNo = 0; lineNo < READ_SIZE_BP; lineNo ++ ) {
                    String[] quals = xrl.next().split(",");
                    int offset = 0;
                    for ( String q : quals ) {
                        if ( offset > 40 ) { break; }
                        qualityScoreHistogram[lineNo][offset++] = Double.parseDouble(q);
                    }
                }
            }  catch (java.io.IOException e ) {
                throw new UserException("File could not be opened",e);
            }
        }
    }

    public Pair<StringBuffer,StringBuffer> reduceInit() {
        return new Pair<StringBuffer, StringBuffer>(new StringBuffer(), new StringBuffer());
    }

    public boolean isReduceByInterval() { return true; }

    public Pair<StringBuffer,StringBuffer> reduce(Pair<Byte,Boolean> map, Pair<StringBuffer,StringBuffer> reduce) {
        if ( map == null || reduce == null ) {
            return reduce;
        }

        if ( map.second ) {
            // if in an exon
            reduce.first.append((char) map.first.byteValue());
        }

        reduce.second.append((char) map.first.byteValue());

        return reduce;
    }

    public Pair<Byte,Boolean> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || ! tracker.hasValues(refSeqRodBinding) ) {
            return new Pair<Byte, Boolean>(ref.getBase(),true);
        }

        //RefSeqFeature geneFeature = getProperFeature(tracker.getValues(refSeqRodBinding));
        RefSeqFeature geneFeature = tracker.getFirstValue(refSeqRodBinding);
        if ( geneFeature == null ) {
            return new Pair<Byte, Boolean>(ref.getBase(),true);
        }

        // todo: probably cache this
        List<GenomeLoc> exons = geneFeature.getExons();
        boolean inExon = false;
        for ( GenomeLoc exon : exons ) {
            inExon |= ref.getLocus().overlapsP(exon);
        }

        return new Pair<Byte, Boolean>(ref.getBase(),inExon);
    }

    private RefSeqFeature getProperFeature(List<RefSeqFeature> refSeqRawFeatures) {
        for ( Feature f : refSeqRawFeatures ) {
            if ( f.getClass().isAssignableFrom(RefSeqFeature.class) ) {
                return (RefSeqFeature) f;
            }
        }

        return null;
    }

    public void onTraversalDone(Pair<StringBuffer,StringBuffer> pair) { }

    public void onTraversalDone(List<Pair<GenomeLoc,Pair<StringBuffer,StringBuffer>>> sequenceByGene) {
        for ( Pair<GenomeLoc,Pair<StringBuffer,StringBuffer>> genePair : sequenceByGene ) {
            String intronLossStr = genePair.second.first.toString();
            String normalStr = genePair.second.second.toString();
            // todo -- generalize random process (uniform/poisson/etc)
            simulateSequencing(intronLossStr, normalStr, nVar, ploidy);
        }
    }

    public void simulateSequencing(String loss, String normal, int nVar, int ploidy) {
        // what is the size of the read pool (normal)
        int nReadPool = normal.length()-READ_SIZE_BP+1;
        // now what is the probability p with which we uniformly select each read to obtain
        // the desired expected coverage, we have that coverage at a specific site is READ_SIZE_BP copies
        // of the random variable, so it's CVG*SEQ_SIZE = p*2*READ_SIZE*READ_POOL -- the factor of 2 is due to pairs
        double pSeq = ( (double) AVG_CVG_PER_CHR_POISSON*normal.length())/( (double) READ_SIZE_BP*nReadPool*2 );

        if ( pSeq > 1.0 ) {
            throw new UserException("Cannot simulate to given coverage without introducing duplication");
        }

        //logger.debug(String.format("Normal length: %d Var Length: %d",normal.length(),loss.length()));

        int outputSequences = 0;

        for ( int offset = 0; offset < normal.length()-READ_SIZE_BP; offset++ ) {
            for ( int normChrom = 0; normChrom < ploidy-nVar; normChrom++ ) {
                if ( GenomeAnalysisEngine.getRandomGenerator().nextDouble() < pSeq ) {
                    FastqRecord[] records = generateFastqFromSequence(normal, offset,outputSequences);
                    if ( records.length > 1 ) {
                        readWriter.write(records[0]);
                        mateWriter.write(records[1]);
                    } else {
                        widowedWriter.write(records[0]);
                    }
                    outputSequences ++;
                }
            }
        }

        nReadPool = loss.length()-READ_SIZE_BP+1;
        pSeq = ( (double) AVG_CVG_PER_CHR_POISSON*loss.length())/( (double) READ_SIZE_BP*nReadPool*2 );
        for ( int offset = 0; offset < loss.length()-READ_SIZE_BP; offset++) {
            for ( int varChrom = 0; varChrom < nVar; varChrom++ ) {
                if ( GenomeAnalysisEngine.getRandomGenerator().nextDouble() > pSeq ) { continue; }
                FastqRecord[] records = generateFastqFromSequence(loss,offset,outputSequences);
                if ( records.length > 1 ) {
                    readWriter.write(records[0]);
                    mateWriter.write(records[1]);
                } else {
                    widowedWriter.write(records[0]);
                }
                outputSequences++;
            }
        }
    }

    private FastqRecord[] generateFastqFromSequence(String seq, int offset, int seqOutput) {
        FastqRecord[] toRet;
        int insert = generateInsertSize();
        boolean widowed = (offset + insert - READ_SIZE_BP ) < 0 || (offset + insert) > seq.length();
        //logger.debug(String.format("seqSize: %d, offset: %d, insert: %d, widowed: %s",seq.length(),offset,insert,widowed));
        String readName = String.format("%s:%d%s",seqPrefix,seqOutput, widowed ? "" : " 1:N:0:ATGCGC");
        Pair<String,String> readSequence = generateReadSequence(seq,offset,insert<0);
        FastqRecord record = new FastqRecord(readName,readSequence.first,readName,readSequence.second);
        if ( widowed ) {
            return new FastqRecord[]{ record };
        }

        String mateName =  String.format("%s:%d%s",seqPrefix,seqOutput, " 2:N:0:ATGCGC");
        Pair<String,String> mateSequence = generateReadSequence(seq,offset+insert-READ_SIZE_BP,!(insert<0));
        FastqRecord mateRecord = new FastqRecord(mateName,mateSequence.first,mateName,mateSequence.second);

        return new FastqRecord[]{ record, mateRecord };
    }

    private int generateInsertSize() {
        if (insertSizeDistribution == null ) {
            int symmetricUnderlying = (int)( 260.5 +  GenomeAnalysisEngine.getRandomGenerator().nextGaussian()*30);
            int chiSquareAdd = (int) ( 0.5 + Math.pow(GenomeAnalysisEngine.getRandomGenerator().nextGaussian()*8,2) );
            return (GenomeAnalysisEngine.getRandomGenerator().nextBoolean() ? 1 : -1 ) * (symmetricUnderlying + chiSquareAdd);
        } else {
            return drawFromCumulativeHistogram(insertSizeDistribution);
        }
    }

    private Pair<String,String> generateReadSequence(String base, int offset, boolean reverseComplement) {
        // todo -- make the random process more real (investigate using site qualities and TT)
        // todo -- incorporate indels as a possibility
        String rawSeq = base.substring(offset,offset+READ_SIZE_BP);
        StringBuffer newSeq = new StringBuffer();
        StringBuffer qualSeq = new StringBuffer();
        int charOffset = 0;
        for ( char b : rawSeq.toCharArray() ) {
            byte qual = generateQualityScore(charOffset);
            char observedB;
            if ( GenomeAnalysisEngine.getRandomGenerator().nextDouble() < QualityUtils.qualToErrorProb(qual) ) {
                observedB = (char) BaseUtils.baseIndexToSimpleBase(BaseUtils.getRandomBaseIndex(BaseUtils.simpleBaseToBaseIndex((byte) b)));
            } else {
                observedB = b;
            }
            newSeq.append(observedB);
            qualSeq.append((char) (qual + 33));
            charOffset++;
        }

        String newSeqStr = newSeq.toString();
        String qualSeqStr = qualSeq.toString();
        if ( reverseComplement ) {
            newSeqStr = BaseUtils.simpleReverseComplement(newSeq.toString());
            qualSeqStr = new StringBuffer(qualSeqStr).reverse().toString();
        }

        return new Pair<String, String>(newSeqStr,qualSeqStr);
    }

    private byte generateQualityScore(int cycle) {
        if ( qualityScoreHistogram != null ) {
            double[] qualHist = qualityScoreHistogram[cycle];
            return (byte) drawFromCumulativeHistogram(qualHist);
        } else {
            // qual: normal with mean 25, SD 10, quantized by rounding to nearest number
            return (byte) ( Math.min(40, Math.max(2, (int) 25.5 + GenomeAnalysisEngine.getRandomGenerator().nextGaussian() * 10)));
        }
    }

    private int drawFromCumulativeHistogram(double[] histogram) {
        double q = GenomeAnalysisEngine.getRandomGenerator().nextDouble();
        int offset = 0;
        double cumsum = 0.0;
        //logger.debug(Arrays.deepToString(org.apache.commons.lang.ArrayUtils.toObject(histogram)));
        do {
           // logger.debug(String.format("%f\t%f",cumsum,q));
            cumsum+= histogram[offset++];
        } while ( cumsum < q && offset < histogram.length);
        return offset;
    }

    private int descendingBinarySearch(byte[] array, byte q) {
        // pick mid
        int left = 0;
        int rt = array.length-1;
        int mid = left + (rt - left)/2;
        while ( array[mid] != q && left < rt ) {
            if ( array[mid] > q ) {
                left = mid+1;
            } else {
                rt = mid-1;
            }
            mid = left + (rt - left)/2;
        }
        return mid;
    }
}