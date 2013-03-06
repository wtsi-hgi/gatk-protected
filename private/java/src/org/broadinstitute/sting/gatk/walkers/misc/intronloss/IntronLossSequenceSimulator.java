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

package org.broadinstitute.sting.gatk.walkers.misc.intronloss;

import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
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
        FastqWriterFactory fastqWriterFactory = new FastqWriterFactory();
        readWriter = fastqWriterFactory.newWriter(m1OutputFastQ);
        mateWriter = fastqWriterFactory.newWriter(m2OutputFastQ);
        widowedWriter = fastqWriterFactory.newWriter(widowedOutputFASTQ);
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
                    GATKReportTable reportTable = report.getTable("InsertSizeDistributionByReadGroup");
                    // rows are insert sizes, columns are read groups
                    for (GATKReportColumn reportColumn : reportTable.getColumnInfo() ) {
                        // annoyingly, the column has no knowledge of its own rows
                        int sum = 0;
                        for ( int row = 0; row < reportTable.getNumRows(); row++ ) {
                            sum += Integer.parseInt( (String) reportTable.get(row,reportColumn.getColumnName()));
                        }
                        byte[] rgHist = new byte[MAX_INSERT_SIZE];
                        for ( int row = 0; row < reportTable.getNumRows(); row++ ) {
                            final int insertSize = Integer.parseInt( (String) reportTable.get(row,0));
                            int val = 1;
                            if ( insertSize < rgHist.length ) {
                                val = Integer.parseInt( (String) reportTable.get(row,reportColumn.getColumnName()));
                            }
                            rgHist[row] = QualityUtils.errorProbToQual((((double) val) / sum), QualityUtils.MAX_SAM_QUAL_SCORE);
                        }

                        insertSizeDistribution = qualToDouble(rgHist);
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

    private double[] qualToDouble(byte[] hist) {
        double[] d = new double[hist.length];
        for ( int i = 0; i < d.length; i++) {
            d[i] = QualityUtils.qualToErrorProb(hist[i]);
        }

        return d;
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
                    FastqRecord[] records = generateFastqFromSequence(normal, offset,outputSequences,false);
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
                FastqRecord[] records = generateFastqFromSequence(loss,offset,outputSequences,true);
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

    private FastqRecord[] generateFastqFromSequence(String seq, int offset, int seqOutput,boolean isLoss) {
        FastqRecord[] toRet;
        int insert = generateInsertSize();
        boolean widowed = (offset + insert - READ_SIZE_BP ) < 0 || (offset + insert) > seq.length();
        //logger.debug(String.format("seqSize: %d, offset: %d, insert: %d, widowed: %s",seq.length(),offset,insert,widowed));
        String readName = String.format("%s:%d%s %s",seqPrefix,seqOutput, widowed ? "" : " 1:N:0:ATGCGC", isLoss ? ":loss" : ":noLoss");
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
            return (byte) ( Math.min(40, Math.max(2, (int) 25.5 + GenomeAnalysisEngine.getRandomGenerator().nextGaussian() * 5)));
        }
    }

    private int drawFromCumulativeHistogram(double[] histogram) {
        double q = GenomeAnalysisEngine.getRandomGenerator().nextDouble();
        int offset = 0;
        double cumsum = 0.0;
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