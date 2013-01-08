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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import com.google.java.contract.Requires;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.variant.utils.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.variant.utils.BaseUtils;

import java.util.*;

/**
 *
 */
public class IntronLossGenotypeLikelihoodCalculationModel {

    public static final int PLOIDY = 2;
    public static final int ALIGNMENT_QUAL_CAP = 500;
    private static final double NEGLOG2 = -Math.log10(2);
    private static final double NEGLOG3 = -Math.log10(3);
    private static final double GARBAGE_COLLECT_THRESHOLD = -20; // Q200
    private GenomeLocParser genomeLocParser;
    private Map<String,byte[]> insertSizeQualHistograms;
    private Set<String> samples;
    private Logger logger;
    private boolean fullSmithWaterman;

    public IntronLossGenotypeLikelihoodCalculationModel(Logger logger, boolean fullSW) {
        this.logger = logger;
        fullSmithWaterman = fullSW;
    }

    public void initialize(RodBinding<RefSeqFeature> refSeqBinding,
                           GenomeLocParser parser, Map<String,byte[]> histograms) {
        this.genomeLocParser = parser;
        this.insertSizeQualHistograms = histograms;
    }

    public void setSamples(Set<String> s) {
        samples =s;
    }

    private void updateLikelihoods(double[] likelihoods, double pH0Lik, double pH1Lik) {
        // (AC/AN)*P(read|H1) + ( (AN-AC)/AN)*P(read|H0)
        // todo -- there is an unaccounted-for alignment bias
        // 0 -- AC0 -- zero retrocopies
        likelihoods[0] += pH0Lik;
        // 1 -- AC1 -- 1 of 3 reads should be retrocopy
        likelihoods[1] += MathUtils.log10sumLog10(new double[]{ NEGLOG3 + pH1Lik,NEGLOG3 - NEGLOG2 +pH0Lik});
        // 2 -- AC2 -- 2 of 4 reads should be retrocopy
        likelihoods[2] += MathUtils.log10sumLog10(new double[]{NEGLOG2 + pH1Lik,NEGLOG2 +pH0Lik});
    }

    public Pair<String,double[]> getLikelihoods(JunctionHypothesis hypothesis, GATKSAMRecord read) {
        // get the probabilities of the pair, if this is the first read (don't double count)
        double[] likelihoods = new double[3];
        double pInsertH0 = 0.0;
        double pInsertH1 = 0.0;
        if ( read.getReadPairedFlag() && ! read.getMateUnmappedFlag() && read.getFirstOfPairFlag() ) {
            // check to see if read and mate overlap one or more exons in the hypothesis
            GenomeLoc readLoc = getGenomeLoc(read);
            GenomeLoc mateLoc = genomeLocParser.createGenomeLoc(read.getMateReferenceName(),read.getMateAlignmentStart(),read.getMateAlignmentStart()+read.getReadLength());
            int insertAdjustment = hypothesis.getInsertAdjustment(readLoc,mateLoc);
            if ( insertAdjustment > 0 ) {
                pInsertH0 += insertSizeProbability(read,0);
                pInsertH1 += insertSizeProbability(read,insertAdjustment);
            }
        }

        GenomeLoc readLoc = read.getReadUnmappedFlag() ? null : genomeLocParser.createGenomeLoc(read);

        if ( read.getReadUnmappedFlag() || hypothesis.overlapsExonIntron(readLoc) ) {
            pInsertH0 += scoreAlignment(read);
            pInsertH1 += scoreRealignment(read,readLoc,hypothesis);
        }

        if ( pInsertH0 > GARBAGE_COLLECT_THRESHOLD || pInsertH1 > GARBAGE_COLLECT_THRESHOLD ) {
            updateLikelihoods(likelihoods,pInsertH0,pInsertH1);
        }

        logger.debug(String.format("%s :: (%f,%f) -- [0] : %f  [1] : %f  [2] : %f",read.getReadName() + (read.getReadUnmappedFlag() ? " (u)" : "(m)"),pInsertH0,pInsertH1,likelihoods[0],likelihoods[1],likelihoods[2]));

        return new Pair<String, double[]>(read.getReadGroup().getSample(),likelihoods);
    }

    private double scoreAlignment(GATKSAMRecord read) {
        // aligned reads are equipped already with necessary metrics
        double logScore = 0.0;
        if ( read == null )
            throw new StingException("Read cannot be null here");
        if ( read.getReadUnmappedFlag() ) {
            return -50.0;
        }
        logScore += ( (Integer) read.getAttribute("NM") )*(-1.2);
        for ( CigarElement e : read.getCigar().getCigarElements() ) {
            if ( e.getOperator().equals(CigarOperator.D) || e.getOperator().equals(CigarOperator.I) ) {
                logScore += -1*(2.0-Math.pow(2,1-e.getLength()));
            } else if ( e.getOperator().equals(CigarOperator.S) ) {
                logScore += -0.50*e.getLength();
            }
        }

        return logScore;
    }

    private double scoreRealignment(GATKSAMRecord read, GenomeLoc readLoc, JunctionHypothesis hypothesis) {
        if ( fullSmithWaterman || read.getReadUnmappedFlag() ) {
            SWPairwiseAlignment alignment = new SWPairwiseAlignment(hypothesis.getSequence().getBytes(),read.getReadBases());
            double score = scoreRealignment(alignment.getCigar().getCigarElements(),alignment.getAlignmentStart2wrt1(),read.getReadBases(),hypothesis.getSequence().getBytes());
            hypothesis.updateBestScore(score,true);
            return score;
        } else {
            // want to use the BWA alignment of the read at its present position to score the alignment
            // this makes the somewhat unfortunate assumption that the cigar operator will just be M for all the
            // overhanging bases. Combined with the GARBAGE_COLLECT_THRESHOLD above, this could easily drop
            // reads from consideration entirely, as every base will mismatch following any indel. Perhaps
            // it would be worth adding in a check for many reads being dropped and recommending that the user
            // run the algorithm with fullSmithWaterman turned on.

            int offset = hypothesis.getBaseOffset(readLoc);
            if ( offset < 0 || offset > hypothesis.size() ) {
                // this occurs if the first exon is <60bp, the read aligned into the intron with
                // small overlap to second exon, so the offset is off the hypothesis
                return -30.0;
            }
            // trailing or starting soft-clips want to become Ms
            Collection<CigarElement> oldElements = read.getCigar().getCigarElements();
            List<CigarElement> newElements = new ArrayList<CigarElement>(oldElements.size());
            for ( CigarElement e : oldElements ) {
                if ( e.getOperator().equals(CigarOperator.S) ) {
                    newElements.add(new CigarElement(e.getLength(),CigarOperator.M));
                } else {
                    newElements.add(e);
                }
            }
            double score =  scoreRealignment(newElements,offset,read.getReadBases(),hypothesis.getSequence().getBytes());
            hypothesis.updateBestScore(score,true);
            return score;
        }
    }

    public GenomeLoc getGenomeLoc(GATKSAMRecord read) {
        if ( read == null )
            return null;
        int start =  read.getSoftStart();
        int end =  Math.max(start, read.getSoftEnd());
        return genomeLocParser.createGenomeLoc(read.getReferenceName(),start,end);
    }

    private double scoreRealignment(Collection<CigarElement> elements, int offset, byte[] readBases, byte[] exonBases) {
        int b = 0;
        int phredBasesAreExonBases = 0;
        int indelPen = 15;
        for ( CigarElement elem : elements) {
            if ( elem.getOperator().equals(CigarOperator.M) ) {
                for ( int ct = elem.getLength(); ct > 0; ct-- ) {
                    int pos = b + offset;
                    if (pos >= exonBases.length || BaseUtils.basesAreEqual(readBases[b], exonBases[pos])) {
                        // bases match, don't do anything
                    } else {
                        // read clipped at Q15, with e/3 this goes to Q12
                        phredBasesAreExonBases += 12;
                    }
                    b++;
                }
            } else if ( elem.getOperator().equals(CigarOperator.D) ) {
                // just increment offset
                offset += elem.getLength();
                phredBasesAreExonBases += indelPen + 15*elem.getLength();
                indelPen += 15;
            } else if (elem.getOperator().equals(CigarOperator.I)) {
                // update the b
                b += elem.getLength();
                offset -= elem.getLength(); // compensate, since pos will move up by elem.getLength()
                phredBasesAreExonBases += indelPen+15*elem.getLength();
                indelPen += 15;
            } else if ( elem.getOperator().equals(CigarOperator.S) ) {
                // just skip the bases and update b; soft clipping occurs at edges, so need to subtract from offset
                b += elem.getLength();
                offset -= elem.getLength();
                phredBasesAreExonBases += 10*elem.getLength();
            } else if ( elem.getOperator().equals(CigarOperator.H) ) {
               // offset -= elem.getLength();
            } else {
                throw new StingException("Unsupported operator: "+elem.getOperator().toString());
            }
        }

        int cappedPhred = phredBasesAreExonBases > ALIGNMENT_QUAL_CAP  ? ALIGNMENT_QUAL_CAP : phredBasesAreExonBases;

        return  ( (double) cappedPhred)/(-10.0);
    }
    private double insertSizeProbability(GATKSAMRecord read, int adjustment) {
        byte[] quals = insertSizeQualHistograms.get(read.getReadGroup().getId());
        if ( quals == null ) {
            logger.warn("No insert size histogram for RGID: " + read.getReadGroup().getId()+" returning an uninformative probability");
            return QualityUtils.qualToErrorProb((byte) 4);
        }
        int size = Math.max(1,Math.min(quals.length-1,Math.abs(read.getInferredInsertSize())-adjustment));
        return Math.log10(QualityUtils.qualToErrorProb(quals[size]));
    }

    enum PostulateClass {
        LOSS, NONE
    }

    class IntronLossPostulate {

        Pair<GenomeLoc,GenomeLoc> exonLocs;
        Pair<Integer,Integer> exonNums;
        PostulateClass clazz;
        Set<Pair<GATKSAMRecord,GATKSAMRecord>> supportingPairs = new HashSet<Pair<GATKSAMRecord,GATKSAMRecord>>(1000);

        public IntronLossPostulate(List<GenomeLoc> exons, GenomeLoc read, GenomeLoc mate) {
            if ( read == null || mate == null ) {
                clazz = PostulateClass.NONE;
                exonLocs = null;
                exonNums = null;
                return;
            }
            int eReadNum = -1, eMateNum = -1;
            int eNum = 0;
            for ( GenomeLoc e : exons ) {
                if ( e.overlapsP(read) ) {
                    eReadNum = eNum;
                }

                if ( e.overlapsP(mate) ) {
                    eMateNum = eNum;
                }
                eNum++;
            }

            if ( eReadNum != eMateNum && eMateNum > -1 && eReadNum > -1 ) {
                clazz = PostulateClass.LOSS;
                exonLocs = new Pair<GenomeLoc, GenomeLoc>(exons.get(eReadNum),exons.get(eMateNum));
                exonNums = new Pair<Integer, Integer>(eReadNum,eMateNum);
            } else {
                clazz = PostulateClass.NONE;
                exonLocs = null;
                exonNums = null;
            }
        }

        public String testConstructor(List<GenomeLoc> exons, GenomeLoc read, GenomeLoc mate) {
            StringBuffer buf = new StringBuffer();
            int eReadNum = -1, eMateNum = -1;
            int eNum = 0;
            buf.append(read);
            buf.append(" :: ");
            buf.append(mate);
            for ( GenomeLoc e : exons ) {
                buf.append("\t");
                buf.append(e);
                if ( e.overlapsP(read) ) {
                    eReadNum = eNum;
                    buf.append("\tREAD_OVERLAP");
                }

                if ( e.overlapsP(mate) ) {
                    eMateNum = eNum;
                    buf.append("\tMATE_OVERLAP");
                }
            }

            return buf.toString();
        }

        public boolean equals(Object other) {
            if ( ! ( other instanceof IntronLossPostulate ) ) {
                return false;
            }
            IntronLossPostulate otherPostulate = (IntronLossPostulate) other;
            return ( (exonNums.first == otherPostulate.exonNums.first) &&
                    exonNums.second == otherPostulate.exonNums.second ) ||
                    ( ( exonNums.second == otherPostulate.exonNums.first) &&
                    ( exonNums.first == otherPostulate.exonNums.second) );
        }

        public int hashCode() {
            return exonNums.first.hashCode() ^ exonNums.second.hashCode();
        }

        public String toString() {
            if ( clazz == PostulateClass.NONE ) {
                return "NONE";
            }

            return String.format("EX%d-EX%d",exonNums.first,exonNums.second);
        }

        public int getEventSize() {
            return exonLocs.first.minDistance(exonLocs.second);
        }

        public boolean isValid() {
            // todo -- this needs to somehow use the insert size distribution
            return getEventSize() > 350;
        }

        public byte[] getSequence(IndexedFastaSequenceFile ref) {
            // BIG TODO -- SPLICE DONOR/ACCEPTOR SITES AND STRAND OF TRANSCRIPT
            StringBuffer buf = new StringBuffer();
            buf.append(new String(ref.getSubsequenceAt(exonLocs.first.getContig(),exonLocs.first.getStart(),exonLocs.first.getStop()).getBases()));
            buf.append(new String(ref.getSubsequenceAt(exonLocs.second.getContig(),exonLocs.second.getStart(),exonLocs.second.getStop()).getBases()));
            return buf.toString().getBytes();
        }

        public boolean isRelevant(GenomeLoc loc1, GenomeLoc loc2) {
            if ( loc1 != null && isRelevant(loc1) ) {
                // loc 2 can't be way too far away
                return loc2 == null || loc2.overlapsP(exonLocs.first) || loc2.overlapsP(exonLocs.second);
            } else if ( loc2 != null && isRelevant(loc2) ) {
                return loc1 == null || loc1.overlapsP(exonLocs.first) || loc1.overlapsP(exonLocs.second);
            }

            return false;
        }

        @Requires("! loc.overlapsP(exonLocs.second)")
        public boolean isRelevant(GenomeLoc loc) {
            // basically if the read hangs partially in the intron
            boolean r1 =  loc.overlapsP(exonLocs.first) && loc.getStop() > exonLocs.first.getStop();
            boolean r2 = loc.overlapsP(exonLocs.second) && loc.getStart() < exonLocs.second.getStart();
            return r1 || r2;
        }
    }
}
