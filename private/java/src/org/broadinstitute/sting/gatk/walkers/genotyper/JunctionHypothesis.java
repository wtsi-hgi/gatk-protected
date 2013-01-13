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

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.sting.utils.codecs.table.TableFeature;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 11/16/11
 * Time: 9:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class JunctionHypothesis implements Comparable, HasGenomeLocation {

    List<Pair<Integer,Integer>> exonJunctionNumbers;
    List<Pair<GenomeLoc,GenomeLoc>> exonJunctionLocs;
    private Map<GenomeLoc,Integer> baseBeforeCache;
    List<GenomeLoc> exons;
    boolean isForward;
    String sequence;
    String name;

    public JunctionHypothesis(String transcript, RefSeqFeature feature, List<Pair<Integer,Integer>> junctions, IndexedFastaSequenceFile refReader) {
        if ( ! transcript.equals(feature.getTranscriptUniqueGeneName()) ) {
            throw new ReviewedStingException("Attempting to form a junction sequence for a RefSeq record whose feature name does not match");
        }
        exons = new ArrayList<GenomeLoc>(feature.getExons()); // ensures same ordering as hypothesis
        Collections.sort(exons);
        exonJunctionNumbers = junctions;
        exonJunctionLocs = new ArrayList<Pair<GenomeLoc,GenomeLoc>>(exonJunctionNumbers.size());
        Integer last = null;
        StringBuilder seqBuilder = new StringBuilder();
        exons = new ArrayList<GenomeLoc>(exonJunctionLocs.size());
        for ( Pair<Integer,Integer> junction : exonJunctionNumbers ) {
            GenomeLoc first = feature.getSortedExonLoc(junction.first);
            GenomeLoc second = feature.getSortedExonLoc(junction.second);
            exonJunctionLocs.add(new Pair<GenomeLoc,GenomeLoc>(first,second));
            if ( last == null || last.equals(junction.first) ) {
                String exonSeq = new String(refReader.getSubsequenceAt(first.getContig(),first.getStart(),first.getStop()).getBases());
                seqBuilder.append(exonSeq);
            } else {
                throw new ReviewedStingException("Attempting to form a junction sequence consisting of discontinous ending and starting exons");
            }
            last = junction.second;
        }
        for ( Pair<GenomeLoc,GenomeLoc> locPair : exonJunctionLocs ) {
            if ( exons.size() == 0 ) {
                exons.add(locPair.first);
            }
            exons.add(locPair.second);
        }
        GenomeLoc lastLoc = exonJunctionLocs.get(exonJunctionLocs.size()-1).second;
        String lastString = new String(refReader.getSubsequenceAt(lastLoc.getContig(),lastLoc.getStart(),lastLoc.getStop()).getBases());
        seqBuilder.append(lastString);
        sequence = seqBuilder.toString();
        isForward = feature.getStrand() == 1;
        name = transcript;
        baseBeforeCache = new HashMap<GenomeLoc,Integer>(exonJunctionLocs.size());
    }

    public String getSequence() {
        return sequence;
    }

    public int size() {
        return sequence.length();
    }

    public int getBaseOffset(GenomeLoc position) {
        // count up the number of bases in exons prior to this position (involved in the hypothesis)
        if ( baseBeforeCache.size() == 0 ) {
            // fill in the cache
            boolean first = true;
            for ( Pair<GenomeLoc,GenomeLoc> jun : exonJunctionLocs ) {
                if ( first ) {
                    baseBeforeCache.put(jun.first,0);
                    first = false;
                }
                baseBeforeCache.put(jun.second,baseBeforeCache.get(jun.first)+ (int)jun.first.size());
            }
        }

        if ( getExonLoc(position) == null ) {
            System.out.println("foo");
        }
            GenomeLoc exon = getExonLoc(position);
        int nBefore = baseBeforeCache.get(exon);

        return nBefore + position.getStart()-exon.getStart();
    }

    public GenomeLoc getLocation() {
        GenomeLoc exonStart = exons.get(0).getStartLocation();
        GenomeLoc exonStop = exons.get(exons.size()-1).getStopLocation();
        return exonStart.endpointSpan(exonStop);
    }

    public int compareTo(Object other) {
        if ( other instanceof JunctionHypothesis ) {
            int locComp = getLocation().compareTo(((JunctionHypothesis) other).getLocation());
            if ( locComp == 0 ) {
                return this.name.compareTo(((JunctionHypothesis) other).name);
            }
        }
        return Integer.MIN_VALUE;
    }

    private double bestUnmappedScore = Double.POSITIVE_INFINITY;
    private double bestMappedScore = Double.POSITIVE_INFINITY;

    public void updateBestScore(double score, boolean unmapped) {
        if ( unmapped ) {
            bestUnmappedScore = score < bestUnmappedScore ? score : bestUnmappedScore;
        } else {
            bestMappedScore = score < bestMappedScore ? score : bestMappedScore;
        }
    }

    public double getScoreDifference() {
        return bestUnmappedScore - bestMappedScore;
    }

    public int getInsertAdjustment(GenomeLoc read, GenomeLoc mate) {
        int dist = -1;
        GenomeLoc ren = getExonLoc(read);
        if ( ren == null ) {
            return dist;
        }
        GenomeLoc men = getExonLoc(mate);
        if ( men == null || men.equals(ren) ) {
            return dist;
        }
        dist = men.distance(ren);
        return dist;
    }

    private GenomeLoc getExonLoc(GenomeLoc loc) {
        for ( GenomeLoc eLoc : exons ) {
            if ( eLoc.overlapsP(loc)) {
                return eLoc;
            } else if ( eLoc.isPast(loc) ) {
                return null;
            }
        }

        return null;
    }

    public boolean overlapsExonIntron(GenomeLoc loc) {
        for ( GenomeLoc eLoc : exons ) {
            if ( eLoc.overlapsP(loc) ) {
                // see if it starts before eloc or ends after eloc
                if ( loc.startsBefore(eLoc) ) {
                    // make sure that eloc isn't hte first
                    if ( ! eLoc.equals(exons.get(0)) ) {
                        return true;
                    } else {
                        return false;
                    }
                } else if ( eLoc.startsBefore(loc) ) {
                    // make sure that eloc isn't hte last
                    if ( ! eLoc.equals(exons.get(exons.size()-1))) {
                        return true;
                    } else {
                        return false;
                    }
                }
            }
        }

        return false;
    }

    public String toString() {
        return String.format("<%s:%s>",name,joinExonNums());
    }

    private String joinExonNums() {
        StringBuilder builder = new StringBuilder();
        if ( exonJunctionNumbers.size() <= 0 ) {
            return "";
        }
        for ( Pair<Integer,Integer> e : exonJunctionNumbers ) {
            builder.append(e.first);
            builder.append('-');
        }
        builder.append(exonJunctionNumbers.get(exonJunctionNumbers.size()-1).second);
        return builder.toString();
    }

    public static List<JunctionHypothesis> generateHypotheses(List<TableFeature> rawHypos, Map<String,RefSeqFeature> refSeqFeaturesByGene, IndexedFastaSequenceFile reader) {
        List<JunctionHypothesis> hypotheses = new ArrayList<JunctionHypothesis>(rawHypos.size());
        for ( TableFeature rhy : rawHypos ) {
            RefSeqFeature feature = refSeqFeaturesByGene.get(rhy.getValue(1));
            Set<List<Pair<Integer,Integer>>> pairs = unwrapPairs(rhy.getValue(3));
            for ( List<Pair<Integer,Integer>> hypo : pairs ) {
                hypotheses.add(new JunctionHypothesis(rhy.getValue(2),feature,hypo,reader));
            }
        }

        return hypotheses;
    }

    public static Set<List<Pair<Integer,Integer>>> unwrapPairs(String rawString) {
        Set<List<Pair<Integer,Integer>>> pairSet = new HashSet<List<Pair<Integer,Integer>>>();
        for ( String hypo : rawString.split("\\|") ) {
            List<Pair<Integer,Integer>> pairs = new ArrayList<Pair<Integer,Integer>>(hypo.split(";").length);
            for ( String rawPair : hypo.split(";")) {
                String[] exons = rawPair.split(",");
                Integer e1 = Integer.parseInt(exons[0]);
                Integer e2 = Integer.parseInt(exons[1]);
                Pair<Integer,Integer> pair = new Pair<Integer,Integer>(e1,e2);
                pairs.add(pair);
            }

            pairSet.add(pairs);
        }

        return pairSet;
    }
}
