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
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.variant.utils.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 11/16/11
 */
@ReadFilters({DuplicateReadFilter.class,FailsVendorQualityCheckFilter.class,MappingQualityZeroFilter.class})
public class ExonJunctionHypothesisGenerator extends ReadWalker<TreeSet<ExonJunctionHypothesisGenerator.IntronLossJunctions>,TreeSet<ExonJunctionHypothesisGenerator.IntronLossJunctions>> {

    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
    @Output
    protected PrintStream out = null;

    @Input(shortName="r",fullName="refSeq",required=true,doc="The RefSeq Gene definition track")
    public RodBinding<RefSeqFeature> refSeqRodBinding;

    @Input(shortName="mis",fullName="minInsertSize",required=false,doc="The minimum insert size for a read pair to consider as evidence of joined exons")
    public int minInferredInsert = 600;


    private IndexedFastaSequenceFile referenceReader;

    private TreeMap<GenomeLoc,byte[]> exonSequences;

    public TreeSet<IntronLossJunctions> reduceInit() { return new TreeSet<IntronLossJunctions>(); }

    public class IntronLossJunctions implements Comparable{

        private Set<Pair<Integer,Integer>> junctions;
        private String uniqueName;
        private GenomeLoc geneLocus;

        private IntronLossJunctions() {
            junctions = new TreeSet<Pair<Integer,Integer>>(new Comparator<Pair<Integer, Integer>>() {
                @Override
                public int compare(Pair<Integer, Integer> integerIntegerPair, Pair<Integer, Integer> integerIntegerPair1) {
                    int first = integerIntegerPair.first - integerIntegerPair1.first;
                    int second = integerIntegerPair.second - integerIntegerPair1.second;
                    return first != 0 ? first : second;
                }
            });
        }

        public IntronLossJunctions(RefSeqFeature feature) {
            uniqueName = feature.getTranscriptUniqueGeneName();
            geneLocus = feature.getLocation();
            junctions = new TreeSet<Pair<Integer,Integer>>(new Comparator<Pair<Integer, Integer>>() {
                @Override
                public int compare(Pair<Integer, Integer> integerIntegerPair, Pair<Integer, Integer> integerIntegerPair1) {
                    int first = integerIntegerPair.first - integerIntegerPair1.first;
                    int second = integerIntegerPair.second - integerIntegerPair1.second;
                    return first != 0 ? first : second;
                }
            });
        }

        public void addJunction(Pair<Integer,Integer> junction) {
            junctions.add(junction);
        }

        public int compareTo(Object other) {
            if ( other instanceof IntronLossJunctions ) {
                int locComp = this.geneLocus.compareTo(((IntronLossJunctions)other).geneLocus);
                int nameComp = this.uniqueName.compareTo(((IntronLossJunctions)other).uniqueName);
                if ( locComp == 0 && nameComp == 0 ) {
                    return 0;
                } else if ( locComp == 0 ) {
                    return nameComp;
                } else {
                    return locComp;
                }
            } else {
                return Integer.MIN_VALUE;
            }
        }

        public void add(IntronLossJunctions other) {
            if ( this.compareTo(other) != 0 ) {
                throw new ReviewedStingException("Attempting to merge junctions across different genes or gene transcripts");
            }
            this.junctions.addAll(other.junctions);
        }

        public String toString() {
            return String.format("%s\t%s\t%s\t%s", geneLocus, uniqueName, joinJunctions(), reconcile().joinJunctions());
        }

        private String joinJunctions() {
            StringBuilder buf = new StringBuilder();
            for ( Pair<Integer,Integer> j : junctions ) {
                buf.append(j.first);
                buf.append(',');
                buf.append(j.second);
                buf.append(';');
            }
            if ( buf.length() > 0 && buf.charAt(buf.length()-1) == ';') {
                buf.deleteCharAt(buf.lastIndexOf(";"));
            }
            return buf.toString();
        }

        private IntronLossJunctions reconcile() {
            IntronLossJunctions reconciled = new IntronLossJunctions();
            reconciled.geneLocus = this.geneLocus;
            reconciled.uniqueName = this.uniqueName;
            // build up an adjacency map, where each entry contains a list of
            // all entries it connects to (before and after)
            HashMap<Integer,TreeSet<Integer>> adjMap = new HashMap<Integer,TreeSet<Integer>>();
            for ( Pair<Integer,Integer> j : this.junctions ) {
                if ( ! adjMap.containsKey(j.first) ) {
                    adjMap.put(j.first,new TreeSet<Integer>());
                }
                adjMap.get(j.first).add(j.second);
                if ( ! adjMap.containsKey(j.second) ) {
                    adjMap.put(j.second, new TreeSet<Integer>());
                }
                adjMap.get(j.second).add(j.first);
            }

            Integer prevStart = null;
            Integer prevStop = null;
            HashSet<Integer> usedPositions = new HashSet<Integer>();
            for ( Pair<Integer,Integer> etry : this.junctions ) {
                Integer potentialStart = etry.first;
                // rule 1: always take the shortest (e.g. first) connection
                boolean isFirstConnection = !usedPositions.contains(etry.first);
                Integer potentialStop = etry.second;
                // rule 2: check to see if the stop position has a shorter (prior) connection
                TreeSet<Integer> stopConnections = adjMap.get(potentialStop);
                boolean targetHasShorter = false;
                for ( Integer p : stopConnections ) {
                    if ( p >= potentialStop ) {
                        break;
                    }
                    targetHasShorter |= ( p > potentialStart && p < potentialStop );
                }
                // rule 3: check to see if the junction has subordinated junctions (hooray n^2)
                boolean targetHasSubordinated = false;
                for ( int i = potentialStart+1; i < potentialStop; i++ ) {
                    if ( adjMap.containsKey(new Integer(i)) ){
                        if ( adjMap.get(new Integer(i)).first() <= potentialStop ) {
                            targetHasSubordinated = true;
                            break;
                        }
                    }
                }
                // if there are no shorter connections for the target and no subordinated junctions
                if ( ! targetHasShorter && ! targetHasSubordinated && isFirstConnection ) {
                    // do we need to fill anything in?
                    if ( prevStart == null || potentialStart.equals(prevStop) ) {
                        // no
                        reconciled.addJunction(new Pair<Integer,Integer>(potentialStart,potentialStop));
                        usedPositions.add(potentialStart);
                    } else {
                        if ( prevStop < potentialStart ) {
                            for ( int j = prevStop; j < potentialStart; j++ ) {
                                reconciled.addJunction(new Pair<Integer,Integer>(j,j+1));
                                usedPositions.add(j);
                            }
                        } else {
                            for ( int j = prevStart; j < potentialStart; j++ ) {
                                reconciled.addJunction(new Pair<Integer,Integer>(j,j+1));
                                usedPositions.add(j);
                            }
                        }
                        reconciled.addJunction(new Pair<Integer,Integer>(potentialStart,potentialStop));
                        usedPositions.add(potentialStart);
                    }

                    prevStart = potentialStart;
                    prevStop = potentialStop;
                }
            }

            return reconciled;
        }
    }

    public TreeSet<IntronLossJunctions> map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
        final List<RefSeqFeature> refSeqFeatures = new ArrayList<RefSeqFeature>(metaDataTracker.getValues(refSeqRodBinding));

        TreeSet<IntronLossJunctions> junctionsSet = new TreeSet<IntronLossJunctions>();
        GenomeLoc readLoc = getToolkit().getGenomeLocParser().createGenomeLoc(read);
        GenomeLoc mateLoc = ( (! read.getReadPairedFlag()) || read.getMateUnmappedFlag() ) ? null : getToolkit().getGenomeLocParser().createGenomeLoc(read.getMateReferenceName(),read.getMateAlignmentStart(),read.getMateAlignmentStart()+read.getReadLength());
        for ( RefSeqFeature refSeqFeature : refSeqFeatures ) {
            int readOverlapInt = refSeqFeature.getSortedOverlapInteger(readLoc);
            if ( mateLoc != null && readOverlapInt > -1 ) {
                int inferredDistance = Math.abs(read.getInferredInsertSize());
                int mateOverlapInt = refSeqFeature.getSortedOverlapInteger(mateLoc);
                Integer mateMapQ = (Integer) read.getAttribute("MQ");
                if ( mateMapQ == null ) {
                    throw new UserException.MalformedFile("Paired reads in the BAM file must have a mate mapping quality tag (MQ) filled in.");
                }
                if ( mateOverlapInt > -1 && mateOverlapInt != readOverlapInt && inferredDistance >= minInferredInsert && mateMapQ > 23 ) {
                    IntronLossJunctions j = new IntronLossJunctions(refSeqFeature);
                    if ( readOverlapInt < mateOverlapInt ) {
                        j.addJunction(new Pair<Integer,Integer>(readOverlapInt,mateOverlapInt));
                    } else {
                        j.addJunction(new Pair<Integer,Integer>(mateOverlapInt,readOverlapInt));
                    }
                    junctionsSet.add(j);
                }
            }

            int unclStart = read.getUnclippedStart();
            int unclEnd = read.getUnclippedEnd();
            if ( readOverlapInt > -1 && readLoc.getStop() > refSeqFeature.getSortedExonLoc(readOverlapInt).getStop() ) {
                for ( int k = readOverlapInt+1; k < refSeqFeature.getNumExons() ; k++ ) {
                    if ( readMatches(read,readLoc,refSeqFeature.getSortedExonLoc(k), unclStart, unclEnd) ) {
                        IntronLossJunctions j = new IntronLossJunctions(refSeqFeature);
                        j.addJunction(new Pair<Integer,Integer>(readOverlapInt,k));
                        junctionsSet.add(j);
                        break;
                    }
                }
            } else if ( readOverlapInt > -1 && readLoc.getStart() < refSeqFeature.getSortedExonLoc(readOverlapInt).getStart() ) {
                for ( int k = readOverlapInt-1; k > -1; k--) {
                    if ( readMatches(read,readLoc,refSeqFeature.getSortedExonLoc(k),unclStart,unclEnd) ) {
                        IntronLossJunctions j = new IntronLossJunctions(refSeqFeature);
                        j.addJunction(new Pair<Integer,Integer>(k,readOverlapInt));
                        junctionsSet.add(j);
                        break;
                    }
                }
            }
        }

        if ( exonSequences != null && exonSequences.size() > 50 ) {
            TreeMap<GenomeLoc,byte[]> keep = new TreeMap<GenomeLoc, byte[]>();
            for ( GenomeLoc l : exonSequences.keySet() ) {
                if ( ! readLoc.isPast(l) || readLoc.distance(l) < 10000 ) {
                    keep.put(l,exonSequences.get(l));
                }
            }
            exonSequences = keep;
        }

        return junctionsSet;
    }


    private boolean readMatches(GATKSAMRecord read, GenomeLoc readLoc, GenomeLoc exonPos, int readUnStart, int readUnEnd) {
        int match = 0;
        if ( readLoc.isPast(exonPos) ) {
            int nBases = read.getAlignmentStart()-readUnStart;
            if (nBases < 10) {
                return false;
            }
            // first bases of read should match the previous exon
            byte[] refBases = getExonBases(exonPos, nBases, false);
            //referenceReader.getSubsequenceAt(exonPos.getContig(),exonPos.getStop()-nBases,exonPos.getStop()).getBases();
            byte[] readBases = read.getReadBases();
            int max = Math.min(refBases.length,10);
            for ( int i = 0; i < max; i++ ) {
                match += (refBases[i] == readBases[i]) ? 1 : 0;
            }
        } else {
            // last bases of read should match the next exon
            int nBases = readUnEnd-read.getAlignmentEnd();
            if ( nBases < 10 ) {
                return false;
            }
            byte[] refBases = getExonBases(exonPos, nBases, true);
            //referenceReader.getSubsequenceAt(exonPos.getContig(),exonPos.getStart(),exonPos.getStart()+nBases).getBases();
            byte[] readBases = read.getReadBases();
            ArrayUtils.reverse(readBases);
            for ( int i = 0; i < 10; i++) {
                match += (refBases[i] == readBases[i]) ? 1 : 0;
            }
        }

        return match >= 8;
    }


    private byte[] getExonBases(GenomeLoc exonPosition, int numBases, boolean startOfExon) {
        if (this.exonSequences == null || ! this.exonSequences.containsKey(exonPosition) ) {
            if ( this.exonSequences == null ) {
                exonSequences = new TreeMap<GenomeLoc,byte[]>();
            }
            exonSequences.put(exonPosition,referenceReader.getSubsequenceAt(exonPosition.getContig(),
                    exonPosition.getStart(),exonPosition.getStop()).getBases());
        }

        byte[] sq = exonSequences.get(exonPosition);
        //logger.debug(String.format("%d: %d %d",sq.length,sq.length,numBases));
        return startOfExon ? Arrays.copyOfRange(sq,0,numBases) : Arrays.copyOfRange(sq,Math.max(0,sq.length-numBases),sq.length);
    }

    public TreeSet<IntronLossJunctions> reduce(TreeSet<IntronLossJunctions> pMap, TreeSet<IntronLossJunctions> pRed) {
        for ( IntronLossJunctions ilj : pMap ) {
            if ( pRed.contains(ilj) ) {
                pRed.tailSet(ilj).first().add(ilj);
            } else {
                pRed.add(ilj);
            }
        }

        return pRed;
    }

    public void onTraversalDone(TreeSet<IntronLossJunctions> junctions) {
        out.printf("HEADER%s\t%s\t%s\t%s%n","Loc","Transcript","Junctions","Potential Event");
        for ( IntronLossJunctions ilj : junctions ) {
            out.printf("%s%n",ilj.toString());
        }
    }

    public void initialize() {
        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile,ex);
        }
    }
}
