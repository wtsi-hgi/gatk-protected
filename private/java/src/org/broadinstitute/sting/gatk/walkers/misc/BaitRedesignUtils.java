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

package org.broadinstitute.sting.gatk.walkers.misc;

import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.alignment.bwa.c.BWACAligner;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/3/12
 * Time: 10:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class BaitRedesignUtils {

    public static final double OPT_GC = 0.47;
    private static final double MIXT_COEF = 0.11;
    private static final double ONE_MINUS_MC = 1.0-MIXT_COEF;
    private static final int POP_SIZE = 100;
    private static final int SELECT_SIZE = 20;
    private static final int GENERATIONS = 50;
    private static final double INIT_MUTATION_RATE = 0.10;
    private static double MUTATION_RATE = INIT_MUTATION_RATE;

    private static double optFunction(byte[] seq, BWACAligner aligner, GenomeLoc initialPosition) {
        // count up possible alignments
        int altMapping = 0;
        boolean found = false;
        for (Alignment[] aliSet : aligner.getAllAlignments(seq)) {
            for ( Alignment a : aliSet ) {
                ++altMapping;
                if (! found && a.getContigIndex() == initialPosition.getContigIndex() &&
                        Math.abs(a.getAlignmentStart()-initialPosition.getStart()) < 50 ) {
                    found = true;
                }
            }
        }
        if ( ! found ) {
            // doesn't align anywhere. Not so good.
            return Double.POSITIVE_INFINITY;
        } else {
            --altMapping; // should be one mapping: the current position
        }

        // todo -- don't need to recalculate GC every time, but can cache it and alter it when anything changes
        return MIXT_COEF*altMapping + ONE_MINUS_MC*Math.pow(OPT_GC-calculateGC(seq),2);
    }

    public static double optFunction(byte[] seq, byte[] init) {
        return MIXT_COEF*editDist(seq,init)+ONE_MINUS_MC*Math.pow(OPT_GC-calculateGC(seq),2);
    }

    public static double editDist(byte[] a, byte[] b) {
        int dist = 0;
        for ( int i = 0; i < a.length; i++ ) {
            dist += b[i] == a[i] ? 0 : 1;
        }

        return ((double)dist)/a.length;
    }

    private static Comparator<byte[]> getSelectionComparator(final BWACAligner aligner, final GenomeLoc position) {
        return new Comparator<byte[]>() {
            @Override
            public int compare(byte[] bytes, byte[] bytes1) {
                int d  = Double.compare(optFunction(bytes,aligner,position),optFunction(bytes1,aligner,position));
                if ( d == 0 ) {
                    return match(bytes,bytes1) ? 0 : GenomeAnalysisEngine.getRandomGenerator().nextInt(1)-1;
                }
                return d;
            }
        };
    }

    private static Comparator<byte[]> getSelectionComparator(final byte[] init) {
        return new Comparator<byte[]>() {
            @Override
            public int compare(byte[] bytes, byte[] bytes1) {
                int d  = Double.compare(optFunction(bytes,init),optFunction(bytes1,init));
                if ( d == 0 ) {
                    return match(bytes,bytes1) ? 0 : GenomeAnalysisEngine.getRandomGenerator().nextInt(1)-1;
                }
                return d;
            }
        };
    }

    private static boolean match(byte[] bytes, byte[] bytes1) {
        return Arrays.equals(bytes,bytes1);
    }

    public static double calculateGC(byte[] seq) {
        int tot = 0;
        int gc = 0;
        for ( byte b : seq ) {
            switch (b) {
                case 'C':
                case 'G':
                    ++gc;
                default:
                    ++tot;
            }
        }

        return ((double) gc)/tot;
    }

    public static byte[] getOptimalBases(BWACAligner aligner, byte[] initialSequence, GenomeLoc position) {
        // this is a quasi-genetic algorithm designed to move the initial sequence towards the optimal GC point
        // optimization function to alignability and GC
        List<byte[]> parents = new ArrayList<byte[]>();
        parents.add(initialSequence.clone());
        MUTATION_RATE = INIT_MUTATION_RATE;
        int generation = 0;
        do {
            List<byte[]> children = reproduce(parents);
            parents = select(children,initialSequence);
            ++generation;
            MUTATION_RATE *= 0.95;
        } while ( generation < GENERATIONS );

        return parents.get(0);
    }

    private static List<byte[]> select(List<byte[]> population, BWACAligner aligner, GenomeLoc position) {
        TreeSet<byte[]> selectionSet = new TreeSet<byte[]>(getSelectionComparator(aligner,position));
        selectionSet.addAll(population);
        List<byte[]> selected = new ArrayList<byte[]>(SELECT_SIZE);
        int num = 0;
        // todo -- when clear this is properly working, break early
        for ( byte[] child : selectionSet ) {
            selected.add(child);
            if ( ++num >= SELECT_SIZE ) {
                break;
            }
        }
        return selected;
    }

    private static List<byte[]> select(List<byte[]> population, byte[] init) {
        TreeSet<byte[]> selectionSet = new TreeSet<byte[]>(getSelectionComparator(init));
        selectionSet.addAll(population);
        List<byte[]> selected = new ArrayList<byte[]>(SELECT_SIZE);
        int num = 0;
        // todo -- when clear this is properly working, break early
        for ( byte[] child : selectionSet ) {
            selected.add(child);
            if ( ++num >= SELECT_SIZE ) {
                break;
            }
        }
        return selected;
    }

    private static List<byte[]> reproduce(List<byte[]> parents) {
        List<byte[]> children = new ArrayList<byte[]>(POP_SIZE+parents.size());
        // perform recombination to generate a population
        while ( children.size() < POP_SIZE ) {
            byte[] p1 = parents.get(GenomeAnalysisEngine.getRandomGenerator().nextInt(parents.size()));
            byte[] p2 = parents.get(GenomeAnalysisEngine.getRandomGenerator().nextInt(parents.size()));
            if ( p1 == p2 ) {
                children.add(mutate(p1.clone()));
            } else {
                // recombine
                int recomLoc = GenomeAnalysisEngine.getRandomGenerator().nextInt(p1.length);
                byte[] child = p1.clone();
                for ( int i = recomLoc; i < child.length; i++) {
                    child[i] = p2[i];
                }
                // mutate
                children.add(mutate(child));
            }
        }

        for ( byte[] b : parents ) {
            children.add(b.clone());
        }

        return children;
    }

    private static byte[] mutate(byte[] offspring) {
        // every base has a chance to mutate to another base, but biased to increase or decrease GC
        // note that A -> {C,G,T}, so on average we expect GC to go up
        // by symmetry, uniform swapping tends to drive GC -> 50%, and additional bias need not factor in
        for ( int i = 0; i < offspring.length; i++ ) {
            if ( GenomeAnalysisEngine.getRandomGenerator().nextDouble() < MUTATION_RATE ) {
                // mutating this base
                offspring[i] = BaseUtils.baseIndexToSimpleBase(
                        BaseUtils.getRandomBaseIndex(
                                BaseUtils.simpleBaseToBaseIndex(offspring[i])));
            }
        }

        return offspring;
    }
}
