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

package org.broadinstitute.sting.walkers.varianteval;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEval;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.utils.report.tags.Analysis;
import org.broadinstitute.sting.utils.report.tags.DataPoint;
import org.broadinstitute.sting.utils.report.utils.TableType;

import java.util.Collection;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Nov 22, 2010
 * Time: 12:22:08 PM
 * To change this template use File | Settings | File Templates.
 */
@Analysis(name = "ACTransitionMatrix", description = "Number of additional genotypes from each new sample; random permutations")
public class ACTransitionTable extends VariantEvaluator {
    private final int NUM_PERMUTATIONS = 50;
    private final double LOW_GQ_PCT = 0.95;
    private final double LOW_GQ_THRSH = 30.0;
    private boolean initialized = false;
    private long skipped = 0l;

    @DataPoint(name="Het transitions",description="AC[s] = AC[s-1]+1 and AC[s] = AC[s-1]+2 transitions")
    TransitionTable transitions = null;
    @DataPoint(name="Private permutations",description="Marginal increase in number of sites per sample")
    PermutationCounts privatePermutations;
    @DataPoint(name="AC2 Permutations",description="Marginal increase in number of AC=2 sites, per sample")
    PermutationCounts doubletonPermutations;
    @DataPoint(name="AC3 Permutations",description="Marginal increase in number of tripleton sites, per sample")
    PermutationCounts tripletonPermutations;

    String[][] permutations;

    public boolean enabled() {
        return true;
    }

    public int getComparisonOrder() {
        return 2;
    }

    public String getName() {
        return "ACTransitionTable";
    }

    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( eval != null && ! initialized ) {
            //this.veWalker.getLogger().warn("Initializing...");
            initialize(eval);
            initialized = true;
        }

        if ( isGood(eval) ) {
            if ( comp != null && ! comp.isFiltered() ) {
                return null;
            }

            int order_offset = 0;
            for ( String[] ordering : permutations ) {
                int sample_offset = 0;
                int variant_ac = 0;
                for ( String sample : ordering ) {
                    if ( eval.getGenotype(sample).isHet() ) {
                        variant_ac++;
                        transitions.hetTransitionCounts[order_offset][variant_ac-1][sample_offset]++;
                    } else if ( eval.getGenotype(sample).isHomVar() ) {
                        variant_ac += 2;
                        transitions.homTransitionCounts[order_offset][variant_ac-1][sample_offset]++;
                    } else {
                        // todo -- note, unclear how to treat no calls. Is the hom in het,ref,ref,nocall,hom sample 4 or 5?
                        // todo -- do we want to tabulate P[sample i is not variant | some variant]? This is just combinatorics so i left it out
                        if ( variant_ac > 0 ) {
                            transitions.stationaryCounts[order_offset][variant_ac-1][sample_offset]++;
                        }
                    }
                    sample_offset ++;
                }
                order_offset++;
            }
        } else {
            skipped++;    
        }

        return null;
    }

    private boolean isGood(VariantContext vc) {
        if ( vc == null || vc.isFiltered() || (vc.getHetCount() + vc.getHomVarCount() == 0) ) { // todo -- should be is variant, but need to ensure no alt alleles at ref sites
            return false;
        } else {
            Collection<Genotype> gtypes = vc.getGenotypes().values();
            int ngood = 0;
            for ( Genotype g : gtypes) {
                if ( g.isCalled() && g.getPhredScaledQual() >= LOW_GQ_THRSH ) {
                    ngood ++;
                }
            }

            return ( (0.0+ngood)/(0.0+gtypes.size()) >= LOW_GQ_PCT );
        }
    }

    public ACTransitionTable(VariantEval parent) {
        //super(parent);
    }

    public void initialize(VariantContext vc) {
        Set<String> permuteSamples = vc.getSampleNames();
        permutations = new String[NUM_PERMUTATIONS][permuteSamples.size()];
        //veWalker.getLogger().warn(String.format("Num samples: %d",permuteSamples.size()));
        int offset = 0;
        for ( String s : permuteSamples ) {
            permutations[0][offset] = s;
            offset ++;
        }
        
        for ( int p = 1; p < NUM_PERMUTATIONS ; p++ ) {
            permutations[p] = permutations[0].clone();
            for ( int o = 0; o < permutations[p].length; o ++ ) {
                int r = (int) Math.floor(Math.random()*(o+1));
                String swap = permutations[p][r];
                permutations[p][r] = permutations[p][o];
                permutations[p][o] = swap;
            }
        }

        transitions = new TransitionTable();
        transitions.hetTransitionCounts = new int[NUM_PERMUTATIONS][permuteSamples.size()*2][permuteSamples.size()];
        transitions.homTransitionCounts = new int[NUM_PERMUTATIONS][permuteSamples.size()*2][permuteSamples.size()];
        transitions.stationaryCounts = new int[NUM_PERMUTATIONS][permuteSamples.size()*2][permuteSamples.size()];
        privatePermutations = new PermutationCounts(1,transitions);
        doubletonPermutations = new PermutationCounts(2,transitions);
        tripletonPermutations = new PermutationCounts(3,transitions);
    }

    public void finalizeEvaluation() { // note: data points are null when this is called (wtf?)
        //veWalker.getLogger().info(String.format("Skipped: %d",skipped));
    }

    class TransitionTable implements TableType {
        int[][][] hetTransitionCounts;
        int[][][] homTransitionCounts;
        int[][][] stationaryCounts;
        String[][] countAverages;
        String[] rowKeys = null;
        String[] colKeys = null;

        public Object[] getRowKeys() {
            if ( rowKeys == null ) {
                rowKeys = new String[3*hetTransitionCounts[0].length];
                for ( int i = 0; i < hetTransitionCounts[0].length; i ++ ) {
                    rowKeys[i] = String.format("%s%d%s","AC_",i,"_(het)");
                }
                for ( int i = 0; i < hetTransitionCounts[0].length; i ++ ) {
                    rowKeys[hetTransitionCounts[0].length+i] = String.format("%s%d%s","AC_",i,"_(hom)");
                }
                for ( int i = 0; i < hetTransitionCounts[0].length; i ++ ) {
                    rowKeys[2*hetTransitionCounts[0].length+i] = String.format("%s%d%s","AC_",i,"_(ref)");
                }
            }


            return rowKeys;
        }

        public String getCell(int x, int y) {
            if ( countAverages == null ) {
                countAverages = new String[hetTransitionCounts[0].length*3][hetTransitionCounts[0][0].length];
                for ( int sam = 0; sam < hetTransitionCounts[0][0].length; sam ++) {
                    for ( int idx = 0 ; idx < hetTransitionCounts[0].length; idx ++ ) {
                        int totalTimesAtACSample = 0;
                        int totalStationary = 0;
                        int totalAC1Shift = 0;
                        int totalAC2Shift = 0;
                        for ( int p = 0; p < hetTransitionCounts.length; p++ ) {
                            totalStationary += stationaryCounts[p][idx][sam];
                            totalAC2Shift += (idx+2 >= hetTransitionCounts[0][0].length) ? 0 : homTransitionCounts[p][idx+2][sam];
                            totalAC1Shift += (idx+1 >= hetTransitionCounts[0][0].length) ? 0 : hetTransitionCounts[p][idx+1][sam];
                        }
                        totalTimesAtACSample = totalStationary+totalAC1Shift+totalAC2Shift;
                        countAverages[idx][sam] = formatProp(totalAC1Shift,totalTimesAtACSample);
                        countAverages[hetTransitionCounts[0].length+idx][sam] = formatProp(totalAC2Shift,totalTimesAtACSample);
                        countAverages[hetTransitionCounts[0].length*2+idx][sam] = formatProp(totalStationary,totalTimesAtACSample);
                    }
                }
            }

            return countAverages[x][y] == null ? "0.00" : countAverages[x][y];
        }

        private String formatProp(int num, int denom) {
            return (denom != 0) ? String.format("%.4f", ((double) num)/denom) : "0.0";
        }

        public String getName() { return "AC Transition Tables"; }

        public Object[] getColumnKeys() {
            if ( colKeys == null ) {
                colKeys = new String[hetTransitionCounts[0][0].length];
                for ( int ac = 0; ac < hetTransitionCounts[0][0].length; ac ++ ) {
                    colKeys[ac] = String.format("Sample_%d",ac);
                }
            }

            return colKeys;
        }
    }


    class PermutationCounts implements TableType {
        int acToExtract;
        TransitionTable table;
        String[] rowNames;
        String[] colNames;

        public PermutationCounts(int ac, TransitionTable tTable) {
            acToExtract = ac;
            table = tTable;
        }

        public String[] getRowKeys() {
            //System.out.printf("%s%n",table);
            if ( rowNames == null ) {
                rowNames = new String[table.stationaryCounts.length];
                for ( int p = 0 ; p < rowNames.length; p ++ ) {
                    rowNames[p] = String.format("Perm%d",p+1);
                }
            }

            return rowNames;
        }

        public String[] getColumnKeys() {
            if ( colNames == null ) {
                colNames = new String[table.stationaryCounts[0][0].length];
                for ( int s = 0 ; s < colNames.length; s ++ ) {
                    colNames[s] = String.format("Sample%d",s+1);
                }
            }

            return colNames;
        }

        public Integer getCell(int x, int y) {
            return table.hetTransitionCounts[x][acToExtract-1][y] +
                    ( (acToExtract > table.homTransitionCounts[0][0].length) ? 0 : table.homTransitionCounts[x][acToExtract-1][y]);
        }

        public String getName() {
            return String.format("PermutationCountsAC%d",acToExtract);
        }

        public void init() {
            getRowKeys();
            getColumnKeys();
            getCell(1,1);
        }
    }


}

