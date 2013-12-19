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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.KmerSequence;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.KmerSequenceGraphMap;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.ReadGraphMap;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.junit.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: valentin
 * Date: 7/25/13
 * Time: 4:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class MLCalculatorsUnitTest extends BaseTest {
    protected final MLLog10PairHMM exactMLE = new MLLog10PairHMM((byte)10);
    protected final MLLog10PairHMM threeDMLE = new MLLog10PairHMM((byte)10);
    public MLCalculatorsUnitTest() {
        exactMLE.initialize(200,400);
        exactMLE.setupForGraphLikelihoodTest();
    }

    //TODO fix this test and move unused component to archive.
    @Test(enabled=false,dataProvider="MLEComparisonDP")
    public void testMLEstimate(final ReadThreadingGraph graph, final String haplotype, final String bases, byte[] bq, byte[] iq, byte[] dp, final String info)  {

        final byte[] overalGCP = new byte[bases.length()];
        Arrays.fill(overalGCP,(byte)10);

        KmerSequenceGraphMap haplotypeMap = new KmerSequenceGraphMap(
                graph,new KmerSequence(haplotype.getBytes(), graph.getKmerSize()));

        //final List<Kmer> missingKmers = haplotypeMap.missingKmers();
        //if (missingKmers.size() > 0) {
        //    throw new IllegalArgumentException("unexpected haplotype missing kmers!!! " + missingKmers.size() + ": " + Utils.join(",", missingKmers) );
        //}



        final ReadGraphMap readMapSingleArray = new ReadGraphMap(graph,bases.getBytes(),bq,iq,dp,100);
        readMapSingleArray.setUseSingleArrayMLImplementation(true);
        final ReadGraphMap readMapMultiArray = new ReadGraphMap(graph,bases.getBytes(),bq,iq,dp,100);
        readMapMultiArray.setUseSingleArrayMLImplementation(false);
        final MLLog10PairHMM exactMLE = new MLLog10PairHMM((byte)10);
        exactMLE.setupForGraphLikelihoodTest();
        if ((bases.equals("CCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTTTCATATATGGG"))
                && haplotype.equals("TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTGTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT")) {
        }
        double graphLn = readMapSingleArray.log10MLE(haplotypeMap);
        double graphLn2 = readMapMultiArray.log10MLE(haplotypeMap);
        exactMLE.initialize(2000,4000);
        double exactLn = exactMLE.computeReadLikelihoodGivenHaplotypeLog10(haplotype.getBytes(),bases.getBytes(),bq,iq,dp,overalGCP,false,null);
        Assert.assertEquals(graphLn, exactLn,0.00001);
        Assert.assertEquals(graphLn,graphLn2,0.0000000001);
    }

    @DataProvider(name="MLEComparisonDP")
    public Iterator<Object[]> mleComparisonDP() {
        List<Object[]> result = new LinkedList<>();
        for (final Object[] tsData : TEST_SETS) {
           final ActiveRegionTestDataSet ts = toTestSet(tsData);
           if (ts == null) {
               continue;
           }
           final ReadThreadingGraph rtg = buildGraphFromHaplotypes(ts);
           final List<String> haplotypes = ts.haplotypesStrings();
           final List<String> reads = ts.readStrings();
           int hIdx = 0;
           for (final String h : haplotypes) {
               int idx = 0;
               for (final String r : reads) {
                   String info = ts.haplotypeCigars[hIdx] + " <> " + ts.readCigars[idx];
                   final Object[] problem = new Object[] { rtg, h, r, Arrays.copyOf(ts.bq,r.length()),
                           Arrays.copyOf(ts.iq,r.length()), Arrays.copyOf(ts.dq, r.length()), info  };
                   result.add(problem);
                   idx++;
               }
               hIdx++;
           }
        }
        return result.iterator();
    }

    protected ReadThreadingGraph buildGraphFromHaplotypes(final ActiveRegionTestDataSet ts) {

        final ReadThreadingGraph rtg = new ReadThreadingGraph(ts.kmerSize());
        rtg.addSequence("REF",ts.getReference().getBytes(),null,true);
        int hapIndex = 1;
        List<String> haplotypes = ts.haplotypesStrings();
        for (final String hap : haplotypes) {
            rtg.addSequence("H" + hapIndex++,hap.getBytes(),null,false);
        }
        rtg.buildGraphIfNecessary();

        return rtg;
    }



    public ActiveRegionTestDataSet toTestSet(final Object[] params) {

        boolean enabled = (boolean) params[0];
        if (!enabled) {
            return null;
        }

        return new ActiveRegionTestDataSet((Integer)params[1], (String)params[2], (String[]) params[3],
                (String[]) params[4], (byte[]) params[5], (byte[]) params[6], (byte[]) params[7]);
    }



    protected Object[][] TEST_SETS = new Object[][] {
            new Object[] {
                    true,//disabled
                    Integer.valueOf(4),
             /* ref 200M*/
                    "GAAG",
                    new String[] {
                            "4="
                    },
                         /* readCigars */
                    new String[] {
                            "0:0:3=1T"
                    },

                    rep((byte)30,1000),
                    rep((byte)60,1000),
                    rep((byte)60,1000)


                    },
            new Object[] {
                    true,
                    Integer.valueOf(2),
             /* ref 200M*/
                    "GCATC",
                    new String[] {
                            "2=1T2=",
                            "2=1D2=",
                            "2=1Ic3=",
                            "4=1V"
                    },
                         /* readCigars */
                    new String[] {
                            "0:0:5="
                    },

                    rep((byte)30,1000),
                    rep((byte)21,1000),
                    rep((byte)21,1000)
            },
            new Object[] {
                    true,
                    Integer.valueOf(11),
             /* ref 200M*/
                    "GAAGATTAGAGAAAAAAGAATAAAAAGAAACGAACAAAGCCTCCAAGAAATATGGGACTA" +
                    "TGTGAAAAGACCAAATCTACGTCTGATTGGTGTACCTGAAAGTGACGGGGAGAATGGAAC" +
                    "CAAGTTGAAAAACACTCTGCAGGATATTATCCAGGAGAACTCCCCCAATCTAGCAAGGCA" +
                    "GGCCAACATTCAGATTCAGG",
                          /* haps */
                    new String[] {
                            "99=1T100=",
                            "50=1V149="
                    },
                         /* readCigars Descriptors, Allele:Offset:Cigar+ */
                    new String[] {
                            "0:30:80=",
                            "0:10:19=1D70=",
                            "0:20:29=1V50=",
                            "1:80:80=",
                            "2:20:20=",
                            "0:10:2V10=2V4D50=",
                            "0:12:5Itcaga30=1V2T4W25=1D",
                            "0:0:5Itcgag30=5Itcgag",
                    },

                    rep((byte)30,1000),
                    rep((byte)60,1000),
                    rep((byte)60,1000)


                    } ,

            new Object[] {
                  true,
                  11,
                  "TAGTGGCGGGCACCTGTAATCCCAGCTACTCGGGAGGCTGAGACAGGAGAATGGCGTGAACCCGGGAGGCGGAGCCTGCAGTGAGCCGAGATAGCGCCCCTGCACTCCAGCCTGGATGACTGAACGAGACCGTCTCAAAAAAAAAAATAAAATAAAAACCAGTCCTGTATTTGAATGGCTAAGATTATGGCTTGAAATATAAAAGAAAAATGCTGACA",
                  new String[] {
                          "TAGTGGCGGGCACCTGTAATCCCAGCTACTCGGGAGGCTGAGACAGGAGAATGGCGTGAACCCGGGAGGCGGAGCCTGCAGTGAGCCGAGATAGCGCCCCTGCACTCCAGCCTGGATGACTGAACGAGACCGTCTCAAAAAAAAAAAAAAAAAAAAAACCAGTCCTGTATTTGAATGGCTAAGATTATGGCTTGAAATATAAAAGAAAAATGCTGACA",
                          "TAGTGGCGGGCACCTGTAATCCCAGCTACTCGGGAGGCTGAGACAGGAGAATGGCGTGAACCCGGGAGGCGGAGCCTGCAGTGAGCCGAGATAGCGCCCCTGCACTCCAGCCTGGATGACTGAACGAGACCGTCTCAAAAAAAAAAAAAACAAAACAAAAACCAGTCCTGTATTTGAATGGCTAAGATTATGGCTTGAAATATAAAAGAAAAATGCTGACA",
                          "TAGTGGCGGGCACCTGTAATCCCAGCTACTCGGGAGGCTGAGACAGGAGAATGGCGTGAACCTGGGAGGCGGAGCCTGCAGTGAGCCGAGATAGCGCCCCTGCACTCCAGCCTGGATGACTGAACGAGACCGTCTCAAAAAAAAAAAAAAAAAAAAAACCAGTCCTGTATTTGAATGGCTAAGATTATGGCTTGAAATATAAAAGAAAAATGCTGACA",
                          "TAGTGGCGGGCACCTGTAATCCCAGCTACTCGGGAGGCTGAGACAGGAGAATGGCGTGAACCTGGGAGGCGGAGCCTGCAGTGAGCCGAGATAGCGCCCCTGCACTCCAGCCTGGATGACTGAACGAGACCGTCTCAAAAAAAAAAAAAACAAAACAAAAACCAGTCCTGTATTTGAATGGCTAAGATTATGGCTTGAAATATAAAAGAAAAATGCTGACA",
                          "TAGTGGCGGGCACCTGTAATCCCAGCTACTCGGGAGGCTGAGACAGGAGAATGGCGTGAACCTGGGAGGCGGAGCCTGCAGTGAGCCGAGATAGCGCCCCTGCACTCCAGCCTGGATGACTGAACGAGACCGTCTCAAAAAAAAAAATAAAATAAAAACCAGTCCTGTATTTGAATGGCTAAGATTATGGCTTGAAATATAAAAGAAAAATGCTGACA",
                  },
                  new String[] {
                          "AAAAAAAAAAAAAACAAAACAAAAACCAGTCCTGTATTTGAATGGCTAAGATTATGGCTTGAAATATAAAAGAAAAATGCTGACA"
                  },
                  rep((byte)30,1000),
                  rep((byte)45,1000),
                  rep((byte)45,1000)
            },

            new Object[] {
                 true,
                 11,
                 "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                 new String[] {
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTGCAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTGCAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTGCAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTGCAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTGCAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCATGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTGCAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCGTGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCGTGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCGTGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCGTGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCGTGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCGTGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTACAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGGGAGGGAGGGAGGGAAGAAGGAAGGAAGGAAGGAAGGGCAGAAAGCAGACCAGACCAGATTTCATAGCTGTATTCTGTCT",
                 } ,
                 new String[] {
                         "CAAGAGTTTGGAAAGATGATTAAAAATGTACCCTCTAAAGAGCAAGCTGGGCATGGTGGCTCGTGCCTGTAGCCCTAGCTACATGGGTGGCTGAGGCAGAAAGATCACTTGAGCCCAGAAGTCCAGAAGTTCAAGGCTGCAGTGAGCTATGATTGTGCCAGTGCCCTCCAGAAGAGAGGAAGAAAGAAAGAGAG"
                 },
                 rep((byte)30,1000),
                 rep((byte)20,1000),
                 rep((byte)20,1000),
            },

            new Object[] {
                true,
                25,
                "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTCTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                new String[] {
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTGTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTGTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTCTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTTTAAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTGTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTGTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTCTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTCTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTTTTTTTTTTTAAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACAATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                        "TAAAGAGAAGCAGGTTTTTAACTGGGTTTTTATTTCTGTTGTGTCTCTAGCCAAAGCCATGCTGATGTTTGATATTTTTTTTTTTTTTTTTTTTTTTTAAGAGATGGGGTTTCACCATGTTGGCCATGCTGGTCGTGAACTCCTGACCTCAATGGCCTCCCAAACTGTTGGGAT",
                },
                new String[] {
                        "CCATGCTGATGTTTGATAGTTTTTTTTTTTTTTTTTTTTTCATATATGGG"
                },
                rep((byte)30,1000),
                rep((byte)45,1000),
                rep((byte)45,1000)
            },
            new Object[] {
                    true,
                    11,
                    "AATAAATGTTTGACTATTAATGCCTGAGAACGGAAGGTGATTATTAATGAGATGAAAAAGTTAATCAGATTCTCCAAGTTAGGAGGGACTTGAAGACCAAATTGATAAAAATAAAAAAAAAGATGTCATAGTAGAATAATCTAGATAATAAGCAATCAATGAGACTGAAAAAATAAAATCAAGTATA",
                    new String[] {"AATAAATGTTTGACTATTAATGCCTGAGAACGGAAGGTGATTATTAATGAGATGAAAAAGTTAATCAGATTCTCCAAGTTAGGAGGGACTTGAAGACCAAATTGATAAAAATAAAAAAAAAAAGATGTCATAGTAGAATAATCTAGATAATAAGCAATCAATGAGACTGAAAAAATAAAATCAAGTATA",
                            "AATAAATGTTTGACTATTAATGCCTGAGAACGGAAGGTGATTATTAATGAGATGAAAAAGTTAATCAGATTCTCCAAGTTAGGAGGGACTTGAAGACCAAATTGATAAAAATAAAAAAAAAAGATGTCATAGTAGAATAATCTAGATAATAAGCAATCAATGAGACTGAAAAAATAAAATCAAGTATA",
                            "AATAAATGTTTGACTATTAATGCCTGAGAACGGAAGGTGATTATTAATGAGATGAAAAAGTTCATCAGATTCTCCAAGTTAGGAGGGACTTGAAGACCAAATTGATAAAAATAAAAAAAAAAAGATGTCATAGTAGAATAATCTAGATAATAAGCAATCAATGAGACTGAAAAAATAAAATCAAGTATA",
                            "AATAAATGTTTGACTATTAATGCCTGAGAACGGAAGGTGATTATTAATGAGATGAAAAAGTTCATCAGATTCTCCAAGTTAGGAGGGACTTGAAGACCAAATTGATAAAAATAAAAAAAAAAGATGTCATAGTAGAATAATCTAGATAATAAGCAATCAATGAGACTGAAAAAATAAAATCAAGTATA",
                            "AATAAATGTTTGACTATTAATGCCTGAGAACGGAAGGTGATTATTAATGAGATGAAAAAGTTCATCAGATTCTCCAAGTTAGGAGGGACTTGAAGACCAAATTGATAAAAATAAAAAAAAAGATGTCATAGTAGAATAATCTAGATAATAAGCAATCAATGAGACTGAAAAAATAAAATCAAGTATA"}
                    , new String[] { "TGAAGACCAAATTGAAAAAAATAAAAAAAAAAGATGTCATAGTAGAATAATCTAGATAATAAGCAATCAATGAGACTGAAAAAATAAAATCAAGTATA"} ,
                    rep((byte)30,1000),
                    rep((byte)45,1000),
                    rep((byte)45,1000)
            }



    };

    private byte[] rep(byte v, int length) {
        final byte[] result = new byte[length];
        Arrays.fill(result, v);
        return result;
    }



}

