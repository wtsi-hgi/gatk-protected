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

package org.broadinstitute.sting.gatk.walkers.HLAcaller;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;

import java.io.PrintStream;
import java.util.*;

/**
 * Calculates likelihood of observing the data given pairs of HLA alleles. NOTE: run CalculateBaseLikelihoods first! Usage: java -jar GenomeAnalysisTK.jar -T CalculateAlleleLikelihoods -I /humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.nuc.imputed.4digit.bam -R /broad/1KG/reference/human_b36_both.fasta -L /humgen/gsa-scr1/GSA/sjia/454_HLA/HAPMAP270/HLA_exons.interval -bl INPUT.baselikelihoods -eth\
nicity Caucasian | grep -v "INFO"  | grep -v "DONE!" > OUTPUT.allelelikelihoods
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CalculateAlleleLikelihoodsWalker extends ReadWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @Argument(fullName = "baseLikelihoods", shortName = "bl", doc = "Base likelihoods file", required = true)
    public String baseLikelihoodsFile = "";

    @Argument(fullName = "debugHLA", shortName = "debugHLA", doc = "Print debug", required = false)
    public boolean DEBUG = false;

    @Argument(fullName = "debugAlleles", shortName = "debugAlleles", doc = "Print likelihood scores for these alleles", required = false)
    public String debugAlleles = "";

    @Argument(fullName = "onlyfrequent", shortName = "onlyfrequent", doc = "Only consider alleles with frequency > 0.0001", required = false)
    public boolean FREQUENT = false;

    @Argument(fullName = "ethnicity", shortName = "ethnicity", doc = "Use allele frequencies for this ethnic group", required = false)
    public String ethnicity = "CaucasianUSA";

    String CaucasianAlleleFrequencyFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_Caucasians.freq";
    String BlackAlleleFrequencyFile     = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_BlackUSA.freq";
    String AlleleFrequencyFile          = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_CaucasiansUSA.freq";
    String UniqueAllelesFile            = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/UniqueAlleles4Digit";
    String HLAdatabaseFile              = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_DICTIONARY.txt";
    String HLA2DigitFile              = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_DICTIONARY_2DIGIT.txt";

    Hashtable AlleleFrequencies,UniqueAlleles,Alleles2Digit;
    
    CigarParser formatter = new CigarParser();
    double[][] baseLikelihoods;
    int[] positions;
    boolean loaded = false;

    String[] HLAnames, HLAreads, HLAnames2, HLAreads2;
    Integer[] HLAstartpos, HLAstoppos, HLAstartpos2, HLAstoppos2;
    ArrayList<String> HLAnamesAL, HLAreadsAL, Loci, AllelesToSearch;
    ArrayList<Integer> HLAstartposAL, HLAstopposAL;

    public Integer reduceInit() {
        if (!loaded){
            loaded = true;
            BaseLikelihoodsFileReader baseLikelihoodsReader = new BaseLikelihoodsFileReader();
            baseLikelihoodsReader.ReadFile(baseLikelihoodsFile, false);
            baseLikelihoods = baseLikelihoodsReader.GetBaseLikelihoods();
            positions = baseLikelihoodsReader.GetPositions();

            HLAnamesAL = new ArrayList<String>();
            HLAreadsAL = new ArrayList<String>();
            HLAstartposAL = new ArrayList<Integer>();
            HLAstopposAL = new ArrayList<Integer>();

            out.printf("INFO  Reading HLA alleles ... ");
            HLAFileReader HLADictionaryReader = new HLAFileReader();
            HLADictionaryReader.ReadFile(HLAdatabaseFile);
            HLAreads = HLADictionaryReader.GetSequences();
            HLAnames = HLADictionaryReader.GetNames();
            HLAstartpos = HLADictionaryReader.GetStartPositions();
            HLAstoppos = HLADictionaryReader.GetStopPositions();

            HLADictionaryReader = new HLAFileReader();
            HLADictionaryReader.ReadFile(HLA2DigitFile);
            HLAreads2 = HLADictionaryReader.GetSequences();
            HLAnames2 = HLADictionaryReader.GetNames();
            HLAstartpos2 = HLADictionaryReader.GetStartPositions();
            HLAstoppos2 = HLADictionaryReader.GetStopPositions();
            out.printf("Done! %s HLA alleles loaded.\n",HLAreads.length);

            //out.printf("INFO Common alleles:\n");
            for (int i = 1; i < UniqueAlleles.size(); i++){
                //out.printf("INFO %s\n",UniqueAlleles.values().toArray()[i]);
            }
            //out.printf("INFO  Reading HLA dictionary ...");


        }
        return 0;
    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        //HLAnamesAL.add(read.getReadName());
        //HLAreadsAL.add(formatter.FormatRead(read.getCigarString(), read.getReadString()));
        //HLAstartposAL.add(read.getAlignmentStart());
        //HLAstopposAL.add(read.getAlignmentEnd());
        //out.printf("INFO\t%s\t%s\t%s\t%s\n",read.getReadName(),read.getAlignmentStart(),read.getAlignmentEnd(),formatter.FormatRead(read.getCigarString(), read.getReadString()));
        return 1;
    }

    private int GenotypeIndex(char a, char b){
        switch(a){
            case 'A':
                switch(b){
                    case 'A': return 0;
                    case 'C': return 1;
                    case 'G': return 2;
                    case 'T': return 3;
                };
            case 'C':
                switch(b){
                    case 'A': return 1;
                    case 'C': return 4;
                    case 'G': return 5;
                    case 'T': return 6;
                };
            case 'G':
                switch(b){
                    case 'A': return 2;
                    case 'C': return 5;
                    case 'G': return 7;
                    case 'T': return 8;
                };
            case 'T':
                switch(b){
                    case 'A': return 3;
                    case 'C': return 6;
                    case 'G': return 8;
                    case 'T': return 9;
                };
            default: return -1;
        }
    }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer numreads) {
        //out.printf("Done! %s alleles found\n", numreads);
        //HLAnames = HLAnamesAL.toArray(new String[numreads]);
        //HLAreads = HLAreadsAL.toArray(new String[numreads]);
        //HLAstartpos = HLAstartposAL.toArray(new Integer[numreads]);
        //HLAstoppos = HLAstopposAL.toArray(new Integer[numreads]);

        double[][] AlleleLikelihoods = new double[numreads][numreads];

        String name1, name2;
        double frq1, frq2;


        double minfrq = 0;
        if (FREQUENT){
            minfrq = 0.0001;
        }
        int numcombinations = 0;
        out.printf("NUM\tAllele1\tAllele2\tSSG\n");

        //debugging specific alleles
        int index1 = -1, index2 = -1;
        if (!debugAlleles.equals("")){
            String s[] = debugAlleles.split(",");
            for (int i = 0; i < numreads; i++){
                if (HLAnames[i].equals(s[0])){
                    index1 = i;
                }
                if (HLAnames[i].equals(s[1])){
                    index2 = i;
                }
                if (index1 > -1 && index2 > -1){
                    out.printf("INFO: debugging %s\t%s\t%s\t%s\n",s[0],s[1],index1,index2);
                    double dl = CalculateLikelihood(index1,index2,HLAreads2,true);
                    break;
                }
            }
        }

        //Pre-process homozygous combinations to determine top possible alleles (for efficiency)
        int numreads2 = HLAnames2.length;
        Alleles2Digit = new Hashtable();
        Loci = new ArrayList<String>();
        double[] AlleleLikelihoods2 = new double[numreads];
        for (int i = 0; i < numreads; i++){
            name1 = HLAnames[i].substring(4);
            String [] n1 = name1.split("\\*");
            numcombinations++;
            AlleleLikelihoods2[i] = CalculateLikelihood(i,i,HLAreads,false);
            if (AlleleLikelihoods2[i] < 0){
                name2 = n1[0] + "*" + n1[1].substring(0, 2);
                if (!Loci.contains(n1[0])){Loci.add(n1[0]);}
                if (!Alleles2Digit.containsKey(name2)){
                    Alleles2Digit.put(name2, AlleleLikelihoods2[i]);
                }else if ((Double) Alleles2Digit.get(name2) < AlleleLikelihoods2[i]){
                    Alleles2Digit.put(name2, AlleleLikelihoods2[i]);
                }
            }
        }

        //Sort alleles at 2 digit resolution for each locus
        AllelesToSearch = new ArrayList<String>();
        for (int i = 0; i < Loci.size(); i++){
            Enumeration k = Alleles2Digit.keys();
            Hashtable AllelesAtLoci = new Hashtable();

            //find alleles at the locus
            while( k.hasMoreElements() ){
                name1 = k.nextElement().toString();
                String [] n1 = name1.split("\\*");
                if (Loci.get(i).equals(n1[0])){
                    AllelesAtLoci.put(-1 * (Double) Alleles2Digit.get(name1), name1);
                }
            }

            //Sort alleles at locus, mark top six 2-digit classes for deep search
            int num = 1;
            Vector v = new Vector(AllelesAtLoci.keySet());
            Collections.sort(v);
            for (Enumeration e = v.elements(); e.hasMoreElements();) {
                Double key = Double.valueOf(e.nextElement().toString());
                String allele = AllelesAtLoci.get(key).toString();
                if (num <= 10){
                    AllelesToSearch.add(allele);
                    
                    num++;
                }
                //out.printf("%s\t%s\n",allele,key);
            }
        }

        //Iterate through allele pairs to calculate likelihoods
        if (true){
            numcombinations = 0;
            for (int i = 0; i < numreads; i++){
                name1 = HLAnames[i].substring(4);
                String [] n1 = name1.split("\\*");
                if (AllelesToSearch.contains(n1[0] + "*" + n1[1].substring(0, 2))){
                    //out.printf("1: %s\n",name1);
                    //frq1 = Double.parseDouble((String) AlleleFrequencies.get(name1).toString());
                    //if (frq1 > minfrq){
                    for (int j = i; j < numreads; j++){
                        name2 = HLAnames[j].substring(4);
                        String [] n2 = name2.split("\\*");
                        if (AllelesToSearch.contains(n2[0] + "*" + n2[1].substring(0, 2))){
                            if ((HLAstartpos[i] < HLAstoppos[j]) && (HLAstartpos[j] < HLAstoppos[i])){
                                numcombinations++;
                                AlleleLikelihoods[i][j] = CalculateLikelihood(i,j,HLAreads,false);
                                if (AlleleLikelihoods[i][j] < 0){
                                    out.printf("%s\t%s\t%s\t%.2f\n",numcombinations,name1,name2,AlleleLikelihoods[i][j]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private double CalculateLikelihood(int a1, int a2, String[] HLAalleles, boolean debug){
        //Calculates likelihood for specific allele pair
        String read1 = HLAalleles[a1];
        String read2 = HLAalleles[a2];
        int start1 = HLAstartpos[a1];
        int start2 = HLAstartpos[a2];
        int stop1 = HLAstoppos[a1];
        int stop2 = HLAstoppos[a2];
        double likelihood = 0;
        int pos, index;
        char c1, c2;
        

        for (int i = 0; i < positions.length; i++){
            pos = positions[i];
            if (pos < stop1 && pos > start1 && pos < stop2 && pos > start2){
                index = GenotypeIndex(read1.charAt(pos-start1),read2.charAt(pos-start2));
                if (index > -1){
                    likelihood = likelihood + baseLikelihoods[i][index];
                    if (debug){
                        c1 = read1.charAt(pos-start1);
                        c2 = read2.charAt(pos-start2);
                        out.printf("INFO: DEBUG %s\t%s\t%s\t%s\t%s\t%s\t%.2f\n",HLAnames[a1],HLAnames[a2],pos,c1,c2,index,likelihood);
                    }
                }
            }
        }
        return likelihood;
    }
}
