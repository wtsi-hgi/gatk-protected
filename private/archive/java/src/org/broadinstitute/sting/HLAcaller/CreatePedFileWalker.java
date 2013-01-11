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
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

/**
 * Creates a ped file of SNPs and amino acids coded as SNPs given an input ped file with 4-digit HLA alleles. Usage: java -jar GenomeAnalysisTK.jar -T CreatePedFile --allelesFile INPUT.ped -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-sc\
r1/GSA/sjia/454_HLA/HLA/HLA.combined.4digitUnique.bam > OUTPUT.log
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CreatePedFileWalker extends ReadWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @Argument(fullName = "allelesFile", shortName = "allelesFile", doc = "Create ped file for HLA alleles named in this file", required = true)
    public String alleleNamesFile = "";

    @Argument(fullName = "pedIntervals", shortName = "pedIntervals", doc = "Create genotypes in these intervals", required = false)
    public String pedIntervalsFile = "";

    @Argument(fullName = "HLAexonIntervals", shortName = "HLAexonIntervals", doc = "HLA exonic intervals", required = false)
    public String exonIntervalsFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_EXON_POSITIONS.txt";

    @Argument(fullName = "DNAcode", shortName = "DNAcode", doc = "Amino acid codes", required = false)
    public String dnaCodesFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/DNA_CODE.txt";

    @Argument(fullName = "PrintDNA", shortName = "PrintDNA", doc = "Print DNA sequences", required = false)
    public boolean PrintDNA = false;

    @Argument(fullName = "PrintAA", shortName = "PrintAA", doc = "Print Amino Acid sequences", required = false)
    public boolean PrintAA = true;
    
    String[] HLAnames, HLAreads, inputFileContents;
    Integer[] HLAstartpos, HLAstoppos;
    ArrayList<String> HLAnamesAL, HLAreadsAL;
    ArrayList<Integer> HLAstartposAL, HLAstopposAL;
    int[][] intervals; String[][] exonIntervals;
    int numIntervals;

    ReadCigarFormatter formatter = new ReadCigarFormatter();
    char c;
    boolean DEBUG = false;
    boolean FilesLoaded = false;
    int HLA_A_start    = 30018310, HLA_A_end    = 30021211;
    int HLA_C_start    = 31344925, HLA_C_end    = 31347827;
    int HLA_B_start    = 31430239, HLA_B_end    = 31432914;
    int HLA_DRB1_start = 32654846, HLA_DRB1_end = 32665497;
    int HLA_DQA1_start = 32713214, HLA_DQA1_end = 32718519;
    int HLA_DQB1_start = 32735991, HLA_DQB1_end = 32742362;
    int HLA_DPA1_start = 33144405, HLA_DPA1_end = 33149325;
    int HLA_DPB1_start = 33151797, HLA_DPB1_end = 33161993;


    String[] SNPnames;
    String SNPname;
    int start, end;
    Integer I;

    Hashtable indexer = new Hashtable();
    Hashtable DNAcode = new Hashtable();

    public Integer reduceInit() {
        if (!FilesLoaded){
            FilesLoaded = true;
            HLAnamesAL = new ArrayList<String>();
            HLAreadsAL = new ArrayList<String>();
            HLAstartposAL = new ArrayList<Integer>();
            HLAstopposAL = new ArrayList<Integer>();

            TextFileReader fileReader = new TextFileReader();
            fileReader.ReadFile(alleleNamesFile);
            inputFileContents = fileReader.GetLines();


            //Determine intervals
            if (!pedIntervalsFile.equals("")){
                fileReader = new TextFileReader();
                fileReader.ReadFile(pedIntervalsFile);
                String[] lines = fileReader.GetLines();
                intervals = new int[lines.length][2];
                for (int i = 0; i < lines.length; i++) {
                    String[] s = lines[i].split(":");
                    String[] intervalPieces = s[0].split("-");
                    intervals[i][0] = Integer.valueOf(intervalPieces[0]);
                    intervals[i][1] = Integer.valueOf(intervalPieces[1]);
                }
                numIntervals = intervals.length;
                for (int i = 0; i < numIntervals; i++){
                    out.printf("INFO  Interval %s: %s-%s\n",i+1,intervals[i][0],intervals[i][1]);
                }
            }

            //load HLA exonic intervals
            if (!exonIntervalsFile.equals("")){
                fileReader = new TextFileReader();
                fileReader.ReadFile(exonIntervalsFile);
                String[] lines = fileReader.GetLines();
                exonIntervals = new String[lines.length][5];
                for (int i = 0; i < lines.length; i++) {
                    String[] s = lines[i].split("\t");
                    String[] intervalPieces = s[1].split("-");
                    exonIntervals[i][1] = intervalPieces[0];
                    exonIntervals[i][2] = intervalPieces[1];
                    exonIntervals[i][0] = s[0]; // Locus
                    exonIntervals[i][3] = s[2]; // Exon number
                    exonIntervals[i][4] = s[3]; // +/- strand
                }
                numIntervals = exonIntervals.length;
                for (int i = 0; i < numIntervals; i++){
                    out.printf("INFO  HLA-%s %s (%s): %s-%s\n",exonIntervals[i][0],exonIntervals[i][3],exonIntervals[i][4],exonIntervals[i][1],exonIntervals[i][2]);
                }
            }

            //load amino-acid coding DNA triplets
            if (!dnaCodesFile.equals("")){
                fileReader = new TextFileReader();
                fileReader.ReadFile(dnaCodesFile);
                String[] lines = fileReader.GetLines();
                for (int i = 0; i < lines.length; i++) {
                    String[] s = lines[i].split("\t");
                    DNAcode.put(s[0],s[1]);
                }

                Enumeration e = DNAcode.keys();
                while( e.hasMoreElements() ){
                    String key = e.nextElement().toString();
                    out.printf("INFO %s encodes %s\n",key,DNAcode.get(key));
                }
            }
        }
        return 0;
    }

    private String[][] GetExonIntervals(String locus, boolean isForwardStrand){
        int numExons = 0; int exonNum;
        for (int i = 0; i < exonIntervals.length; i++){
            if (exonIntervals[i][0].equals(locus)){
                numExons++;
            }
        }
        String[][] ExonIntervals = new String[numExons][5];
        if (isForwardStrand){exonNum = 1;}else{exonNum = ExonIntervals.length;}
        for (int i = 0; i < exonIntervals.length; i++){
            if (exonIntervals[i][0].equals(locus)){
                ExonIntervals[exonNum-1]=exonIntervals[i];
                if (isForwardStrand){
                    exonNum++;
                }else{
                    exonNum--;
                }
            }
        }
        return ExonIntervals;
    }

    private int BaseCharToInt(char c){
        switch(c){
            case 'A': return 1;
            case 'C': return 2;
            case 'G': return 3;
            case 'T': return 4;
            default: return -1;
        }
    }

    private char Complement(char c){
        switch(c){
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            default: return '0';
        }
    }

    private char GetAminoAcid(String codon){
        if (DNAcode.containsKey(codon)){
            return DNAcode.get(codon).toString().charAt(0);
        }else{
            return '0';
        }
    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        HLAnamesAL.add(read.getReadName());
        HLAreadsAL.add(formatter.FormatRead(read.getCigarString(), read.getReadString()));
        HLAstartposAL.add(read.getAlignmentStart());
        HLAstopposAL.add(read.getAlignmentEnd());
        return 1;
    }

    private String PrintGenotypes(String ID, String alleleName1, String alleleName2, int startpos, int stoppos){

        String error = "";
        //prints genotypes for allele1 and allele2 at given interval
        int i1 = GetAlleleIndex(alleleName1);
        int i2 = GetAlleleIndex(alleleName2);
        String s1, s2;
        int start1, start2, stop1, stop2;
        char c1, c2;

        if (i1 > -1){
            s1 = HLAreads[i1];
            start1 = HLAstartpos[i1];
            stop1 = HLAstoppos[i1];
        }else{
            error = error + "INFO  " + alleleName1 + " for " + ID + " not found in HLA dictionary\n";
            s1 = "";
            start1 = -1;
            stop1 = -1;
        }

        if (i2 > -1){
            s2 = HLAreads[i2];
            start2 = HLAstartpos[i2];
            stop2 = HLAstoppos[i2];
        }else{
            error = error + "INFO  " + alleleName2 + " for " + ID + " not found in HLA dictionary\n";
            s2 = "";
            start2 = -1;
            stop2 = -1;
        }
        
        for (int pos = startpos; pos <= stoppos; pos++){
            c1 = GetBase(pos,s1,start1,stop1);
            c2 = GetBase(pos,s2,start2,stop2);
            out.printf("\t%s %s",c1,c2);
        }
        return error;
    }

private String PrintAminoAcids(String ID, String alleleName1, String alleleName2, String[][] ExonIntervals){

        String error = "";
        //prints genotypes for allele1 and allele2 at given interval
        int i1 = GetAlleleIndex(alleleName1);
        int i2 = GetAlleleIndex(alleleName2);
        String s1, s2;
        int start1, start2, stop1, stop2;
        char c1, c2;
        boolean isForwardStrand = false;
        if (ExonIntervals[0][4].equals("+")){isForwardStrand=true;}

        int AAcount=0;
        int baseCount=0;
        String codon1 = ""; String codon2 = "";

        if (i1 > -1){
            s1 = HLAreads[i1];
            start1 = HLAstartpos[i1];
            stop1 = HLAstoppos[i1];
        }else{
            s1 = "";
            start1 = -1;
            stop1 = -1;
            error = error + "INFO  " + alleleName1 + " for " + ID + " not found in HLA dictionary\n";
        }

        if (i2 > -1){
            s2 = HLAreads[i2];
            start2 = HLAstartpos[i2];
            stop2 = HLAstoppos[i2];
        }else{
            s2 = "";
            start2 = -1;
            stop2 = -1;
            error = error + "INFO  " + alleleName2 + " for " + ID + " not found in HLA dictionary\n";
        }

        int i;
        for (int exonNum = 1; exonNum <= ExonIntervals.length; exonNum++){
            if (isForwardStrand){i=exonNum-1;}else{i=ExonIntervals.length-exonNum;}
            int exonStart = Integer.parseInt(ExonIntervals[i][1]);
            int exonStop = Integer.parseInt(ExonIntervals[i][2]);
            for (int pos = exonStart; pos <= exonStop; pos++){
                c1 = GetBase(pos,s1,start1,stop1);
                c2 = GetBase(pos,s2,start2,stop2);
                if (!isForwardStrand){
                    c1 = Complement(c1);
                    c2 = Complement(c2);
                }
                if (baseCount < 3){
                    if (isForwardStrand){
                        codon1 = codon1 + c1;
                        codon2 = codon2 + c2;
                    }else{
                        codon1 = c1 + codon1;
                        codon2 = c2 + codon2;
                    }
                    baseCount++;
                }

                if (baseCount == 3){
                    out.printf("\t%s %s",GetAminoAcid(codon1),GetAminoAcid(codon2));
                    baseCount = 0;
                    AAcount++;
                    codon1 = "";
                    codon2 = "";
                }
            }
        }
        if (baseCount > 0){
            //Print stop or start codon depending on strandedness
            if (isForwardStrand){out.printf("\tO O");}else{out.printf("\tM M");}
        }
        
        return error;
    }

    private char GetBase(int pos, String str, int start, int stop){
        char base;
        if (pos >= start && pos <= stop){
            base = str.charAt(pos-start);
            if (base == 'D'){base = '0';}
        }else{
            base = '0';
        }
        return base;
    }

    private int GetAlleleIndex(String alleleName){
        //Find first allele that matches name, or matches part of name for 2-digit allele
        int i;
        for (i = 0; i < HLAnames.length; i++){
            if (HLAnames[i].indexOf(alleleName) > -1){
                return i;
            }
        }
        return -1;
        
    }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    private String GetAlleleName(String locus, String sep, String allele){
        if (allele.length() > 1){
            return locus + sep + allele;
        }else{
            return locus + sep + "0000";
        }
    }

    public void onTraversalDone(Integer numreads) {
        HLAnames = HLAnamesAL.toArray(new String[numreads]);
        HLAreads = HLAreadsAL.toArray(new String[numreads]);
        HLAstartpos = HLAstartposAL.toArray(new Integer[numreads]);
        HLAstoppos = HLAstopposAL.toArray(new Integer[numreads]);
        String star = "*";
        String error = "";

        //out.printf("INFO %s alleles in dictionary\n",HLAnames.length);
        String[][] A_exons = GetExonIntervals("A",true);
        String[][] B_exons = GetExonIntervals("B",false);
        String[][] C_exons = GetExonIntervals("C",false);
        String[][] DRB1_exons = GetExonIntervals("DRB1",false);
        String[][] DQB1_exons = GetExonIntervals("DQB1",false);
        String[][] DQA1_exons = GetExonIntervals("DQA1",true);
        String[][] DPB1_exons = GetExonIntervals("DPB1",true);
        String[][] DPA1_exons = GetExonIntervals("DPA1",false);
        //Print individual info and genotypes
        for (int i = 0; i < inputFileContents.length; i++){
            String[] s = inputFileContents[i].split(" ");
            //out.printf("%s\t%s\n",inputFileContents[i],s.length);
            if (s.length > 10){
                error = "";
                out.printf("%s\t%s\t%s\t%s\t%s\t%s",s[0],s[1],s[2],s[3],s[4],s[5]);
                String HLA_A_1 = GetAlleleName("HLA_A",star,s[6]);
                String HLA_A_2 = GetAlleleName("HLA_A",star,s[7]);
                String HLA_B_1 = GetAlleleName("HLA_B",star,s[8]);
                String HLA_B_2 = GetAlleleName("HLA_B",star,s[9]);
                String HLA_C_1 = GetAlleleName("HLA_C",star,s[10]);
                String HLA_C_2 = GetAlleleName("HLA_C",star,s[11]);
                String HLA_DPA1_1 = GetAlleleName("HLA_DPA1",star,s[12]);
                String HLA_DPA1_2 = GetAlleleName("HLA_DPA1",star,s[13]);
                String HLA_DPB1_1 = GetAlleleName("HLA_DPB1",star,s[14]);
                String HLA_DPB1_2 = GetAlleleName("HLA_DPB1",star,s[15]);
                String HLA_DQA1_1 = GetAlleleName("HLA_DQA1",star,s[16]);
                String HLA_DQA1_2 = GetAlleleName("HLA_DQA1",star,s[17]);
                String HLA_DQB1_1 = GetAlleleName("HLA_DQB1",star,s[18]);
                String HLA_DQB1_2 = GetAlleleName("HLA_DQB1",star,s[19]);
                String HLA_DRB1_1 = GetAlleleName("HLA_DRB1",star,s[20]);
                String HLA_DRB1_2 = GetAlleleName("HLA_DRB1",star,s[21]);

                

                if (true) {
                    if (PrintDNA){
                        error = error + PrintGenotypes(s[1], HLA_A_1,HLA_A_2, HLA_A_start,HLA_A_end);
                        error = error + PrintGenotypes(s[1], HLA_C_1,HLA_C_2, HLA_C_start,HLA_C_end);
                        error = error + PrintGenotypes(s[1], HLA_B_1,HLA_B_2, HLA_B_start,HLA_B_end);
                        error = error + PrintGenotypes(s[1], HLA_DRB1_1,HLA_DRB1_2, HLA_DRB1_start,HLA_DRB1_end);
                        error = error + PrintGenotypes(s[1], HLA_DQA1_1,HLA_DQA1_2, HLA_DQA1_start,HLA_DQA1_end);
                        error = error + PrintGenotypes(s[1], HLA_DQB1_1,HLA_DQB1_2, HLA_DQB1_start,HLA_DQB1_end);
                        error = error + PrintGenotypes(s[1], HLA_DPA1_1,HLA_DPA1_2, HLA_DPA1_start,HLA_DPA1_end);
                        error = error + PrintGenotypes(s[1], HLA_DPB1_1,HLA_DPB1_2, HLA_DPB1_start,HLA_DPB1_end);
                    }
                    if (PrintAA){
                        error = error + PrintAminoAcids(s[1], HLA_A_1,HLA_A_2, A_exons);
                        error = error + PrintAminoAcids(s[1], HLA_C_1,HLA_C_2, C_exons);
                        error = error + PrintAminoAcids(s[1], HLA_B_1,HLA_B_2, B_exons);
                        error = error + PrintAminoAcids(s[1], HLA_DRB1_1,HLA_DRB1_2, DRB1_exons);
                        error = error + PrintAminoAcids(s[1], HLA_DQA1_1,HLA_DQA1_2, DQA1_exons);
                        error = error + PrintAminoAcids(s[1], HLA_DQB1_1,HLA_DQB1_2, DQB1_exons);
                        error = error + PrintAminoAcids(s[1], HLA_DPA1_1,HLA_DPA1_2, DPA1_exons);
                        error = error + PrintAminoAcids(s[1], HLA_DPB1_1,HLA_DPB1_2, DPB1_exons);
                    }
                    out.printf("\n");
                    out.printf("%s",error);
                }
            }
        }

        //Prints SNP names for each site
        if (true){
            if (PrintDNA){
                PrintSNPS(HLA_A_start,HLA_A_end);
                PrintSNPS(HLA_C_start,HLA_C_end);
                PrintSNPS(HLA_B_start,HLA_B_end);
                PrintSNPS(HLA_DRB1_start,HLA_DRB1_end);
                PrintSNPS(HLA_DQA1_start,HLA_DQA1_end);
                PrintSNPS(HLA_DQB1_start,HLA_DQB1_end);
                PrintSNPS(HLA_DPA1_start,HLA_DPA1_end);
                PrintSNPS(HLA_DPB1_start,HLA_DPB1_end);
            }

            if (PrintAA){
                PrintAminoAcidSites(A_exons,"A",true);
                PrintAminoAcidSites(C_exons,"C",false);
                PrintAminoAcidSites(B_exons,"B",false);
                PrintAminoAcidSites(DRB1_exons,"DRB1",false);
                PrintAminoAcidSites(DQA1_exons,"DQA1",true);
                PrintAminoAcidSites(DQB1_exons,"DQB1",false);
                PrintAminoAcidSites(DPA1_exons,"DPA1",false);
                PrintAminoAcidSites(DPB1_exons,"DPB1",true);
            }
        }

    }

    private void PrintSNPS(int startpos, int stoppos){
        for (int pos = startpos; pos <= stoppos; pos++){
            SNPname = "CHR6_POS" + String.valueOf(pos);
            out.printf("6\t%s\t0\t%s\n",SNPname,pos);
        }
    }

    private void PrintAminoAcidSites(String[][] ExonIntervals, String locus, boolean isForwardStrand){
        int AAcount=1; int baseCount = 1; int exonNum;

        if (!isForwardStrand){
            for (int i = 1; i <= ExonIntervals.length; i++){
                int exonStart = Integer.parseInt(ExonIntervals[i-1][1]);
                int exonStop = Integer.parseInt(ExonIntervals[i-1][2]);
                for (int pos = exonStart; pos <= exonStop; pos++){
                    if (baseCount == 3){
                        AAcount++;
                        baseCount = 1;
                    }else{
                        baseCount++;
                    }
                }
            }
        }

        for (int i = 1; i <= ExonIntervals.length; i++){
            if (isForwardStrand){exonNum = i;}else{exonNum = ExonIntervals.length - i + 1;}
            int exonStart = Integer.parseInt(ExonIntervals[exonNum-1][1]);
            int exonStop = Integer.parseInt(ExonIntervals[exonNum-1][2]);
            for (int pos = exonStart; pos <= exonStop; pos++){
                if (baseCount == 2){
                    SNPname = locus + "_AA" + String.valueOf(AAcount) + "_E" + exonNum + "_" + String.valueOf(pos);
                    out.printf("6\t%s\t0\t%s\n",SNPname,pos);
                }
                if (baseCount == 3){
                    if (isForwardStrand){AAcount++;}else{AAcount--;}
                    baseCount = 1;
                }else{
                    baseCount++;
                }
            }
        }
    }
}
