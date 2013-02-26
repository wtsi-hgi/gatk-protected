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
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;

import java.io.*;
import java.util.ArrayList;
import java.util.Hashtable;
/**
 * ImputeAllelesWalker fills in missing intronic info for HLA alleles based on the the most similar HLA allele per read
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class ImputeAllelesWalker extends ReadWalker<Integer, Integer> {
    @Output
    PrintStream out;           

    String HLAdatabaseFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA_DICTIONARY.sam";
//    String ClosestAllelesFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.CLASS1.closest";
    String ClosestAllelesFile = "/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.CLASS2.closest";

    boolean DatabaseLoaded = false;
    boolean DEBUG = false;

    ArrayList<String> HLAreads = new ArrayList<String>();
    ArrayList<String> HLAcigars = new ArrayList<String>();
    ArrayList<String> HLAnames = new ArrayList<String>();
    ArrayList<String> HLApositions = new ArrayList<String>();
    double[] SingleAlleleFrequencies;

    int numHLAlleles = 0;
    int[] HLAstartpos;
    int[] HLAstoppos;
    int minstartpos = 0;
    int maxstoppos = 0;

    int HLA_A_start = 30018310;
    int HLA_A_end = 30021211;
    int HLA_B_start = 31430239;
    int HLA_B_end = 31432914;
    int HLA_C_start = 31344925;
    int HLA_C_end = 31347827;
    int HLA_DQA1_start = 32713161;
    int HLA_DQA1_end = 32719407;
    int HLA_DQB1_start = 32735635;
    int HLA_DQB1_end = 32742444;
    int HLA_DPA1_start = 33140772;
    int HLA_DPA1_end = 33149356;
    int HLA_DPB1_start = 33151738;
    int HLA_DPB1_end = 33162954;
    int HLA_DRB1_start = 32654525;
    int HLA_DRB1_end = 32665540;


    ArrayList<String> PolymorphicSites = new ArrayList<String>();

    Hashtable ClosestAllele = new Hashtable();
    int iAstart = -1, iAstop = -1, iBstart = -1, iBstop = -1, iCstart = -1, iCstop = -1, iDRBstart = -1, iDRBstop = -1, iDQAstart = -1, iDQAstop = -1, iDQBstart = -1, iDQBstop = -1, iDPAstart = -1, iDPAstop = -1, iDPBstart = -1, iDPBstop = -1;
    CigarParser formatter = new CigarParser();
    
    public Integer reduceInit() {
    if (!DatabaseLoaded){
            try{
                out.printf("Reading HLA database ...\n");
                FileInputStream fstream = new FileInputStream(HLAdatabaseFile);
                DataInputStream in = new DataInputStream(fstream);
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine; String [] s = null;
                //Read File Line By Line
                int i = 0;
                while ((strLine = br.readLine()) != null)   {
                    s = strLine.split("\\t");

                    if (s.length>=10){
                        //Parse the reads with cigar parser
                        HLAreads.add(formatter.FormatRead(s[5],s[9]));
                        HLAcigars.add(s[5]);
                        HLAnames.add(s[0]);

                        HLApositions.add(s[3]);
                        if (s[0].indexOf("HLA_A") > -1){
                            if (iAstart < 0){iAstart=i;}
                            iAstop = i; i++;
                        }else if (s[0].indexOf("HLA_B") > -1){
                            if (iBstart < 0){iBstart=i;}
                            iBstop = i; i++;
                        }else if (s[0].indexOf("HLA_C") > -1){
                            if (iCstart < 0){iCstart=i;}
                            iCstop = i; i++;
                        }else if (s[0].indexOf("HLA_DRB1") > -1){
                            if (iDRBstart < 0){iDRBstart=i;}
                            iDRBstop = i; i++;
                        }else if (s[0].indexOf("HLA_DQA1") > -1){
                            if (iDQAstart < 0){iDQAstart=i;}
                            iDQAstop = i; i++;
                        }else if (s[0].indexOf("HLA_DQB1") > -1){
                            if (iDQBstart < 0){iDQBstart=i;}
                            iDQBstop = i; i++;
                        }else if (s[0].indexOf("HLA_DPA1") > -1){
                            if (iDPAstart < 0){iDPAstart=i;}
                            iDPAstop = i; i++;
                        }else if (s[0].indexOf("HLA_DPB1") > -1){
                            if (iDPBstart < 0){iDPBstart=i;}
                            iDPBstop = i; i++;
                        }
                    }
                }
                in.close();
                int n = HLApositions.size(); numHLAlleles = n;
                HLAstartpos = new int[n]; HLAstoppos = new int[n];
                SingleAlleleFrequencies = new double[n];


                for (i = 0; i < n; i++){
                    //Find start and stop positions for each allele
                    HLAstartpos[i]=Integer.parseInt(HLApositions.get(i));
                    HLAstoppos[i]=HLAstartpos[i]+HLAreads.get(i).length()-1;
                    if (minstartpos == 0){minstartpos = HLAstartpos[i];}
                    minstartpos = Math.min(minstartpos, HLAstartpos[i]);
                    maxstoppos = Math.max(maxstoppos, HLAstoppos[i]);
                    SingleAlleleFrequencies[i]=0.0;
                    //Initialize matrix of probabilities / likelihoods

                }
                out.printf("DONE! Read %s alleles\n",HLAreads.size());
            }catch (Exception e){//Catch exception if any
              System.err.println("ImputeAllelsWalker Error: " + e.getMessage());
            }

            try{
                out.printf("Reading closest allele file ...");
                FileInputStream fstream = new FileInputStream(ClosestAllelesFile);
                DataInputStream in = new DataInputStream(fstream);
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine; String [] s = null;
                //Read File Line By Line
                int count = 0;
                while ((strLine = br.readLine()) != null)   {
                    s = strLine.split("\\t");
                    ClosestAllele.put(s[0], s[2]);
//                    out.printf("loading: %s\t%s\n",s[0],s[2]);
                    count++;
                }
                in.close();
                out.printf("Done! Read %s alleles\n",count);
            }catch (Exception e){//Catch exception if any
              System.err.println("ImputeAllelsWalker Error: " + e.getMessage());
            }

            char c;
            DatabaseLoaded = true;
            
            out.printf("Imputing alleles ...\n");

            if (DEBUG){
                //out.printf("Astart[%s]\tAstop[%s]\tBstart[%s]\tBstop[%s]\tCstart[%s]\tCstop[%s]\tnumAlleles[%s]\n",iAstart,iAstop,iBstart,iBstop,iCstart,iCstop,numHLAlleles);
            }
        }
        return 0;
    }


    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        int readstart = read.getAlignmentStart();
        int readstop = read.getAlignmentEnd();
        int startimputation = 0, stopimputation = 0;

        String s1 = formatter.FormatRead(read.getCigarString(), read.getReadString());
        char c;
        String readstring = "", name = "", cigar = "", qualitystring = "";
        int numM = 0, numI = 0, numD = 0;

        name = read.getReadName();
        
        String matchedAllele = (String) ClosestAllele.get(name);

        //out.printf("%s\t%s\n",name,matchedAllele);
        int index = HLAnames.indexOf(matchedAllele);
        
        String matchedRead = HLAreads.get(index);

        if (name.indexOf("HLA_A") > -1){
            startimputation = HLA_A_start;
            stopimputation = HLA_A_end;
        } else if (name.indexOf("HLA_B") > -1){
            startimputation = HLA_B_start;
            stopimputation = HLA_B_end;
        } else if (name.indexOf("HLA_C") > -1){
            startimputation = HLA_C_start;
            stopimputation = HLA_C_end;
        } else if (name.indexOf("HLA_DRB1") > -1){
            startimputation = HLA_DRB1_start;
            stopimputation = HLA_DRB1_end;
        } else if (name.indexOf("HLA_DQA1") > -1){
            startimputation = HLA_DQA1_start;
            stopimputation = HLA_DQA1_end;
        } else if (name.indexOf("HLA_DQB1") > -1){
            startimputation = HLA_DQB1_start;
            stopimputation = HLA_DQB1_end;
        } else if (name.indexOf("HLA_DPA1") > -1){
            startimputation = HLA_DPA1_start;
            stopimputation = HLA_DPA1_end;
        } else if (name.indexOf("HLA_DPB1") > -1){
            startimputation = HLA_DPB1_start;
            stopimputation = HLA_DPB1_end;
        }

        //out.printf("DEBUG %s\t%s\t%s\t%s\t%s\n",name,matchedAllele,index,startimputation,stopimputation);
        for (int i = startimputation; i <= stopimputation; i++){
            //if position is within read
            if (i >= readstart && i <= readstop){
                c = s1.charAt(i-readstart);
                //if position is not missing
                if (c != 'D'){
                    readstring = readstring + c;
                    qualitystring = qualitystring + 'I';
                    numM++;
                    if (numD > 0){
                        cigar = cigar + String.valueOf(numD) + "D";
                        numD = 0;
                    } else if (numI > 0){
                        cigar = cigar + String.valueOf(numI) + "I";
                        numI = 0;
                    }
                //if position is missing, get base from matched allele
                }else{
                    c = matchedRead.charAt(i-HLAstartpos[index]);
                    //if matched allele is also missing / deleted at position
                    if (c == 'D'){
                        numD++;
                        if (numM > 0){
                            cigar = cigar + String.valueOf(numM) + "M";
                            numM = 0;
                        }
                    //if matched allele is not missing / deleted at position
                    }else{
                        readstring = readstring + c;
                        qualitystring = qualitystring + 'I';
                        numM++;
                        if (numD > 0){
                            cigar = cigar + String.valueOf(numD) + "D";
                            numD = 0;
                        } else if (numI > 0){
                            cigar = cigar + String.valueOf(numI) + "I";
                            numI = 0;
                        }
                    }
                }
            //if position is outside of range of read, look at matched allele
            }else{
                //if within range of matched allele
                if (i >= HLAstartpos[index] && i <= HLAstoppos[index]){
                    c = matchedRead.charAt(i-HLAstartpos[index]);
                    //if matched allele is also missing / deleted at position
                    if (c == 'D'){
                        numD++;
                        if (numM > 0){
                            cigar = cigar + String.valueOf(numM) + "M";
                            numM = 0;
                        }
                    //if matched allele is not missing / deleted at position
                    }else{
                        readstring = readstring + c;
                        qualitystring = qualitystring + 'I';
                        numM++;
                        if (numD > 0){
                            cigar = cigar + String.valueOf(numD) + "D";
                            numD = 0;
                        } else if (numI > 0){
                            cigar = cigar + String.valueOf(numI) + "I";
                            numI = 0;
                        }
                    }
                }else{
                    numD++;
                    if (numM > 0){
                        cigar = cigar + String.valueOf(numM) + "M";
                        numM = 0;
                    }
                }
            }
        }

        if (numM > 0){
            cigar = cigar + String.valueOf(numM) + "M";
        }else if(numD > 0){
            cigar = cigar + String.valueOf(numD) + "D";
        }else if(numI > 0){
            cigar = cigar + String.valueOf(numI) + "I";
        }
        
        out.printf("%s\t0\t6\t%s\t99\t%s\t*\t0\t0\t%s\t%s\n",name,startimputation,cigar,readstring,qualitystring);
        
        
        return 1;
    }




    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}

