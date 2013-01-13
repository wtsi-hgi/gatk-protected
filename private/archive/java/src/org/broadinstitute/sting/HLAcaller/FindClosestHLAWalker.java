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
import java.util.Hashtable;

/**
 * Finds the most similar HLA allele for each read (helps detect misalignments). Usage: java -jar GenomeAnalysisTK.jar -T FindClosestHLA -I INPUT.bam -R /broad/1KG/reference/human_b36_both.fasta -L INPUT.interval | grep -v INFO | sort -k1 > OUTPUT
 * @author shermanjia
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class FindClosestHLAWalker extends ReadWalker<Integer, Integer> {
    @Output
    protected PrintStream out;

    @Argument(fullName = "debugRead", shortName = "debugRead", doc = "Print match score for read", required = false)
    public String debugRead = "";

    @Argument(fullName = "findFirst", shortName = "findFirst", doc = "For each read, stop when first HLA allele is found with concordance = 1", required = false)
    public boolean findFirst = false;

    @Argument(fullName = "DEBUG", shortName = "DEBUG", doc = "Debug walker", required = false)
    public boolean DEBUG = false;
    
    @Argument(fullName = "debugAllele", shortName = "debugAllele", doc = "Print match score for allele", required = false)
    public String debugAllele = "";
    
    @Argument(fullName = "useInterval", shortName = "useInterval", doc = "Use only these intervals", required = false)
    public String intervalFile = "";

    @Argument(fullName = "dictionary", shortName = "dictionary", doc = "bam file of HLA ditionary", required = false)
    public String HLAdictionaryFile ="/humgen/gsa-scr1/GSA/sjia/454_HLA/HLA/HLA.nuc.sam";

    @Argument(fullName = "onlyfrequent", shortName = "onlyfrequent", doc = "Only consider alleles with frequency > 0.0001", required = false)
    public boolean ONLYFREQUENT = false;
    
    @Argument(fullName = "HLAdictionary", shortName = "HLAdictionary", doc = "HLA dictionary file", required = true)
    public String HLAdatabaseFile = "HLA_DICTIONARY.txt";

    @Argument(fullName = "PolymorphicSites", shortName = "PolymorphicSites", doc = "file containing polymorphic sites within the HLA", required = true)
    public String PolymorphicSitesFile = "HLA_POLYMORPHIC_SITES.txt";
    
    HLAFileReader HLADictionaryReader = new HLAFileReader();    

    boolean DatabaseLoaded = false;
    ArrayList<String> ClosestAlleles = new ArrayList<String>();

    String[] HLAnames, HLAreads;
    Integer[] HLAstartpos, HLAstoppos, PolymorphicSites,NonPolymorphicSites;
    double[] SingleAlleleFrequencies;

    double[] nummatched, concordance, numcompared;
    int numHLAlleles = 0;
    int minstartpos = 0;
    int maxstoppos = 0;
    int numpolymorphicsites = 0, numnonpolymorphicsites = 0, pos =0;

    Hashtable AlleleFrequencies = new Hashtable();
    int iAstart = -1, iAstop = -1, iBstart = -1, iBstop = -1, iCstart = -1, iCstop = -1;
    CigarParser formatter = new CigarParser();
    int [][] intervals; int numIntervals;

    public Integer reduceInit() { 
        if (!DatabaseLoaded){
            DatabaseLoaded = true;

            //Load HLA dictionary
            out.printf("INFO  Loading HLA dictionary ... ");

            HLADictionaryReader.ReadFile(HLAdatabaseFile);
            HLAreads = HLADictionaryReader.GetSequences();
            HLAnames = HLADictionaryReader.GetNames();
            HLAstartpos = HLADictionaryReader.GetStartPositions();
            HLAstoppos = HLADictionaryReader.GetStopPositions();
            minstartpos = HLADictionaryReader.GetMinStartPos();
            maxstoppos = HLADictionaryReader.GetMaxStopPos();
            
            out.printf("Done! %s HLA alleles loaded.\n",HLAreads.length);

            nummatched = new double[HLAreads.length];
            concordance = new double[HLAreads.length];
            numcompared = new double[HLAreads.length];

            //Load list of polymorphic sites
            PolymorphicSitesFileReader siteFileReader = new PolymorphicSitesFileReader();
            siteFileReader.ReadFile(PolymorphicSitesFile);
            PolymorphicSites = siteFileReader.GetPolymorphicSites();
            NonPolymorphicSites = siteFileReader.GetNonPolymorphicSites();
            numpolymorphicsites = PolymorphicSites.length;
            numnonpolymorphicsites = NonPolymorphicSites.length;

            if (!intervalFile.equals("")){
                TextFileReader fileReader = new TextFileReader();
                fileReader.ReadFile(intervalFile);
                String[] lines = fileReader.GetLines();
                intervals = new int[lines.length][2];
                for (int i = 0; i < lines.length; i++) {
                    String[] s = lines[i].split(":");
                    String[] intervalPieces = s[1].split("-");
                    intervals[i][0] = Integer.valueOf(intervalPieces[0]);
                    intervals[i][1] = Integer.valueOf(intervalPieces[1]);
                }
                numIntervals = intervals.length;
            }
            
            out.printf("INFO  %s polymorphic and %s non-polymorphic sites found in HLA dictionary\n",numpolymorphicsites,numnonpolymorphicsites);
            out.printf("INFO  Comparing reads to database ...\n");

            if (DEBUG){
                //out.printf("Astart[%s]\tAstop[%s]\tBstart[%s]\tBstop[%s]\tCstart[%s]\tCstop[%s]\tnumAlleles[%s]\n",iAstart,iAstop,iBstart,iBstop,iCstart,iCstop,numHLAlleles);
            }
        }
        return 0;
    }

    private double CalculateConcordance(SAMRecord read){
        int readstart = read.getAlignmentStart();
        int readstop = read.getAlignmentEnd();
        char c1, c2;
        double maxConcordance = 0.0, freq = 0.0, minFreq = 0.0;
        String s1 = formatter.FormatRead(read.getCigarString(), read.getReadString());
        String s2;
        int allelestart, allelestop;

        if (ONLYFREQUENT){
            minFreq = 0.0001;
        }

        for (int i = 0; i < HLAreads.length; i++){
            nummatched[i] = 0; concordance[i] = 0; numcompared[i] = 0;
            freq = GetAlleleFrequency(HLAnames[i]);
            //Get concordance between read and specific allele
            if (readstart <= HLAstoppos[i] && readstop >= HLAstartpos[i] && freq > minFreq){
                s2 = HLAreads[i];
                
                allelestart = HLAstartpos[i];
                allelestop = HLAstoppos[i];

                //Polymorphic sites: always increment denominator, increment numerator when bases are concordant
                for (int j = 0; j < numpolymorphicsites; j++){
                    pos = PolymorphicSites[j];
                    if (DEBUG == true){
                        out.printf("DEBUG\tPOS\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(),HLAnames[i],pos,allelestart,allelestop,IsWithin(pos,readstart,readstop), IsWithin(pos,allelestart,allelestop),IsWithinInterval(pos));
                    }
                    if (pos >= readstart && pos <= readstop && pos >= allelestart && pos <= allelestop && IsWithinInterval(pos)){
                        c1 = s1.charAt(pos-readstart);
                        c2 = s2.charAt(pos-allelestart);
                        if (c1 != 'D' && c2 != 'D'){//allow for deletions (sequencing errors)
                            numcompared[i]++;
                            if (c1 == c2){
                                nummatched[i]++;
                            }else{
                                if (debugRead.equals(read.getReadName()) && debugAllele.equals(HLAnames[i])){
                                    out.printf("DEBUG\t%s\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(), HLAnames[i], j, pos,c1,c2);
                                }
                            }
                        }
                    }
                }

                //Non-polymorphic sites: increment denominator only when bases are discordant
                if (numcompared[i] > 0){
                    for (int j = 0; j < numnonpolymorphicsites; j++){
                        pos = NonPolymorphicSites[j];
                        if (pos >= readstart && pos <= readstop && pos >= allelestart && pos <= allelestop && IsWithinInterval(pos)){
                            c1 = s1.charAt(pos-readstart);
                            c2 = s2.charAt(pos-allelestart);
                            if (c1 != c2 && c1 != 'D' && c2 != 'D'){//allow for deletions (sequencing errors)
                                numcompared[i]++;
                                if (debugRead.equals(read.getReadName()) && debugAllele.equals(HLAnames[i])){
                                    out.printf("DEBUG\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(), HLAnames[i], j, c1,c2);
                                }
                            }
                        }
                    }
                }
            
                //Update concordance array
                concordance[i]=nummatched[i]/numcompared[i];
                if (concordance[i] > maxConcordance){maxConcordance = concordance[i];}
                if (DEBUG == true){
                    out.printf("DEBUG\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(),HLAnames[i],concordance[i],numcompared[i],numcompared[i]-nummatched[i]);
                }
                if (debugRead.equals(read.getReadName()) && debugAllele.equals(HLAnames[i])){
                    out.printf("DEBUG\t%s\t%s\t%s\t%s\t%s\n",read.getReadName(),HLAnames[i],concordance[i],numcompared[i],numcompared[i]-nummatched[i]);
                }
                if (findFirst && (concordance[i] == 1)){
                    break;
                }
            }
        
        }

        return maxConcordance;
    }

    private double FindMaxAlleleFrequency(double maxConcordance){
        //finds the max frequency of the alleles that share the maximum concordance with the read of interest
        double freq, maxFreq = 0.0;
        for (int i = 0; i < HLAreads.length; i++){
            if (concordance[i] == maxConcordance && maxConcordance > 0){
                freq = GetAlleleFrequency(HLAnames[i]);
                if (freq > maxFreq){maxFreq = freq;}
            }
        }
        return maxFreq;
    }

    private boolean IsWithin(int pos, int start, int stop){
        return pos >= start && pos <= stop;
    }

    private boolean IsWithinInterval(int pos){
        boolean isWithinInterval = false;
        for (int i = 0; i < numIntervals; i++){
            if (pos >= intervals[i][0] && pos <= intervals[i][1]){
                isWithinInterval = true;
                break;
            }
        }
        return isWithinInterval;
    }
    
    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        //Calculate concordance for this read and all overlapping reads
        if (DEBUG == true){
            out.printf("%s\t%s\n",read.getReadName(),read.getMappingQuality());
        }
        if (read.getMappingQuality() > 0 || DEBUG == true){
            double maxConcordance = CalculateConcordance(read);
            String stats = "", topAlleles = "";
            if (maxConcordance > 0 || DEBUG == true){
                String readname = read.getReadName(), allelename = ""; double freq;
                //For input bam files that contain HLA alleles, find and print allele frequency
                out.printf("%s\t%s-%s", readname,read.getAlignmentStart(),read.getAlignmentEnd());

                //Print concordance statistics between this read and the most similar HLA allele(s)

                for (int i = 0; i < HLAreads.length; i++){
                    if (concordance[i] == maxConcordance){
                        //freq = GetAlleleFrequency(HLAnames[i]);
                        if (topAlleles.equals("")){
                            topAlleles = HLAnames[i];
                        }else{
                            topAlleles = topAlleles + "," + HLAnames[i];
                        }
                        stats = String.format("%.1f\t%.3f\t%.0f\t%.0f",1.0,concordance[i],numcompared[i],numcompared[i]-nummatched[i]);
                        
                    }
                }
                out.printf("\t%s\t%s\t%s\n",stats,topAlleles,maxConcordance);
            }
        }
        return 1;
    }

    private double GetAlleleFrequency(String allelename){
        double frequency = 0.0;
        //Truncate names to 4-digit "A*0101" format
        if (allelename.length() >= 10){
            allelename=allelename.substring(4,10);
        }else{
            allelename=allelename.substring(4);
        }
        if (AlleleFrequencies.containsKey(allelename)){
            frequency = Double.parseDouble((String) AlleleFrequencies.get(allelename).toString());
        }else{
            frequency=0.0001;
        }
        return frequency;
    }


    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    @Override
    public void onTraversalDone(Integer result) {
        // Double check traversal result to make count is the same.
        // TODO: Is this check necessary?
        out.println("[REDUCE RESULT] Traversal result is: " + result);
    }    
}

