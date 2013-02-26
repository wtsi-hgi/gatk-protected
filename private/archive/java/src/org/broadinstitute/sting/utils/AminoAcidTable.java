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

package org.broadinstitute.sting.utils;

import java.util.HashMap;

/**
 * A simple {codon -> amino acid name} lookup table.
 * Handles differences between mitochondrial and nuclear genes.
 */
public class AminoAcidTable {


    protected static final AminoAcid UNKNOWN = new AminoAcid("X" , "Unknown", "Unk");
    protected static final AminoAcid ISOLEUCINE = new AminoAcid("I" , "Isoleucine", "Ile");
    protected static final AminoAcid LEUCINE = new AminoAcid("L" , "Leucine", "Leu");
    protected static final AminoAcid VALINE =    new AminoAcid("V" , "Valine", "Val");
    protected static final AminoAcid PHENYLALANINE =    new AminoAcid("F" , "Phenylalanine", "Phe");
    protected static final AminoAcid METHIONINE = new AminoAcid("M" , "Methionine", "Met");
    protected static final AminoAcid CYSTEINE = new AminoAcid("C" , "Cysteine", "Cys");
    protected static final AminoAcid ALANINE = new AminoAcid("A" , "Alanine", "Ala");
    protected static final AminoAcid STOP_CODON = new AminoAcid("*" , "Stop Codon", "Stop");
    protected static final AminoAcid GLYCINE = new AminoAcid("G" , "Glycine", "Gly");
    protected static final AminoAcid PROLINE = new AminoAcid("P" , "Proline", "Pro");
    protected static final AminoAcid THEONINE = new AminoAcid("T" , "Threonine", "Thr");
    protected static final AminoAcid SERINE = new AminoAcid("S" , "Serine", "Ser");
    protected static final AminoAcid TYROSINE = new AminoAcid("Y" , "Tyrosine", "Tyr");
    protected static final AminoAcid TRYPTOPHAN = new AminoAcid("W" , "Tryptophan", "Trp");
    protected static final AminoAcid GLUTAMINE = new AminoAcid("Q" , "Glutamine", "Gln");
    protected static final AminoAcid ASPARAGINE = new AminoAcid("N" , "Asparagine", "Asn");
    protected static final AminoAcid HISTIDINE = new AminoAcid("H" , "Histidine", "His");
    protected static final AminoAcid GLUTAMIC_ACID = new AminoAcid("E" , "Glutamic acid", "Glu");
    protected static final AminoAcid ASPARTIC_ACID = new AminoAcid("D" , "Aspartic acid", "Asp");
    protected static final AminoAcid LYSINE = new AminoAcid("K" , "Lysine", "Lys");
    protected static final AminoAcid ARGININE = new AminoAcid("R" , "Arginine", "Arg");

    protected static HashMap<String, AminoAcid> aminoAcidTable = new HashMap<String, AminoAcid>();
    protected static HashMap<String, AminoAcid> mitochondrialAminoAcidTable = new HashMap<String, AminoAcid>();

    static {
        //populate the tables
        aminoAcidTable.put("ATT", ISOLEUCINE);
        aminoAcidTable.put("ATC", ISOLEUCINE);
        aminoAcidTable.put("ATA", ISOLEUCINE);


        aminoAcidTable.put("CTT", LEUCINE);
        aminoAcidTable.put("CTC", LEUCINE);
        aminoAcidTable.put("CTA", LEUCINE);
        aminoAcidTable.put("CTG", LEUCINE);
        aminoAcidTable.put("TTA", LEUCINE);
        aminoAcidTable.put("TTG", LEUCINE);


        aminoAcidTable.put("GTT", VALINE);
        aminoAcidTable.put("GTC", VALINE);
        aminoAcidTable.put("GTA", VALINE);
        aminoAcidTable.put("GTG", VALINE);


        aminoAcidTable.put("TTT", PHENYLALANINE);
        aminoAcidTable.put("TTC", PHENYLALANINE);


        aminoAcidTable.put("ATG", METHIONINE);

        aminoAcidTable.put("TGT", CYSTEINE);
        aminoAcidTable.put("TGC", CYSTEINE);

        aminoAcidTable.put("GCT", ALANINE);
        aminoAcidTable.put("GCC", ALANINE);
        aminoAcidTable.put("GCA", ALANINE);
        aminoAcidTable.put("GCG", ALANINE);


        aminoAcidTable.put("GGT", GLYCINE);
        aminoAcidTable.put("GGC", GLYCINE);
        aminoAcidTable.put("GGA", GLYCINE);
        aminoAcidTable.put("GGG", GLYCINE);


        aminoAcidTable.put("CCT", PROLINE);
        aminoAcidTable.put("CCC", PROLINE);
        aminoAcidTable.put("CCA", PROLINE);
        aminoAcidTable.put("CCG", PROLINE);




        aminoAcidTable.put("ACT", THEONINE);
        aminoAcidTable.put("ACC", THEONINE);
        aminoAcidTable.put("ACA", THEONINE);
        aminoAcidTable.put("ACG", THEONINE);



        aminoAcidTable.put("TCT", SERINE);
        aminoAcidTable.put("TCC", SERINE);
        aminoAcidTable.put("TCA", SERINE);
        aminoAcidTable.put("TCG", SERINE);
        aminoAcidTable.put("AGT", SERINE);
        aminoAcidTable.put("AGC", SERINE);

        aminoAcidTable.put("TAT", TYROSINE);
        aminoAcidTable.put("TAC", TYROSINE);



        aminoAcidTable.put("TGG", TRYPTOPHAN);


        aminoAcidTable.put("CAA", GLUTAMINE);
        aminoAcidTable.put("CAG", GLUTAMINE);


        aminoAcidTable.put("AAT", ASPARAGINE);
        aminoAcidTable.put("AAC", ASPARAGINE);


        aminoAcidTable.put("CAT", HISTIDINE);
        aminoAcidTable.put("CAC", HISTIDINE);


        aminoAcidTable.put("GAA", GLUTAMIC_ACID);
        aminoAcidTable.put("GAG", GLUTAMIC_ACID);



        aminoAcidTable.put("GAT", ASPARTIC_ACID);
        aminoAcidTable.put("GAC", ASPARTIC_ACID);


        aminoAcidTable.put("AAA", LYSINE);
        aminoAcidTable.put("AAG", LYSINE);


        aminoAcidTable.put("CGT", ARGININE);
        aminoAcidTable.put("CGC", ARGININE);
        aminoAcidTable.put("CGA", ARGININE);
        aminoAcidTable.put("CGG", ARGININE);
        aminoAcidTable.put("AGA", ARGININE);
        aminoAcidTable.put("AGG", ARGININE);


        aminoAcidTable.put("TAA", STOP_CODON );
        aminoAcidTable.put("TAG", STOP_CODON);
        aminoAcidTable.put("TGA", STOP_CODON);


        //populate the mitochondrial AA table
        mitochondrialAminoAcidTable.putAll(aminoAcidTable);
        mitochondrialAminoAcidTable.put("AGA", STOP_CODON);
        mitochondrialAminoAcidTable.put("AGG", STOP_CODON);
        mitochondrialAminoAcidTable.put("ATA", METHIONINE);
        mitochondrialAminoAcidTable.put("TGA", TRYPTOPHAN);
    }

    /**
     * Returns the amino acid encoded by the given codon in a eukaryotic genome.
     *
     * @param codon The 3-letter mRNA nucleotide codon 5' to 3'. Expects T's instead of U's. Not case sensitive.
     *
     * @return The amino acid matching the given codon, or the UNKNOWN amino acid if the codon string doesn't match anything
     */
    public static AminoAcid getEukaryoticAA(String codon) {
        codon = codon.toUpperCase();
        final AminoAcid aa = aminoAcidTable.get(codon);
        return aa == null ? UNKNOWN : aa;
    }


    /**
     * Returns the amino acid encoded by the given codon in a mitochondrial genome.
     *
     * @param codon The 3-letter mRNA nucleotide codon 5' to 3'. Expects T's instead of U's. Not case sensitive.
     * @param isFirstCodon If this is the 1st codon in the gene, then "ATT" encodes Methyonine
     *
     * @return The amino acid matching the given codon in mitochondrial genes, or the UNKNOWN amino acid if the codon string doesn't match anything
     */
    public static AminoAcid getMitochondrialAA(String codon, boolean isFirstCodon) {
        codon = codon.toUpperCase();
        final AminoAcid aa = mitochondrialAminoAcidTable.get(codon);
        if(aa == null) {
            return UNKNOWN;
        } else if(isFirstCodon && codon.equals("ATT")) {
            return METHIONINE; //special case - 'ATT' in the first codon of a mitochondrial gene codes for methionine instead of isoleucine
        } else {
            return aa;
        }
    }
}
