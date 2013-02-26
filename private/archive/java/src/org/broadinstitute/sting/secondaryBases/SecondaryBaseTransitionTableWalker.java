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

package org.broadinstitute.sting.gatk.walkers.secondaryBases;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.utils.NamedTable;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.commandline.Output;

import java.util.HashMap;
import java.io.PrintStream;

/**
 * Given a secondary base annotated .bam file and a reference, this walker generates a table of secondary base counts
 * for all called loci in the .bam. Each base call made is an instance included in the table. Specifically, the walker
 * maps the following vector to a count of secondary bases:
 * <HomRef/Het/HomVar, Read Group, Called Genotype, AlleleBalance, Reference Base, Primary Base, Previous Read Base, Secondary Base>.
 */

@Reference(window=@Window(start=-1,stop=1))
public class SecondaryBaseTransitionTableWalker extends LocusWalker<Integer, Integer> {
    @Output
    PrintStream out;

    HashMap<String,Long> counts = new HashMap<String,Long>();
    private UnifiedGenotyperEngine ug;
    private NamedTable altTable;

    public void initialize() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.STANDARD_CONFIDENCE_FOR_CALLING = uac.STANDARD_CONFIDENCE_FOR_EMITTING = 50.0;
        uac.ALL_BASES_MODE = true;
        ug = new UnifiedGenotyperEngine(getToolkit(), uac);

        altTable = new NamedTable();
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        char refBase = Character.toUpperCase(ref.getBaseAsChar());
        ReadBackedPileup pileup = context.getBasePileup();
        int[] baseCounts = pileup.getBaseCounts();
        int length = 0;
        for (int i : baseCounts) {length += i;}
        byte[] contextBases = ref.getBases();
        byte prevBase = (byte)Character.toUpperCase(contextBases[0]);
        byte nextBase = (byte)Character.toUpperCase(contextBases[contextBases.length - 1]);

        if (contextBases.length == 3 && refBase != 'N' && pileup.getBases() != null && pileup.getSecondaryBases() != null) {
            VariantCallContext ugResult = ug.runGenotyper(tracker,ref,context);
            if (ugResult != null && ugResult.vc != null) {
                Genotype res = ugResult.vc.getGenotype(0);
                String call = res.getGenotypeString();
                String type;
                String alleleBalance = "N/A";
                if (res.isHomRef()) {
                    type = "homref";
                }
                else if (!res.isHet()) {type = "homvar";}
                else if (call.contains(Character.toString(refBase))) {
                    type = "het";
                    char alt;
                    if (call.charAt(0) == refBase) {alt = call.charAt(1);}
                    else {alt = call.charAt(0);}
                    double refCount = baseCounts[BaseUtils.simpleBaseToBaseIndex(refBase)];
                    double altCount = baseCounts[BaseUtils.simpleBaseToBaseIndex(alt)];
                    alleleBalance = Double.toString(Math.round(100.0*refCount/(refCount + altCount))/100.0);
                }
                else {type = "bad";}
                if (!type.equals("bad")) {
                    for (PileupElement element : pileup) {
                        char primaryBase = Character.toUpperCase((char)element.getBase());
                        char secondaryBase = Character.toUpperCase((char)element.getSecondBase());
                        String RG = element.getRead().getReadGroup().getReadGroupId();
                        if (secondaryBase != 'N' && secondaryBase != '.' && primaryBase != 'N') {
                            String strandRef;
                            String strandPrimary;
                            String strandPrev;
                            String strandSecondary;
                            String strandCall;
                            if (!element.getRead().getReadNegativeStrandFlag()) {
                                strandRef = Character.toString(refBase);
                                strandPrimary = Character.toString(primaryBase);
                                strandPrev = Character.toString((char)prevBase);
                                strandSecondary = Character.toString(secondaryBase);
                                strandCall = call;
                            }
                            else {
                                strandRef = Character.toString(BaseUtils.simpleComplement(refBase));
                                strandPrimary = Character.toString(BaseUtils.simpleComplement(primaryBase));
                                strandPrev = Character.toString(BaseUtils.simpleComplement((char)nextBase));
                                strandSecondary = Character.toString(BaseUtils.simpleComplement(secondaryBase));
                                strandCall = BaseUtils.simpleReverseComplement(call);
                            }
                            if (strandPrev.charAt(0) != 'N') {
                                String key = type+' '+RG+' '+strandCall+' '+alleleBalance+' '+strandRef+' '+strandPrimary+' '+strandPrev+' '+strandSecondary;
                                if (counts.containsKey(key)) {
                                    counts.put(key, counts.get(key) + Long.valueOf(1));
                                }
                                else {
                                    counts.put(key, Long.valueOf(1));
                                }
                            }
                        }
                    }
                }
            }
        }
    return 1;
    }

    public Integer reduceInit() {return 0;}

    public Integer reduce(Integer value, Integer sum) {return sum + value;}

    public void onTraversalDone(Integer result) {
        out.println(">>>");
        out.println("Type ReadGroup CalledGenotype AlleleBalance ReferenceBase PrimaryBase PreviousBase SecondaryBase Count");
        for (String key : counts.keySet()) {
            out.println(key + ' ' + counts.get(key).toString());
        }
        out.println("Processed " + result.toString() + " loci.");
    }
}