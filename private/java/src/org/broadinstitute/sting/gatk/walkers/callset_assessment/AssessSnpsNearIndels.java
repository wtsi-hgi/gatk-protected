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

package org.broadinstitute.sting.gatk.walkers.callset_assessment;

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.ArrayList;

/**
 * Assesses distance of SNP calls to nearest indel call in Exomes.
 * Use --snps and --indels
 */
public class AssessSnpsNearIndels extends RodWalker<Integer, Integer> {

    @Input(fullName="snps", shortName = "snps", doc="Input VCF snps file", required=true)
    public RodBinding<VariantContext> snpTrack;

    @Input(fullName="indels", shortName = "indels", doc="Input VCF indels file", required=true)
    public RodBinding<VariantContext> indelTrack;

    private class TiTv {
        int TiCount = 0, TvCount = 0;

        public TiTv() {};
    }

    private static final int Bin0to5 = 0;
    private static final int Bin6to10 = 1;
    private static final int Bin11to15 = 2;
    private static final int Bin16to20 = 3;
    private static final int Bin21to25 = 4;
    private static final int Bin26to30 = 5;
    private static final int BinMoreThan30 = 6;

    private GenomeLoc previousIndel = null;
    private ArrayList<VariantContext> snpQueue = new ArrayList<VariantContext>();
    private TiTv[] counts = new TiTv[7];
    private GenomeLocParser GLparser = null;

    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out = null;

    public void initialize() {
        GLparser = getToolkit().getGenomeLocParser();

        for (int i = 0; i < 7; i++)
            counts[i] = new TiTv();
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        VariantContext snp = tracker.getFirstValue(snpTrack, context.getLocation());
        VariantContext indel = tracker.getFirstValue(indelTrack, context.getLocation());

        // first add the snp if available
        if ( snp != null && !snp.isFiltered() ) {

            // flush the queue on a new contig
            if ( !snpQueue.isEmpty() && !snpQueue.get(0).getChr().equals(snp.getChr()) ) {
                for ( VariantContext vc : snpQueue )
                    calculateDistance(vc, previousIndel, null);
                snpQueue.clear();
            }

            snpQueue.add(snp);
        }

        // then look for the indel
        if ( indel != null && !indel.isFiltered() ) {

            GenomeLoc loc = GLparser.createGenomeLoc(indel.getChr(), indel.getStart(), indel.getEnd());
            boolean sameContig = !snpQueue.isEmpty() && snpQueue.get(0).getChr().equals(indel.getChr());

            // flush the queue
            for ( VariantContext vc : snpQueue )
                calculateDistance(vc, previousIndel, sameContig ? loc : null);
            snpQueue.clear();

            previousIndel = loc;
        }

        return 1;
    }

    private void calculateDistance(VariantContext snp, GenomeLoc previousIndel, GenomeLoc nextIndel) {

        GenomeLoc loc = GLparser.createGenomeLoc(snp.getChr(), snp.getStart(), snp.getEnd());

        int previousDistance = -1, nextDistance = -1;

        if ( previousIndel != null ) {
            // watch out for spanning deletions
            if ( previousIndel.getStop() > snp.getStart() )
                previousDistance = 0;
            else
                previousDistance = snp.getStart() - previousIndel.getStop();
        }

        if ( nextIndel != null )
            nextDistance = nextIndel.getStart() - loc.getStart();

        if ( previousDistance == -1 && nextDistance == -1 )
            return;

        int distance;
        if ( previousDistance == -1 )
            distance = nextDistance;
        else if ( nextDistance == -1 )
            distance = previousDistance;
        else
            distance = Math.min(previousDistance, nextDistance);

        TiTv obj;
        if ( distance < 0 )
             throw new IllegalStateException("Found a negative distance at " + loc);
        else if ( distance < 6 )
            obj = counts[Bin0to5];
        else if ( distance < 11 )
            obj = counts[Bin6to10];
        else if ( distance < 16 )
            obj = counts[Bin11to15];
        else if ( distance < 21 )
            obj = counts[Bin16to20];
        else if ( distance < 26 )
            obj = counts[Bin21to25];
        else if ( distance < 31 )
            obj = counts[Bin26to30];
        else
            obj = counts[BinMoreThan30];

        if ( BaseUtils.isTransition(snp.getReference().getBases()[0], snp.getAlternateAllele(0).getBases()[0]) )
            obj.TiCount++;
        else
            obj.TvCount++;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {
        // flush the queue
        for ( VariantContext vc : snpQueue )
            calculateDistance(vc, previousIndel, null);

        out.println("Bin\tnumTi\tnumTv\tTi/Tv");
        printLine(counts[Bin0to5], "0to5");
        printLine(counts[Bin6to10], "6to10");
        printLine(counts[Bin11to15], "11to15");
        printLine(counts[Bin16to20], "16to20");
        printLine(counts[Bin21to25], "21to25");
        printLine(counts[Bin26to30], "26to30");
        printLine(counts[BinMoreThan30], ">30");
    }

    private void printLine(TiTv obj, String s) {
        out.println(String.format("%s\t%d\t%d\t%.2f", s, obj.TiCount, obj.TvCount, ((double)obj.TiCount/(double)obj.TvCount)));
    }
}
