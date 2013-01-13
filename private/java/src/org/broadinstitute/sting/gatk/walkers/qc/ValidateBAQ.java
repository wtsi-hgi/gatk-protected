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

package org.broadinstitute.sting.gatk.walkers.qc;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;

/**
 * Walks over the input data set, calculating the number of reads seen for diagnostic purposes.
 * Can also count the number of reads matching a given criterion using read filters (see the
 * --read-filter command line argument).  Simplest example of a read-backed analysis.
 */
@BAQMode(QualityMode = BAQ.QualityMode.DONT_MODIFY, ApplicationTime = ReadTransformer.ApplicationTime.HANDLED_IN_WALKER)
@Reference(window=@Window(start=-5,stop=5))
@Requires({DataSource.READS, DataSource.REFERENCE})
public class ValidateBAQ extends ReadWalker<Integer, Integer> {
    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    @Argument(doc="maximum read length to apply the BAQ calculation too",required=false)
    protected int maxReadLen = 1000;

    @Argument(doc="",required=false)
    protected int bw = 7;

    @Argument(doc="",required=false)
    protected boolean samtoolsMode = false;

    @Argument(doc="only operates on reads with this name",required=false)
    protected String readName = null;

    @Argument(doc="If true, all differences are errors", required=false)
    protected boolean strict = false;

    @Argument(doc="prints info for each read", required=false)
    protected boolean printEachRead = false;

    @Argument(doc="Also prints out detailed comparison information when for known calculation differences", required=false)
    protected boolean alsoPrintWarnings = false;

    @Argument(doc="Include reads without BAQ tag", required=false)
    protected boolean includeReadsWithoutBAQTag = false;

    @Argument(doc="x each read is processed", required=false)
    protected int magnification = 1;

    @Argument(doc="Profile performance", required=false)
    protected boolean profile = false;

    int counter = 0;

    BAQ baqHMM = null;         // matches current samtools parameters

    public void initialize() {
        if ( samtoolsMode )
            baqHMM = new BAQ(1e-3, 0.1, bw, (byte)0, true);
        else
            baqHMM = new BAQ();
    }

    long goodReads = 0, badReads = 0;

    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {

        if ( (readName == null || readName.equals(read.getReadName())) && read.getReadLength() <= maxReadLen && (includeReadsWithoutBAQTag || BAQ.hasBAQTag(read) ) ) {
            if ( baqHMM.excludeReadFromBAQ(read) )
                return 0;

            if ( profile ) {
                profileBAQ(ref, read);
            } else {
                validateBAQ(ref, read);
            }

            return 1;
        }

        return 0;
    }

    SimpleTimer tagTimer = new SimpleTimer("from.tag");
    SimpleTimer baqReadTimer = new SimpleTimer("baq.read");
    SimpleTimer glocalTimer = new SimpleTimer("hmm.glocal");

    private void profileBAQ(ReferenceContext ref, SAMRecord read) {
        IndexedFastaSequenceFile refReader = this.getToolkit().getReferenceDataSource().getReference();
        BAQ.BAQCalculationResult baq = null;

        tagTimer.restart();
        for ( int i = 0; i < magnification; i++ ) { BAQ.calcBAQFromTag(read, false, includeReadsWithoutBAQTag); }
        tagTimer.stop();

        baqReadTimer.restart();
        for ( int i = 0; i < magnification; i++ ) { baqHMM.baqRead(read, refReader, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.DONT_MODIFY ); }
        baqReadTimer.stop();

        glocalTimer.restart();
        for ( int i = 0; i < magnification; i++ )
            baqHMM.baqRead(read, refReader, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.DONT_MODIFY);
        glocalTimer.stop();
    }


    private void validateBAQ(ReferenceContext ref, SAMRecord read) {
        IndexedFastaSequenceFile refReader = this.getToolkit().getReferenceDataSource().getReference();
        byte[] baqFromTag = BAQ.calcBAQFromTag(read, false, includeReadsWithoutBAQTag);
        if (counter++ % 1000 == 0 || printEachRead) out.printf("Checking read %s (%d)%n", read.getReadName(), counter);

        BAQ.BAQCalculationResult baq = baqHMM.calcBAQFromHMM(read, refReader);

        boolean fail = false;
        boolean print = false;
        int badi = 0;

        if ( BAQ.hasBAQTag(read) ) {
            for ( badi = 0; badi < baqFromTag.length; badi++ ) {
                if ( baqFromTag[badi] != baq.bq[badi] ) {
                    if ( cigarLength(read) != read.getReadLength() ) {
                        print = true;
                        fail = false;
                        out.printf("  different, but cigar length != read length%n");
                        break;
                    }
                    if (MathUtils.arrayMin(read.getBaseQualities()) == 0) {
                        print = true;
                        fail = strict;
                        out.printf("  different, but Q0 base detected%n");
                        break;
                    }
                    else if (readHasSoftClip(read) && ! samtoolsMode) {
                        print = true;
                        fail = strict;
                        out.printf("  different, but soft clip detected%n");
                        break;
                    } else if (readHasDeletion(read) ) { // && ! samtoolsMode) {
                        print = true;
                        fail = strict;
                        out.printf("  different, but deletion detected%n");
                        break;
                    } else if ( baq.bq[badi] < baqHMM.getMinBaseQual() ) {
                        print = fail = true;
                        out.printf("  Base quality %d < min %d", baq.bq[badi], baqHMM.getMinBaseQual());
                        break;
                    } else {
                        print = fail = true;
                        break;
                    }
                }
            }
            if ( fail || print )
                badReads++;
            else
                goodReads++;
        }

        if ( fail || printEachRead || ( print && alsoPrintWarnings ) ) {
            byte[] pos = new byte[baq.bq.length];
            for ( int i = 0; i < pos.length; i++ ) pos[i] = (byte)i;

            out.printf("  read length   : %d%n", read.getReadLength());
            out.printf("  read start    : %d (%d unclipped)%n", read.getAlignmentStart(), read.getUnclippedStart());
            out.printf("  cigar         : %s%n", read.getCigarString());
            out.printf("  ref bases     : %s%n", new String(baq.refBases));
            out.printf("  read bases    : %s%n", new String(read.getReadBases()));
            out.printf("  ref length    : %d%n", baq.refBases.length);
            out.printf("  BQ tag        : %s%n", read.getStringAttribute(BAQ.BAQ_TAG));
            if ( BAQ.hasBAQTag(read) ) printQuals("  BQ deltas     : ", getBAQDeltas(read), true);
            printQuals("  original quals: ", read.getBaseQualities(), true);
            printQuals("  baq      quals: ", baq.bq, true);
            printQuals("  positions     : ", pos, true);
            printQuals("  original quals: ", read.getBaseQualities());
            if ( BAQ.hasBAQTag(read) ) printQuals("  tag      quals: ", baqFromTag);
            printQuals("  hmm      quals: ", baq.bq);
            out.printf("  read bases    : %s%n", new String(read.getReadBases()));
            out.println(Utils.dupString('-', 80));
        }


        if ( fail )
            throw new StingException(String.format("BAQ from read and from HMM differ in read %s at position %d: tag qual = %d, hmm qual = %d",
                    read.getReadName(), badi, baqFromTag[badi], baq.bq[badi]));
    }

    private final static boolean readHasSoftClip(SAMRecord read) {
        for (CigarElement e : read.getCigar().getCigarElements()) {
            if ( e.getOperator() == CigarOperator.SOFT_CLIP )
                return true;
        }

        return false;
    }

    private final static boolean readHasDeletion(SAMRecord read) {
        for (CigarElement e : read.getCigar().getCigarElements()) {
            if ( e.getOperator() == CigarOperator.DELETION )
                return true;
        }

        return false;
    }

    public final void printQuals( String prefix, byte[] quals ) {
        printQuals(prefix, quals, false);
    }

    public final void printQuals( String prefix, byte[] quals, boolean asInt ) {
        printQuals(out, prefix, quals, asInt);
    }

    public final static void printQuals( PrintStream out, String prefix, byte[] quals, boolean asInt ) {
        out.print(prefix);
        for ( int i = 0; i < quals.length; i++) {
            if ( asInt ) {
                out.printf("%2d", (int)quals[i]);
                if ( i+1 != quals.length ) out.print(",");
            } else
                out.print((char)(quals[i]+33));
        }
        out.println();
    }

    /**
     * Get the BAQ delta bytes from the tag in read.  Returns null if no BAQ tag is present.
     * @param read
     * @return
     */
    public static byte[] getBAQDeltas(SAMRecord read) {
        byte[] baq = BAQ.getBAQTag(read);
        if ( baq != null ) {
            byte[] deltas = new byte[baq.length];
            for ( int i = 0; i < deltas.length; i++)
                deltas[i] = (byte)(-1 * (baq[i] - 64));
            return deltas;
        } else
            return null;
    }

    private int cigarLength(SAMRecord read) {
        int readI = 0;
        for ( CigarElement elt : read.getCigar().getCigarElements() ) {
            int l = elt.getLength();
            switch (elt.getOperator()) {
                case N: // cannot handle these
                    return 0;
                case H : case P : // ignore pads and hard clips
                    break;
                case S :
                case I :
                    readI += l;
                    break;
                case D : break;
                case M :
                    readI += l;
                    break;
                default:
                    throw new ReviewedStingException("BUG: Unexpected CIGAR element " + elt + " in read " + read.getReadName());
            }
        }
        return readI;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer nreads) {
        if ( profile ) {
            out.printf("n.reads baq.per.read calculation time.in.secs%n");
            printTimer(nreads, tagTimer);
            printTimer(nreads, glocalTimer);
            printTimer(nreads, baqReadTimer);
        } else {
            out.printf("total reads BAQ'd %d; concordant BAQ reads %d %.4f; discordant BAQ reads %d %.4f%n", nreads,
                    goodReads, (100.0 * goodReads) / nreads,
                    badReads, (100.0 * badReads) / nreads);
        }
    }

    private final void printTimer(int nreads, SimpleTimer timer) {
        out.printf("%d %d %s %.2f%n", nreads, magnification, timer.getName(), timer.getElapsedTime());
    }
}

