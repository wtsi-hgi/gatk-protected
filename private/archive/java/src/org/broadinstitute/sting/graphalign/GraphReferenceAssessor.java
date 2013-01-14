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

package org.broadinstitute.sting.gatk.walkers.graphalign;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.util.*;
import java.io.*;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.variant.utils.BaseUtils;

/**
 * A completely experimental read walker that consumes a graphical reference emitted by GraphReferenceBuilder as a
 * serialized java object and evaluates the number of mismatches to both the flat reference and the graphical
 * reference for each read [Not for public use and will change drastically in the future].
 */
public class GraphReferenceAssessor extends ReadWalker<Integer, Integer> {
    @Output
    PrintStream out;

    @Argument(fullName="graphFile", shortName="GF", doc="", required=true)
    String graphFile = null;
    ObjectInputStream graphSerialStream = null;

    @Argument(fullName="MAX", shortName="MAX", doc="", required=false)
    int MAXREADS = -1;

    @Argument(fullName="ignore0MM", shortName="I0", doc="", required=false)
    boolean IGNORE_0_MM = false;

    @Argument(fullName="DEBUG", shortName="DB", doc="", required=false)
    int DEBUG_LEVEL = 0;

    static boolean DEBUG = false;
    static boolean DEBUG2 = false; // crazy level

    @Argument(fullName="read", doc="", required=false)
    String onlyDoRead = null;

    ReferenceGraph graphRef = null;

    public void initialize() {
        super.initialize();

        DEBUG = DEBUG_LEVEL > 0;
        DEBUG2 = DEBUG_LEVEL > 1; // crazy level
        
        try {
            logger.info("Reading graph reference " + graphFile );
            graphSerialStream = new ObjectInputStream( new FileInputStream( graphFile ) );
            graphRef = (ReferenceGraph)graphSerialStream.readObject();
            graphRef.setDebugPrinting(DEBUG);
            graphRef.validateGraph();            
            logger.info(graphRef.toBriefString());
        } catch ( FileNotFoundException e ) {
            throw new StingException("Couldn't open file " + graphFile, e);
        } catch ( IOException e ) {
            throw new StingException("Couldn't write to file " + graphFile, e);
        } catch ( ClassNotFoundException e ) {
            throw new StingException("Couldn't read ReferenceGraph from file " + graphFile, e);
        }
    }

    private static MismatchCounter countMismatches(byte[] ref, int refOffset, byte[] bases, byte[] quals, int basesOffset, int length) {
        MismatchCounter mm = new MismatchCounter();

        for ( int i = 0; i < length; i++ ) {
            byte rawRefBase = ref[i + refOffset];
            byte rawReadBase = bases[i + basesOffset];
            int fragBase = BaseUtils.simpleBaseToBaseIndex((char) rawRefBase);
            int readBase = BaseUtils.simpleBaseToBaseIndex((char)rawReadBase);

            boolean mmP = fragBase != -1 && readBase != -1 && fragBase != readBase;
            if ( mmP ) {
                mm.nMM++;
                mm.qSum += quals != null ? quals[i + basesOffset] : 0;
            }

            if ( GraphReferenceAssessor.DEBUG2 )
                System.out.printf("%s%d %c %c %s %b%n", Utils.dupString(' ', basesOffset + 2), basesOffset, (char)rawRefBase, (char)rawReadBase, mm, mmP);
        }

        return mm;
    }

    private static MismatchCounter countMismatches(byte[] ref, byte[] bases, byte[] quals) {
        return countMismatches(ref, 0, bases, quals, 0, bases.length);
    }

    private static MismatchCounter countMismatchesOnGraph( ReferenceGraph graph, Collection<Fragment> frags, int fragOffset, byte[] bases, byte[] quals, int readOffset ) {
        if ( frags.size() == 0 )
            throw new RuntimeException("Fragment list is empty!");

        MismatchCounter minNMM = MismatchCounter.MAX_VALUE;

        for ( Fragment next : frags ) {
            MismatchCounter recNMM = countMismatchesOnGraph( graph, next, 0, bases, quals, readOffset );
            minNMM = minNMM.min( recNMM );
        }

        return minNMM;
    }

    private static MismatchCounter countMismatchesOnGraph( ReferenceGraph graph, Fragment frag, int fragOffset, byte[] bases, byte[] quals, int readOffset ) {
        if ( GraphReferenceAssessor.DEBUG )System.out.printf("%sfrag %s -> %d%n", Utils.dupString(' ', readOffset + 2), frag, readOffset);

        MismatchCounter mm = new MismatchCounter();

        if ( readOffset < bases.length ) {
            int nRemainingBases = bases.length - readOffset;
            int cmpLength = frag.getBaseLengthFrom(fragOffset, nRemainingBases);       // how many bases over in the fragment are we from the offset
            MismatchCounter fragMM = countMismatches(frag.getUnderlyingBases(), frag.getUnderlyingOffset() + fragOffset, bases, quals, readOffset, cmpLength);
            mm.add(fragMM);

//            // still have some counting to do
//            for ( int i = 0; i < baseLength; i++ ) {
//                int fragBaseOffset = fragOffset + i;
//                int readBaseOffset = readOffset + i;
//
//                byte rawFragBase = frag.getBase(fragBaseOffset);
//                byte rawReadBase = bases[readBaseOffset];
//                int fragBase = BaseUtils.simpleBaseToBaseIndex((char)rawFragBase);
//                int readBase = BaseUtils.simpleBaseToBaseIndex((char)rawReadBase);
//
//                boolean mmP = fragBase != -1 && readBase != -1 && fragBase != readBase;
//                if ( mmP ) nMM++;

            if ( nRemainingBases > cmpLength ) {
                MismatchCounter recMM = countMismatchesOnGraph( graph, graph.outgoingFragments(frag), 0, bases, quals, readOffset + cmpLength );
                mm.add(recMM);
            }
        }

        if ( GraphReferenceAssessor.DEBUG ) System.out.printf("%s=> %s%n", Utils.dupString(' ', readOffset + 2), mm);
        return mm;
    }

    private static MismatchCounter countMismatchesOnGraph(ReferenceGraph graph, SAMRecord read) {
        if ( GraphReferenceAssessor.DEBUG ) System.out.printf("countMismatchesOnGraph( read=%s%n", read.getReadName());
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(read);
        MismatchCounter minNMM = MismatchCounter.MAX_VALUE;

        for ( Fragment frag : graph.getStartingFragment(loc) ) {
            int fragOffset = frag.getFragOffsetFrom(loc);       // how many bases over in the fragment are we from the offset

            if ( GraphReferenceAssessor.DEBUG )
                System.out.printf("  countMismatchesOnGraph frag=%s loc=%s bases=%s offset=%d%n", frag, loc, read.getReadString(), fragOffset);

            MismatchCounter recNMM = countMismatchesOnGraph(graph, frag, fragOffset, read.getReadBases(), read.getBaseQualities(), 0);
            minNMM = minNMM.min( recNMM );
        }

        return minNMM;
    }

    public Integer map(ReferenceContext refArg, SAMRecord read, ReadMetaDataTracker metaDataTracker) {

        if ( MAXREADS-- == 0 ) {
            System.exit(0);
        } else if ( onlyDoRead != null && ! read.getReadName().equals(onlyDoRead) ) {
            ;
        } else if ( ! read.getReadUnmappedFlag() && read.getCigar().numCigarElements() == 1 ) {
            try {
                byte[] ref = refArg.getBases();
                // we're all XM
                int nMMFromRead = (Short)read.getAttribute("NM");
                MismatchCounter nAlignedMM = countMismatches(ref, read.getReadBases(), read.getBaseQualities());
                if ( ! IGNORE_0_MM || nAlignedMM.nMM > 0 ) {
                    MismatchCounter nGraphMM = countMismatchesOnGraph(graphRef, read);
                    MismatchCounter deltaMM = nAlignedMM.minus(nGraphMM);

                    out.printf("%50s with %5s at %10s: mismatches: %3d (delta %3d) -- %3d %3d -- %3d %3d -- delta %3d %3d%n",
                            read.getReadName(), read.getCigarString(), GenomeLocParser.createGenomeLoc(read),
                            nMMFromRead, nMMFromRead - nAlignedMM.nMM,
                            nAlignedMM.nMM, nAlignedMM.qSum,
                            nGraphMM.nMM, nGraphMM.qSum,
                            deltaMM.nMM, deltaMM.qSum);

                    if ( deltaMM.nMM < 0 || deltaMM.qSum < 0 )
                        throw new StingException(read.getReadName() + " is miscalculated");
                }
            } catch ( Exception e ) {
                System.out.printf("Exception at %s at %s%n", read.getReadName(), GenomeLocParser.createGenomeLoc(read));
                throw new RuntimeException(e);
            }
        } else {
                ; // don't do anything
        }

        return 0;
    }

    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     *
     * @return
     */
    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer readScore, Integer data) {
        return data + readScore;
    }

    public void onTraversalDone(Integer data) {
        out.printf(data.toString());
    }
}

class MismatchCounter {
    int nMM = 0;
    int qSum = 0;

    public static MismatchCounter MAX_VALUE = new MismatchCounter(Integer.MAX_VALUE, Integer.MAX_VALUE );

    public MismatchCounter() {}

    public MismatchCounter(int nMM, int qSum) {
        this.nMM = nMM;
        this.qSum = qSum;
    }

    public void add(MismatchCounter that) {
        this.nMM += that.nMM;
        this.qSum += that.qSum;
    }

    public MismatchCounter min(MismatchCounter that) {
        int cmpQSum = Integer.valueOf(this.qSum).compareTo(that.qSum);
        if ( cmpQSum < 0 ) { return this; }
        else if ( cmpQSum > 0 ) { return that; }
        else if ( this.nMM < that.nMM ) { return this; }
        else if ( this.nMM > that.nMM ) { return that; }
        else { return this; }
    }

    public MismatchCounter minus(MismatchCounter that) {
        return new MismatchCounter(this.nMM - that.nMM, this.qSum - that.qSum);
    }

    public String toString() { return String.format("[MM %d %d]", nMM, qSum); }
}