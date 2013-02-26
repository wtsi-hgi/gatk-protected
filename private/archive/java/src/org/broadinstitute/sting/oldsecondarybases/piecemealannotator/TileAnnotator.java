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

package org.broadinstitute.sting.piecemealannotator;

import net.sf.samtools.*;
import org.broadinstitute.sting.secondarybase.*;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.containers.BoundedScoringSet;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;

import java.io.File;
import java.util.ArrayList;
import java.util.Date;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class TileAnnotator extends CommandLineProgram {
    public static TileAnnotator instance = null;

    private String currentContig = "";
    private byte[] refbases;

    @Argument(fullName="sam_tile_in", shortName="STI", doc="SAM tile file", required=false) public File SAM_TILE_IN;
    @Argument(fullName="sam_tile_out", shortName="STO", doc="Annotated SAM tile output file") public File SAM_TILE_OUT;
    @Argument(fullName="reference", shortName="R", doc="The fasta reference") public File REFERENCE;
    @Argument(fullName="bustard_dir", shortName="D", doc="The Bustard directory") public File BUSTARD_DIR;
    @Argument(fullName="training_limit", shortName="TL", doc="Number of reads to train from", required=false) public int TRAINING_LIMIT = 10000;
    @Argument(fullName="run_barcode", shortName="RB", doc="Illumina run barcode") public String RUN_BARCODE;
    @Argument(fullName="cycle_ranges", shortName="CR", doc="Cycle ranges for single-end or paired reads (i.e. '0-75,76-151') (0-based, inclusive)") public String CYCLE_RANGES;
    @Argument(fullName="lane", shortName="L", doc="The lane to process (if not specified, this will be read from the 'sam_tile_in' file)", required=false) public Integer lane;
    @Argument(fullName="tile", shortName="T", doc="The tile to process (if not specified, this will be read from the 'sam_tile_in' file)", required=false) public Integer tile;

    public static void main(String[] argv) {
        instance = new TileAnnotator();
        start(instance, argv);
    }

    protected int execute() {
        ArrayList<Pair<Integer, Integer>> cycleRanges = getCycleRanges(CYCLE_RANGES);
        SecondaryBaseAnnotator sba = new SecondaryBaseAnnotator();

        System.out.printf("%s: Loading training set...\n", (new Date()).toString());
        loadTrainingData(sba);

        System.out.printf("%s: Calling bases...\n", (new Date()).toString());
        callBases(sba, cycleRanges);

        System.out.println("Done.");

        return 0;
    }

    private ArrayList<Pair<Integer, Integer>> getCycleRanges(String cycleRangesString) {
        ArrayList< Pair<Integer, Integer> > cycleRanges = new ArrayList< Pair<Integer, Integer> >();

        String[] pieces = cycleRangesString.split(",");

        Pattern p = Pattern.compile("(\\d+)-(\\d+)");

        for (String piece : pieces) {
            Matcher m = p.matcher(piece);

            if (m.find()) {
                Integer cycleStart = new Integer(m.group(1));
                Integer cycleStop = new Integer(m.group(2));

                cycleRanges.add(new Pair<Integer, Integer>(cycleStart, cycleStop));
            }
        }

        if (cycleRanges.size() == 0) {
            throw new StingException("At least one cycle range must be specified.");
        }

        if (cycleRanges.size() > 2) {
            throw new StingException(cycleRanges.size() + " specified, but we're unable to handle more than 2.");
        }

        return cycleRanges;
    }

    private void loadTrainingData(SecondaryBaseAnnotator sba) {
        IlluminaTile tileParser = new IlluminaTile(BUSTARD_DIR, lane, tile);

        for (RawRead rr : tileParser) { sba.addTrainingRead(rr); }

        tileParser.close();

        sba.doneTraining();
    }

    private void callBases(SecondaryBaseAnnotator sba, ArrayList<Pair<Integer, Integer>> cycleRanges) {
        SAMFileHeader sheader = new SAMFileHeader();
        sheader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        SAMFileWriter swriter =  new SAMFileWriterFactory().makeSAMOrBAMWriter(sheader, true, SAM_TILE_OUT);

        IlluminaTile tileParser = new IlluminaTile(BUSTARD_DIR, lane, tile);

        BasecallingStats bstats = new BasecallingStats();

        for (RawRead rr : tileParser) {
            bstats.update(rr, sba.getFourProbRead(rr));

            byte[] sqtag = sba.getSqTagValue(rr);

            SAMRecord sr = constructSAMRecord(rr, sqtag, sheader, false, false);

            swriter.addAlignment(sr);
        }

        bstats.notifyNow();

        tileParser.close();
        swriter.close();
    }

    private SAMRecord constructSAMRecord(RawRead rr, byte[] sqtag, SAMFileHeader sfh, boolean isPaired, boolean isSecondEndOfPair) {
        SAMRecord sr = new SAMRecord(sfh);

        sr.setReadName(String.format("%s:%d:%d:%d:%d#0", RUN_BARCODE, lane, tile, rr.getXCoordinate(), rr.getYCoordinate()));
        sr.setReadUmappedFlag(true);
        sr.setReadString(rr.getSequenceAsString());
        sr.setBaseQualities(rr.getQuals());
        sr.setAttribute("SQ", sqtag);

        sr.setReadPairedFlag(isPaired);
        if (isPaired) {
            sr.setMateUnmappedFlag(true);
            sr.setFirstOfPairFlag(!isSecondEndOfPair);
            sr.setSecondOfPairFlag(isSecondEndOfPair);
        }

        return sr;
    }

    /*
    private SAMRecord constructSAMRecord(RawRead rr, FourProbRead fpr, SAMFileHeader sfh, boolean isPaired, boolean isSecondEndOfPair) {
        SAMRecord sr = new SAMRecord(sfh);

        sr.setReadName(String.format("%s:%d:%d:%d:%d#0", RUN_BARCODE, lane, tile, rr.getXCoordinate(), rr.getYCoordinate()));
        sr.setReadUmappedFlag(true);
        sr.setReadString(rr.getSequenceAsString());
        sr.setBaseQualities(rr.getQuals());
        sr.setAttribute("SQ", fpr.getSQTag(rr));

        sr.setReadPairedFlag(isPaired);
        if (isPaired) {
            sr.setMateUnmappedFlag(true);
            sr.setFirstOfPairFlag(!isSecondEndOfPair);
            sr.setSecondOfPairFlag(isSecondEndOfPair);
        }

        return sr;
    }

    private void callBases(BasecallingReadModel model, ArrayList<Pair<Integer, Integer>> cycleRanges) {
        SAMFileHeader sheader = new SAMFileHeader();
        sheader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        SAMFileWriter swriter =  new SAMFileWriterFactory().makeSAMOrBAMWriter(sheader, true, SAM_TILE_OUT);

        IlluminaTile tileParser = new IlluminaTile(BUSTARD_DIR, lane, tile);

        BasecallingStats bstats = new BasecallingStats();

        for (RawRead rr : tileParser) {
            FourProbRead fpr = model.call(rr);

            for (int rangeIndex = 0; rangeIndex < cycleRanges.size(); rangeIndex++) {
                FourProbRead fprEnd = fpr.getSubset(cycleRanges.get(rangeIndex).getFirst(), cycleRanges.get(rangeIndex).getSecond());
                RawRead rrEnd = rr.getSubset(cycleRanges.get(rangeIndex).getFirst(), cycleRanges.get(rangeIndex).getSecond());

                SAMRecord sr = constructSAMRecord(rrEnd, fprEnd, sheader, cycleRanges.size() > 1, rangeIndex == 1);

                swriter.addAlignment(sr);
            }

            bstats.update(rr, fpr);
            bstats.notifyOnInterval(10000);
        }

        bstats.notifyNow();

        tileParser.close();
        swriter.close();
    }

    private ArrayList<RawRead> loadTrainingData() {
        FastaSequenceFile2 ref = new FastaSequenceFile2(REFERENCE);
        HashMap<String, SAMRecord> srs = loadTileAlignments(ref);
        return loadGoodReads(srs, BUSTARD_DIR);
    }

    private HashMap<String, SAMRecord> loadTileAlignments(FastaSequenceFile2 ref) {
        HashMap<String, SAMRecord> srs = new HashMap<String, SAMRecord>();
        HashSet<String> seenEnds = new HashSet<String>();
        
        int numPerfect = 0;

        if (SAM_TILE_IN != null && SAM_TILE_IN.exists()) {
            SAMFileReader sreader = new SAMFileReader(SAM_TILE_IN);

            for (SAMRecord sr : sreader) {
                Pattern p = Pattern.compile(":(\\d+):(\\d+):(\\d+):(\\d+)#");
                Matcher m = p.matcher(sr.getReadName());

                if (m.find()) {
                    this.lane = Integer.valueOf(m.group(1));
                    this.tile = Integer.valueOf(m.group(2));
                    int x = Integer.valueOf(m.group(3));
                    int y = Integer.valueOf(m.group(4));
                    boolean end = sr.getReadPairedFlag() && sr.getSecondOfPairFlag();

                    String otherKey = String.format("%d:%d:%b", x, y, !end);
                    String currentKey = String.format("%d:%d:%b", x, y, end);

                    seenEnds.add(currentKey);

                    if (isWellAligned(sr, ref)) {
                        if (srs.containsKey(otherKey) || !seenEnds.contains(otherKey)) {
                            srs.put(currentKey, sr);
                        }

                        if (srs.containsKey(currentKey) && srs.containsKey(otherKey)) {
                            numPerfect++;
                            if (numPerfect % (TRAINING_LIMIT < 1000 ? TRAINING_LIMIT : 1000) == 0) {
                                System.out.println("  " + numPerfect + " well-aligned reads");
                            }
                        }
                    } else {
                        if (srs.containsKey(otherKey)) {
                            srs.remove(otherKey);
                        }
                    }
                }

                if (numPerfect >= TRAINING_LIMIT) { break; }
            }

            sreader.close();
        }

        return srs;
    }

    private boolean isWellAligned(SAMRecord sr, FastaSequenceFile2 ref) {
        boolean valid = false;
        int mismatches = 0;

        if (!sr.getReadUnmappedFlag() && sr.getCigar().numCigarElements() == 1) {
            if (!currentContig.matches(sr.getReferenceName())) {
                ref.seekToContig(sr.getReferenceName());
                currentContig = sr.getReferenceName();

                refbases = ref.nextSequence().getBases();
            }

            byte[] readbases = sr.getReadBases();
            int offset = sr.getAlignmentStart();

            if (offset + readbases.length < refbases.length) {
                valid = true;
    
                for (int i = offset, j = 0; i < offset + readbases.length; i++, j++) {
                    int refbase = BaseUtils.simpleBaseToBaseIndex((char) refbases[i - 1]);
                    int readbase = BaseUtils.simpleBaseToBaseIndex((char) readbases[j]);

                    mismatches += (refbase >= 0 && readbase >= 0 && refbase != readbase) ? 1 : 0;
                }
            }
        }

        return (valid && mismatches == 0);
    }

    private ArrayList<RawRead> loadGoodReads(HashMap<String, SAMRecord> srs, File bustardDir) {
        ArrayList<RawRead> trainingData = new ArrayList<RawRead>();
        BoundedScoringSet<RawRead> additionalData = new BoundedScoringSet<RawRead>(TRAINING_LIMIT);

        IlluminaTile tileParser = new IlluminaTile(bustardDir, lane, tile);

        int correlatedReads = 0;
        for (RawRead rr : tileParser) {
            String key1 = String.format("%d:%d:%b", rr.getXCoordinate(), rr.getYCoordinate(), false);
            String key2 = String.format("%d:%d:%b", rr.getXCoordinate(), rr.getYCoordinate(), true);

            if (srs.containsKey(key1) && srs.containsKey(key2)) {
                byte[] quals = new byte[rr.getReadLength()];
                for (int cycle = 0; cycle < rr.getReadLength(); cycle++) {
                    quals[cycle] = (byte) (BaseUtils.simpleBaseToBaseIndex((char) rr.getSequence()[cycle]) >= 0 ? 50 : 0);
                }
                rr.setQuals(quals);

                trainingData.add(rr);

                correlatedReads++;
                if (correlatedReads % (TRAINING_LIMIT < 1000 ? TRAINING_LIMIT : 1000) == 0) {
                    System.out.println("  " + correlatedReads + " intensity-correlated reads");
                }
            } else {
                additionalData.add(rr);
            }
        }
        
        tileParser.close();

        System.out.printf("  found %d perfect reads with an optional reservoir of %d good reads\n", trainingData.size(), additionalData.size());

        RawRead[] qrs = additionalData.toArray(new RawRead[0]);
        int limit = (TRAINING_LIMIT - trainingData.size() < additionalData.size()) ? (TRAINING_LIMIT - trainingData.size()) : additionalData.size();
        for (int i = 0; i < limit; i++) {
            trainingData.add(qrs[i]);
        }

        return trainingData;
    }
    */
}
