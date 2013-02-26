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

package org.broadinstitute.sting.gatk.traversals;

import com.google.java.contract.PreconditionError;
import net.sf.samtools.*;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.datasources.providers.ActiveRegionShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.LocusShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ReadShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.reads.*;
import org.broadinstitute.sting.gatk.resourcemanagement.ThreadAllocation;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.activeregion.ActiveRegionReadState;
import org.broadinstitute.sting.utils.activeregion.ExperimentalActiveRegionShardType;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.executive.WindowMaker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActivityProfileResult;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 11/13/12
 * Time: 2:47 PM
 *
 * Test the Active Region Traversal Contract
 * http://iwww.broadinstitute.org/gsa/wiki/index.php/Active_Region_Traversal_Contract
 */
public class TraverseActiveRegionsUnitTest extends BaseTest {

    private class DummyActiveRegionWalker extends ActiveRegionWalker<Integer, Integer> {
        private final double prob;
        private EnumSet<ActiveRegionReadState> states = super.desiredReadStates();

        protected List<GenomeLoc> isActiveCalls = new ArrayList<GenomeLoc>();
        protected Map<GenomeLoc, ActiveRegion> mappedActiveRegions = new HashMap<GenomeLoc, ActiveRegion>();

        public DummyActiveRegionWalker() {
            this.prob = 1.0;
        }

        public DummyActiveRegionWalker(double constProb) {
            this.prob = constProb;
        }

        public DummyActiveRegionWalker(EnumSet<ActiveRegionReadState> wantStates) {
            this.prob = 1.0;
            this.states = wantStates;
        }

        @Override
        public EnumSet<ActiveRegionReadState> desiredReadStates() {
            return states;
        }

        @Override
        public ActivityProfileResult isActive(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
            isActiveCalls.add(ref.getLocus());
            return new ActivityProfileResult(ref.getLocus(), prob);
        }

        @Override
        public Integer map(ActiveRegion activeRegion, RefMetaDataTracker metaDataTracker) {
            mappedActiveRegions.put(activeRegion.getLocation(), activeRegion);
            return 0;
        }

        @Override
        public Integer reduceInit() {
            return 0;
        }

        @Override
        public Integer reduce(Integer value, Integer sum) {
            return 0;
        }
    }

    private final TraverseActiveRegions<Integer, Integer> traverse = new TraverseActiveRegions<Integer, Integer>();
    private final ExperimentalReadShardTraverseActiveRegions<Integer, Integer> readShardTraverse = new ExperimentalReadShardTraverseActiveRegions<Integer, Integer>();
    private final ExperimentalActiveRegionShardTraverseActiveRegions<Integer, Integer> activeRegionShardTraverse = new ExperimentalActiveRegionShardTraverseActiveRegions<Integer, Integer>();

    private IndexedFastaSequenceFile reference;
    private SAMSequenceDictionary dictionary;
    private GenomeLocParser genomeLocParser;

    private List<GenomeLoc> intervals;

    private static final String testBAM = "TraverseActiveRegionsUnitTest.bam";
    private static final String testBAI = "TraverseActiveRegionsUnitTest.bai";

    private static final ExperimentalActiveRegionShardType shardType = ExperimentalActiveRegionShardType.LOCUSSHARD;

    @BeforeClass
    private void init() throws FileNotFoundException {
        reference = new CachingIndexedFastaSequenceFile(new File(hg19Reference));
        dictionary = reference.getSequenceDictionary();
        genomeLocParser = new GenomeLocParser(dictionary);

        // TODO: reads with indels
        // TODO: reads which span many regions
        // TODO: reads which are partially between intervals (in/outside extension)
        // TODO: duplicate reads
        // TODO: read at the end of a contig
        // TODO: reads which are completely outside intervals but within extension
        // TODO: test the extension itself
        // TODO: unmapped reads

        intervals = new ArrayList<GenomeLoc>();
        intervals.add(genomeLocParser.createGenomeLoc("1", 10, 20));
        intervals.add(genomeLocParser.createGenomeLoc("1", 1, 999));
        intervals.add(genomeLocParser.createGenomeLoc("1", 1000, 1999));
        intervals.add(genomeLocParser.createGenomeLoc("1", 2000, 2999));
        intervals.add(genomeLocParser.createGenomeLoc("1", 10000, 20000));
        intervals.add(genomeLocParser.createGenomeLoc("2", 1, 100));
        intervals.add(genomeLocParser.createGenomeLoc("20", 10000, 10100));
        intervals = IntervalUtils.sortAndMergeIntervals(genomeLocParser, intervals, IntervalMergingRule.OVERLAPPING_ONLY).toList();

        List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
        reads.add(buildSAMRecord("simple", "1", 100, 200));
        reads.add(buildSAMRecord("overlap_equal", "1", 10, 20));
        reads.add(buildSAMRecord("overlap_unequal", "1", 10, 21));
        reads.add(buildSAMRecord("boundary_equal", "1", 1990, 2009));
        reads.add(buildSAMRecord("boundary_unequal", "1", 1990, 2008));
        reads.add(buildSAMRecord("boundary_1_pre", "1", 1950, 2000));
        reads.add(buildSAMRecord("boundary_1_post", "1", 1999, 2050));
        reads.add(buildSAMRecord("extended_and_np", "1", 990, 1990));
        reads.add(buildSAMRecord("outside_intervals", "1", 5000, 6000));
        reads.add(buildSAMRecord("shard_boundary_1_pre", "1", 16300, 16385));
        reads.add(buildSAMRecord("shard_boundary_1_post", "1", 16384, 16400));
        reads.add(buildSAMRecord("shard_boundary_equal", "1", 16355, 16414));
        reads.add(buildSAMRecord("simple20", "20", 10025, 10075));

        createBAM(reads);
    }

    private void createBAM(List<GATKSAMRecord> reads) {
        File outFile = new File(testBAM);
        outFile.deleteOnExit();
        File indexFile = new File(testBAI);
        indexFile.deleteOnExit();

        SAMFileWriter out = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reads.get(0).getHeader(), true, outFile);
        for (GATKSAMRecord read : ReadUtils.sortReadsByCoordinate(reads)) {
            out.addAlignment(read);
        }
        out.close();
    }

    @Test
    public void testAllBasesSeen() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker();

        List<GenomeLoc> activeIntervals = getIsActiveIntervals(walker, intervals);
        // Contract: Every genome position in the analysis interval(s) is processed by the walker's isActive() call
        verifyEqualIntervals(intervals, activeIntervals);
    }

    private List<GenomeLoc> getIsActiveIntervals(DummyActiveRegionWalker walker, List<GenomeLoc> intervals) {
        List<GenomeLoc> activeIntervals = new ArrayList<GenomeLoc>();
        for (ShardDataProvider dataProvider : createDataProviders(intervals, testBAM)) {
            traverse(walker, dataProvider, 0);
            activeIntervals.addAll(walker.isActiveCalls);
        }

        return activeIntervals;
    }

    @Test (expectedExceptions = PreconditionError.class)
    public void testIsActiveRangeLow () {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(-0.1);
        getActiveRegions(walker, intervals).values();
    }

    @Test (expectedExceptions = PreconditionError.class)
    public void testIsActiveRangeHigh () {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(1.1);
        getActiveRegions(walker, intervals).values();
    }

    @Test
    public void testActiveRegionCoverage() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker();

        Collection<ActiveRegion> activeRegions = getActiveRegions(walker, intervals).values();
        verifyActiveRegionCoverage(intervals, activeRegions);
    }

    private void verifyActiveRegionCoverage(List<GenomeLoc> intervals, Collection<ActiveRegion> activeRegions) {
        List<GenomeLoc> intervalStarts = new ArrayList<GenomeLoc>();
        List<GenomeLoc> intervalStops = new ArrayList<GenomeLoc>();

        for (GenomeLoc interval : intervals) {
            intervalStarts.add(interval.getStartLocation());
            intervalStops.add(interval.getStopLocation());
        }

        Map<GenomeLoc, ActiveRegion> baseRegionMap = new HashMap<GenomeLoc, ActiveRegion>();

        for (ActiveRegion activeRegion : activeRegions) {
            for (GenomeLoc activeLoc : toSingleBaseLocs(activeRegion.getLocation())) {
                // Contract: Regions do not overlap
                Assert.assertFalse(baseRegionMap.containsKey(activeLoc), "Genome location " + activeLoc + " is assigned to more than one region");
                baseRegionMap.put(activeLoc, activeRegion);
            }

            GenomeLoc start = activeRegion.getLocation().getStartLocation();
            if (intervalStarts.contains(start))
                intervalStarts.remove(start);

            GenomeLoc stop = activeRegion.getLocation().getStopLocation();
            if (intervalStops.contains(stop))
                intervalStops.remove(stop);
        }

        for (GenomeLoc baseLoc : toSingleBaseLocs(intervals)) {
            // Contract: Each location in the interval(s) is in exactly one region
            // Contract: The total set of regions exactly matches the analysis interval(s)
            Assert.assertTrue(baseRegionMap.containsKey(baseLoc), "Genome location " + baseLoc + " is not assigned to any region");
            baseRegionMap.remove(baseLoc);
        }

        // Contract: The total set of regions exactly matches the analysis interval(s)
        Assert.assertEquals(baseRegionMap.size(), 0, "Active regions contain base(s) outside of the given intervals");

        // Contract: All explicit interval boundaries must also be region boundaries
        Assert.assertEquals(intervalStarts.size(), 0, "Interval start location does not match an active region start location");
        Assert.assertEquals(intervalStops.size(), 0, "Interval stop location does not match an active region stop location");
    }

    @Test
    public void testActiveRegionExtensionOnContig() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker();

        Collection<ActiveRegion> activeRegions = getActiveRegions(walker, intervals).values();
        for (ActiveRegion activeRegion : activeRegions) {
            GenomeLoc loc = activeRegion.getExtendedLoc();

            // Contract: active region extensions must stay on the contig
            Assert.assertTrue(loc.getStart() > 0, "Active region extension begins at location " + loc.getStart() + ", past the left end of the contig");
            int refLen = dictionary.getSequence(loc.getContigIndex()).getSequenceLength();
            Assert.assertTrue(loc.getStop() <= refLen, "Active region extension ends at location " + loc.getStop() + ", past the right end of the contig");
        }
    }

    @Test
    public void testPrimaryReadMapping() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker();

        // Contract: Each read has the Primary state in a single region (or none)
        // This is the region of maximum overlap for the read (earlier if tied)

        // simple: Primary in 1:1-999
        // overlap_equal: Primary in 1:1-999
        // overlap_unequal: Primary in 1:1-999
        // boundary_equal: Non-Primary in 1:1000-1999, Primary in 1:2000-2999
        // boundary_unequal: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_pre: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_post: Non-Primary in 1:1000-1999, Primary in 1:2000-2999
        // extended_and_np: Non-Primary in 1:1-999, Primary in 1:1000-1999, Extended in 1:2000-2999
        // outside_intervals: none
        // shard_boundary_1_pre: Primary in 1:14908-16384, Non-Primary in 1:16385-16927
        // shard_boundary_1_post: Non-Primary in 1:14908-16384, Primary in 1:16385-16927
        // shard_boundary_equal: Non-Primary in 1:14908-16384, Primary in 1:16385-16927
        // simple20: Primary in 20:10000-10100

        Map<GenomeLoc, ActiveRegion> activeRegions = getActiveRegions(walker, intervals);
        ActiveRegion region;

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1, 999));
        verifyReadMapping(region, "simple", "overlap_equal", "overlap_unequal");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1000, 1999));
        verifyReadMapping(region, "boundary_unequal", "extended_and_np", "boundary_1_pre");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 2000, 2999));
        verifyReadMapping(region, "boundary_equal", "boundary_1_post");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 14908, 16384));
        verifyReadMapping(region, "shard_boundary_1_pre");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 16385, 16927));
        verifyReadMapping(region, "shard_boundary_1_post", "shard_boundary_equal");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("20", 10000, 10100));
        verifyReadMapping(region, "simple20");
    }

    @Test
    public void testNonPrimaryReadMapping() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(
                EnumSet.of(ActiveRegionReadState.PRIMARY, ActiveRegionReadState.NONPRIMARY));

        // Contract: Each read has the Primary state in a single region (or none)
        // This is the region of maximum overlap for the read (earlier if tied)

        // Contract: Each read has the Non-Primary state in all other regions it overlaps

        // simple: Primary in 1:1-999
        // overlap_equal: Primary in 1:1-999
        // overlap_unequal: Primary in 1:1-999
        // boundary_equal: Non-Primary in 1:1000-1999, Primary in 1:2000-2999
        // boundary_unequal: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_pre: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_post: Non-Primary in 1:1000-1999, Primary in 1:2000-2999
        // extended_and_np: Non-Primary in 1:1-999, Primary in 1:1000-1999, Extended in 1:2000-2999
        // outside_intervals: none
        // shard_boundary_1_pre: Primary in 1:14908-16384, Non-Primary in 1:16385-16927
        // shard_boundary_1_post: Non-Primary in 1:14908-16384, Primary in 1:16385-16927
        // shard_boundary_equal: Non-Primary in 1:14908-16384, Primary in 1:16385-16927
        // simple20: Primary in 20:10000-10100

        Map<GenomeLoc, ActiveRegion> activeRegions = getActiveRegions(walker, intervals);
        ActiveRegion region;

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1, 999));
        verifyReadMapping(region, "simple", "overlap_equal", "overlap_unequal", "extended_and_np");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1000, 1999));
        verifyReadMapping(region, "boundary_equal", "boundary_unequal", "extended_and_np", "boundary_1_pre", "boundary_1_post");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 2000, 2999));
        verifyReadMapping(region, "boundary_equal", "boundary_unequal", "boundary_1_pre", "boundary_1_post");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 14908, 16384));
        verifyReadMapping(region, "shard_boundary_1_pre", "shard_boundary_1_post", "shard_boundary_equal");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 16385, 16927));
        verifyReadMapping(region, "shard_boundary_1_pre", "shard_boundary_1_post", "shard_boundary_equal");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("20", 10000, 10100));
        verifyReadMapping(region, "simple20");
    }

    @Test
    public void testExtendedReadMapping() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(
                EnumSet.of(ActiveRegionReadState.PRIMARY, ActiveRegionReadState.NONPRIMARY, ActiveRegionReadState.EXTENDED));

        // Contract: Each read has the Primary state in a single region (or none)
        // This is the region of maximum overlap for the read (earlier if tied)

        // Contract: Each read has the Non-Primary state in all other regions it overlaps
        // Contract: Each read has the Extended state in regions where it only overlaps if the region is extended

        // simple: Primary in 1:1-999
        // overlap_equal: Primary in 1:1-999
        // overlap_unequal: Primary in 1:1-999
        // boundary_equal: Non-Primary in 1:1000-1999, Primary in 1:2000-2999
        // boundary_unequal: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_pre: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_post: Non-Primary in 1:1000-1999, Primary in 1:2000-2999
        // extended_and_np: Non-Primary in 1:1-999, Primary in 1:1000-1999, Extended in 1:2000-2999
        // outside_intervals: none
        // shard_boundary_1_pre: Primary in 1:14908-16384, Non-Primary in 1:16385-16927
        // shard_boundary_1_post: Non-Primary in 1:14908-16384, Primary in 1:16385-16927
        // shard_boundary_equal: Non-Primary in 1:14908-16384, Primary in 1:16385-16927
        // simple20: Primary in 20:10000-10100

        Map<GenomeLoc, ActiveRegion> activeRegions = getActiveRegions(walker, intervals);
        ActiveRegion region;

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1, 999));
        verifyReadMapping(region, "simple", "overlap_equal", "overlap_unequal", "extended_and_np");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1000, 1999));
        verifyReadMapping(region, "boundary_equal", "boundary_unequal", "extended_and_np", "boundary_1_pre", "boundary_1_post");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 2000, 2999));
        verifyReadMapping(region, "boundary_equal", "boundary_unequal", "extended_and_np", "boundary_1_pre", "boundary_1_post");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 14908, 16384));
        verifyReadMapping(region, "shard_boundary_1_pre", "shard_boundary_1_post", "shard_boundary_equal");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 16385, 16927));
        verifyReadMapping(region, "shard_boundary_1_pre", "shard_boundary_1_post", "shard_boundary_equal");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("20", 10000, 10100));
        verifyReadMapping(region, "simple20");
    }

    @Test
    public void testUnmappedReads() {
        // TODO
    }

    private void verifyReadMapping(ActiveRegion region, String... reads) {
        Collection<String> wantReads = new ArrayList<String>(Arrays.asList(reads));
        for (SAMRecord read : region.getReads()) {
            String regionReadName = read.getReadName();
            Assert.assertTrue(wantReads.contains(regionReadName), "Read " + regionReadName + " assigned to active region " + region);
            wantReads.remove(regionReadName);
        }

        Assert.assertTrue(wantReads.isEmpty(), "Reads missing in active region " + region);
    }

    private Map<GenomeLoc, ActiveRegion> getActiveRegions(DummyActiveRegionWalker walker, List<GenomeLoc> intervals) {
        for (ShardDataProvider dataProvider : createDataProviders(intervals, testBAM))
            traverse(walker, dataProvider, 0);

        endTraversal(walker, 0);

        return walker.mappedActiveRegions;
    }

    private Collection<GenomeLoc> toSingleBaseLocs(GenomeLoc interval) {
        List<GenomeLoc> bases = new ArrayList<GenomeLoc>();
        if (interval.size() == 1)
            bases.add(interval);
        else {
            for (int location = interval.getStart(); location <= interval.getStop(); location++)
                bases.add(genomeLocParser.createGenomeLoc(interval.getContig(), location, location));
        }

        return bases;
    }

    private Collection<GenomeLoc> toSingleBaseLocs(List<GenomeLoc> intervals) {
        Set<GenomeLoc> bases = new TreeSet<GenomeLoc>();    // for sorting and uniqueness
        for (GenomeLoc interval : intervals)
            bases.addAll(toSingleBaseLocs(interval));

        return bases;
    }

    private void verifyEqualIntervals(List<GenomeLoc> aIntervals, List<GenomeLoc> bIntervals) {
        Collection<GenomeLoc> aBases = toSingleBaseLocs(aIntervals);
        Collection<GenomeLoc> bBases = toSingleBaseLocs(bIntervals);

        Assert.assertTrue(aBases.size() == bBases.size(), "Interval lists have a differing number of bases: " + aBases.size() + " vs. " + bBases.size());

        Iterator<GenomeLoc> aIter = aBases.iterator();
        Iterator<GenomeLoc> bIter = bBases.iterator();
        while (aIter.hasNext() && bIter.hasNext()) {
            GenomeLoc aLoc = aIter.next();
            GenomeLoc bLoc = bIter.next();
            Assert.assertTrue(aLoc.equals(bLoc), "Interval locations do not match: " + aLoc + " vs. " + bLoc);
        }
    }

    // copied from LocusViewTemplate
    protected GATKSAMRecord buildSAMRecord(String readName, String contig, int alignmentStart, int alignmentEnd) {
        SAMFileHeader header = ArtificialSAMUtils.createDefaultReadGroup(new SAMFileHeader(), "test", "test");
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        GATKSAMRecord record = new GATKSAMRecord(header);

        record.setReadName(readName);
        record.setReferenceIndex(dictionary.getSequenceIndex(contig));
        record.setAlignmentStart(alignmentStart);

        Cigar cigar = new Cigar();
        int len = alignmentEnd - alignmentStart + 1;
        cigar.add(new CigarElement(len, CigarOperator.M));
        record.setCigar(cigar);
        record.setReadString(new String(new char[len]).replace("\0", "A"));
        record.setBaseQualities(new byte[len]);

        return record;
    }

    private List<ShardDataProvider> createDataProviders(List<GenomeLoc> intervals, String bamFile) {
        GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        engine.setGenomeLocParser(genomeLocParser);
        GATKArgumentCollection arguments = new GATKArgumentCollection();
        arguments.activeRegionShardType = shardType;     // make explicit
        engine.setArguments(arguments);

        Collection<SAMReaderID> samFiles = new ArrayList<SAMReaderID>();
        SAMReaderID readerID = new SAMReaderID(new File(bamFile), new Tags());
        samFiles.add(readerID);

        SAMDataSource dataSource = new SAMDataSource(samFiles, new ThreadAllocation(), null, genomeLocParser);

        List<ShardDataProvider> providers = new ArrayList<ShardDataProvider>();

        switch (shardType) {
            case LOCUSSHARD:
                traverse.initialize(engine);
                for (Shard shard : dataSource.createShardIteratorOverIntervals(new GenomeLocSortedSet(genomeLocParser, intervals), new LocusShardBalancer())) {
                    for (WindowMaker.WindowMakerIterator window : new WindowMaker(shard, genomeLocParser, dataSource.seek(shard), shard.getGenomeLocs())) {
                        providers.add(new LocusShardDataProvider(shard, shard.getReadProperties(), genomeLocParser, window.getLocus(), window, reference, new ArrayList<ReferenceOrderedDataSource>()));
                    }
                }
                break;
            case READSHARD:
                readShardTraverse.initialize(engine);
                for (Shard shard : dataSource.createShardIteratorOverIntervals(new GenomeLocSortedSet(genomeLocParser, intervals), new ReadShardBalancer())) {
                    providers.add(new ReadShardDataProvider(shard, genomeLocParser, shard.iterator(), reference, new ArrayList<ReferenceOrderedDataSource>()));
                }
                break;
            case ACTIVEREGIONSHARD:
                activeRegionShardTraverse.initialize(engine);
                for (Shard shard : dataSource.createShardIteratorOverIntervals(new GenomeLocSortedSet(genomeLocParser, intervals), new ActiveRegionShardBalancer())) {
                    for (WindowMaker.WindowMakerIterator window : new WindowMaker(shard, genomeLocParser, dataSource.seek(shard), shard.getGenomeLocs())) {
                        providers.add(new ActiveRegionShardDataProvider(shard, shard.getReadProperties(), genomeLocParser, shard.iterator(), window.getLocus(), window, reference, new ArrayList<ReferenceOrderedDataSource>()));
                    }
                }
                break;
            default: throw new TestException("Invalid shard type");
        }

        return providers;
    }

    private void traverse(DummyActiveRegionWalker walker, ShardDataProvider dataProvider, int i) {
        switch (shardType) {
            case LOCUSSHARD:
                traverse.traverse(walker, (LocusShardDataProvider) dataProvider, i);
                break;
            case READSHARD:
                readShardTraverse.traverse(walker, (ReadShardDataProvider) dataProvider, i);
                break;
            case ACTIVEREGIONSHARD:
                activeRegionShardTraverse.traverse(walker, (ActiveRegionShardDataProvider) dataProvider, i);
                break;
            default: throw new TestException("Invalid shard type");
        }
    }

    private void endTraversal(DummyActiveRegionWalker walker, int i) {
        switch (shardType) {
            case LOCUSSHARD:
                traverse.endTraversal(walker, i);
                break;
            case READSHARD:
                readShardTraverse.endTraversal(walker, i);
                break;
            case ACTIVEREGIONSHARD:
                activeRegionShardTraverse.endTraversal(walker, i);
                break;
            default: throw new TestException("Invalid shard type");
        }
    }

}
