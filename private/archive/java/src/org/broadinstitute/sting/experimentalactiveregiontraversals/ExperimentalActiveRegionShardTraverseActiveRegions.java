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

import net.sf.samtools.SAMFileHeader;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.executive.WindowMaker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionExtension;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActivityProfile;
import org.broadinstitute.sting.utils.activeregion.ActivityProfileResult;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

public class ExperimentalActiveRegionShardTraverseActiveRegions <M,T> extends TraversalEngine<M,T,ActiveRegionWalker<M,T>,ActiveRegionShardDataProvider> {
    /**
     * our log, which we want to capture anything from this class
     */
    protected final static Logger logger = Logger.getLogger(TraversalEngine.class);

    private final LinkedList<ActiveRegion> workQueue = new LinkedList<ActiveRegion>();
    private final LinkedList<GATKSAMRecord> myReads = new LinkedList<GATKSAMRecord>();

    @Override
    public String getTraversalUnits() {
        return "active regions";
    }

    @Override
    public T traverse( final ActiveRegionWalker<M,T> walker,
                       final ActiveRegionShardDataProvider dataProvider,
                       T sum) {
        logger.debug(String.format("ExperimentalActiveRegionShardTraverseActiveRegions.traverse: Shard is %s", dataProvider));

        ReadShardDataProvider readDataProvider = dataProvider.getReadShardDataProvider();

        final int activeRegionExtension = walker.getClass().getAnnotation(ActiveRegionExtension.class).extension();
        final int maxRegionSize = walker.getClass().getAnnotation(ActiveRegionExtension.class).maxRegion();

        final ReadView readView = new ReadView(readDataProvider);

        final List<ActiveRegion> activeRegions = new LinkedList<ActiveRegion>();
        ActivityProfile profile = new ActivityProfile(engine.getGenomeLocParser(), walker.hasPresetActiveRegions());

        Shard readShard = readDataProvider.getShard();
        SAMFileHeader header = readShard.getReadProperties().getHeader();
        WindowMaker windowMaker = new WindowMaker(readShard, engine.getGenomeLocParser(),
                readView.iterator(), readShard.getGenomeLocs(), SampleUtils.getSAMFileSamples(header));

        for(WindowMaker.WindowMakerIterator iterator: windowMaker) {
            LocusShardDataProvider locusDataProvider = dataProvider.getLocusShardDataProvider(iterator);
            final LocusView locusView = new AllLocusView(locusDataProvider);
            final LocusReferenceView referenceView = new LocusReferenceView( walker, locusDataProvider );
            ReferenceOrderedView referenceOrderedDataView = getReferenceOrderedView(walker, locusDataProvider, locusView);

            // We keep processing while the next reference location is within the interval
            GenomeLoc prevLoc = null;
            while( locusView.hasNext() ) {
                final AlignmentContext locus = locusView.next();
                final GenomeLoc location = locus.getLocation();

                if ( prevLoc != null && location.getStart() != prevLoc.getStop() + 1 ) {
                    // we've move across some interval boundary, restart profile
                    profile = incorporateActiveRegions(profile, activeRegions, activeRegionExtension, maxRegionSize);
                }

                readDataProvider.getShard().getReadMetrics().incrementNumIterations();

                // create reference context. Note that if we have a pileup of "extended events", the context will
                // hold the (longest) stretch of deleted reference bases (if deletions are present in the pileup).
                final ReferenceContext refContext = referenceView.getReferenceContext(location);

                // Iterate forward to get all reference ordered data covering this location
                final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(locus.getLocation(), refContext);

                // Call the walkers isActive function for this locus and add them to the list to be integrated later
                profile.add(walkerActiveProb(walker, tracker, refContext, locus, location));

                prevLoc = location;

                printProgress(locus.getLocation());
            }

            locusDataProvider.close();
        }

        windowMaker.close();

        updateCumulativeMetrics(readDataProvider.getShard());

        if ( ! profile.isEmpty() )
            incorporateActiveRegions(profile, activeRegions, activeRegionExtension, maxRegionSize);

        // add active regions to queue of regions to process
        // first check if can merge active regions over shard boundaries
        if( !activeRegions.isEmpty() ) {
            if( !workQueue.isEmpty() ) {
                final ActiveRegion last = workQueue.getLast();
                final ActiveRegion first = activeRegions.get(0);
                if( last.isActive == first.isActive && last.getLocation().contiguousP(first.getLocation()) && last.getLocation().size() + first.getLocation().size() <= maxRegionSize ) {
                    workQueue.removeLast();
                    activeRegions.remove(first);
                    workQueue.addLast(new ActiveRegion(last.getLocation().union(first.getLocation()), first.isActive, this.engine.getGenomeLocParser(), activeRegionExtension));
                }
            }
            workQueue.addAll( activeRegions );
        }

        logger.debug("Integrated " + profile.size() + " isActive calls into " + activeRegions.size() + " regions." );

        // now process the active regions, where possible
        boolean emptyQueue = false;
        sum = processActiveRegions(walker, sum, emptyQueue);

        return sum;
    }

    /**
     * Take the individual isActive calls and integrate them into contiguous active regions and
     * add these blocks of work to the work queue
     * band-pass filter the list of isActive probabilities and turn into active regions
     *
     * @param profile
     * @param activeRegions
     * @param activeRegionExtension
     * @param maxRegionSize
     * @return
     */
    private ActivityProfile incorporateActiveRegions(final ActivityProfile profile,
                                                     final List<ActiveRegion> activeRegions,
                                                     final int activeRegionExtension,
                                                     final int maxRegionSize) {
        if ( profile.isEmpty() )
            throw new IllegalStateException("trying to incorporate an empty active profile " + profile);

        final ActivityProfile bandPassFiltered = profile.bandPassFilter();
        activeRegions.addAll(bandPassFiltered.createActiveRegions( activeRegionExtension, maxRegionSize ));
        return new ActivityProfile( engine.getGenomeLocParser(), profile.hasPresetRegions() );
    }


    // --------------------------------------------------------------------------------
    //
    // simple utility functions
    //
    // --------------------------------------------------------------------------------

    private final ActivityProfileResult walkerActiveProb(final ActiveRegionWalker<M,T> walker,
                                          final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                          final AlignmentContext locus, final GenomeLoc location) {
        if ( walker.hasPresetActiveRegions() ) {
            return new ActivityProfileResult(location, walker.presetActiveRegions.overlaps(location) ? 1.0 : 0.0);
        } else {
            return walker.isActive( tracker, refContext, locus );
        }
    }

    private ReferenceOrderedView getReferenceOrderedView( final ActiveRegionWalker<M,T> walker,
                                                          final LocusShardDataProvider dataProvider,
                                                          final LocusView locusView) {
        if ( WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE_ORDERED_DATA )
            return new ManagingReferenceOrderedView( dataProvider );
        else
            return (RodLocusView)locusView;
    }

    // --------------------------------------------------------------------------------
    //
    // code to handle processing active regions
    //
    // --------------------------------------------------------------------------------

    private T processActiveRegions( final ActiveRegionWalker<M,T> walker, T sum, boolean emptyQueue ) {
        if( walker.activeRegionOutStream != null ) {
            writeActiveRegionsToStream(walker);
            return sum;
        } else {
            return callWalkerMapOnActiveRegions(walker, sum, emptyQueue);
        }
    }

    /**
     * Write out each active region to the walker activeRegionOutStream
     *
     * @param walker
     */
    private void writeActiveRegionsToStream( final ActiveRegionWalker<M,T> walker ) {
        // Just want to output the active regions to a file, not actually process them
        for( final ActiveRegion activeRegion : workQueue ) {
            if( activeRegion.isActive ) {
                walker.activeRegionOutStream.println( activeRegion.getLocation() );
            }
        }
    }

    private T callWalkerMapOnActiveRegions( final ActiveRegionWalker<M,T> walker, T sum, boolean emptyQueue ) {
        final int lastRegionStart = workQueue.getLast().getLocation().getStart();
        final String lastRegionContig = workQueue.getLast().getLocation().getContig();

        // If we've traversed sufficiently past the beginning of the workQueue we can unload those regions and process them
        // TODO can implement parallel traversal here
        while( workQueue.peekFirst() != null ) {
            ActiveRegion firstRegion = workQueue.getFirst();
            final String firstRegionContig = firstRegion.getLocation().getContig();
            if (emptyQueue || firstRegionContig != lastRegionContig) {
                sum = processFirstActiveRegion(sum, walker);
            }
            else {
                final int firstRegionMaxReadStop = walker.wantsExtendedReads() ? firstRegion.getMaxReadStop() : firstRegion.getExtendedMaxReadStop();
                if (lastRegionStart > firstRegionMaxReadStop) {
                    sum = processFirstActiveRegion( sum, walker );
                }
                else {
                    break;
                }
            }
        }

        return sum;
    }

    /**
     * Process the first active region and all remaining reads which overlap
     *
     * Remove the first active region from the queue
     * (NB: some reads associated with this active region may have already been processed)
     *
     * Remove all of these reads from the queue
     * (NB: some may be associated with other active regions)
     *
     * @param sum
     * @param walker
     * @return
     */
    private T processFirstActiveRegion( final T sum, final ActiveRegionWalker<M,T> walker ) {
        final ActiveRegion firstRegion = workQueue.removeFirst();

        GATKSAMRecord firstRead = myReads.peekFirst();     // don't remove because it may not be placed here
        GenomeLoc firstReadLoc = this.engine.getGenomeLocParser().createGenomeLoc( firstRead );

        while ( firstRegion.getLocation().overlapsP( firstReadLoc ) ||
                (walker.wantsExtendedReads() && firstRegion.getExtendedLoc().overlapsP( firstReadLoc ))) {
            if( firstRegion.getLocation().overlapsP( firstReadLoc ) ) {
                // The region which the highest amount of overlap is chosen as the primary region for the read (tie breaking is done as right most region)
                long maxOverlap = firstRegion.getLocation().sizeOfOverlap( firstReadLoc );
                ActiveRegion bestRegion = firstRegion;
                for( final ActiveRegion otherRegionToTest : workQueue ) {
                    if( otherRegionToTest.getLocation().sizeOfOverlap(firstReadLoc) >= maxOverlap ) {
                        maxOverlap = otherRegionToTest.getLocation().sizeOfOverlap( firstReadLoc );
                        bestRegion = otherRegionToTest;
                    }
                }
                bestRegion.add( firstRead );

                // The read is also added to all other regions in which it overlaps but marked as non-primary
                if( walker.wantsNonPrimaryReads() ) {
                    if( !bestRegion.equals(firstRegion) ) {
                        firstRegion.add(firstRead);
                    }
                    for( final ActiveRegion otherRegionToTest : workQueue ) {
                        if( !bestRegion.equals(otherRegionToTest) ) {
                            // check for non-primary vs. extended
                            if ( otherRegionToTest.getLocation().overlapsP( firstReadLoc ) ) {
                                otherRegionToTest.add( firstRead );
                            } else if ( walker.wantsExtendedReads() && otherRegionToTest.getExtendedLoc().overlapsP( firstReadLoc ) ) {
                                otherRegionToTest.add( firstRead );
                            }
                        }
                    }
                }

                // check for non-primary vs. extended
            } else if( firstRegion.getLocation().overlapsP( firstReadLoc ) ) {
                if ( walker.wantsNonPrimaryReads() ) {
                    firstRegion.add( firstRead );
                }
            } else if( walker.wantsExtendedReads() && firstRegion.getExtendedLoc().overlapsP( firstReadLoc )) {
                firstRegion.add( firstRead );
            }

            myReads.removeFirst();
            firstRead = myReads.peekFirst();
            firstReadLoc = this.engine.getGenomeLocParser().createGenomeLoc( firstRead );
        }

        logger.debug(">> Map call with " + firstRegion.getReads().size() + " " + (firstRegion.isActive ? "active" : "inactive") + " reads @ " + firstRegion.getLocation() + " with full extent: " + firstRegion.getReferenceLoc());
        final M x = walker.map( firstRegion, null );
        return walker.reduce(x, sum);
    }

    /**
     * Special function called in LinearMicroScheduler to empty out the work queue.
     * Ugly for now but will be cleaned up when we push this functionality more into the engine
     */
    public T endTraversal( final Walker<M,T> walker, T sum) {
        boolean emptyQueue = true;
        return processActiveRegions((ActiveRegionWalker<M,T>)walker, sum, emptyQueue);
    }
}
