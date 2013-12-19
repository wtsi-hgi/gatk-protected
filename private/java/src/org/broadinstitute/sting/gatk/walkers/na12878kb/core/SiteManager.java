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

package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.InvalidRecordHandler;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;

import java.util.*;

/**
 * Manages the iteration of sites from the MongoDB, ensuring the records are
 * queried in the appropriate contig order.
 *
 * User: ebanks
 * Date: 6/18/13
 * Time: 10:30 AM
 */
public class SiteManager {

    // the chunks to process
    private final List<SmartIteratorChunk> chunks = new ArrayList<>();

    private final GenomeLocParser parser;

    private final class SmartIteratorChunk {

        final SiteSelector selector;
        final List<GenomeLoc> locs;
        final DBCursor cursor;

        public SmartIteratorChunk(final SiteSelector selector, final List<GenomeLoc> locs) {
            this.selector = selector;
            this.locs = locs;
            cursor = null;
        }

        public SmartIteratorChunk(final DBCursor cursor) {
            this.cursor = cursor;
            selector = null;
            locs = null;
        }
    }

    /**
     * Create a new SiteManager that iterates over the entire KB
     *
     * @param parser  genome loc parser
     */
    public SiteManager(final GenomeLocParser parser) {
        this.parser = parser;
        final SiteSelector selector = new SiteSelector(parser);
        chunks.add(new SmartIteratorChunk(selector, null));
    }

    /**
     * Create a new SiteManager that uses the given cursor.
     *
     * @param parser  genome loc parser
     * @param cursor the database cursor
     */
    protected SiteManager(final GenomeLocParser parser, final DBCursor cursor) {
        this.parser = parser;
        chunks.add(new SmartIteratorChunk(cursor));
    }

    /**
     * Create a new SiteManager that uses the given selector.  For debugging purposes only.
     *
     * @param parser  genome loc parser
     * @param selector the site selector
     */
    protected SiteManager(final GenomeLocParser parser, final SiteSelector selector) {
        this.parser = parser;
        chunks.add(new SmartIteratorChunk(selector, null));
    }

    /**
     * Create a new SiteManager that iterates over the given intervals in the order specified by the dictionary.
     *
     * @param parser  genome loc parser
     * @param allIntervals  the intervals to use, can be null (in which case iteration occurs over the entire KB)
     * @param dictionary    the master sequence dictionary to use
     */
    public SiteManager(final GenomeLocParser parser, final GenomeLocSortedSet allIntervals, final SAMSequenceDictionary dictionary) {
        this.parser = parser;

        for ( final GenomeLoc contig : GenomeLocSortedSet.createSetFromSequenceDictionary(dictionary) ) {

            final SiteSelector selector = new SiteSelector(parser);

            // if allIntervals is null then we want the entire contig
            if ( allIntervals == null ) {
                selector.addInterval(contig);
                chunks.add(new SmartIteratorChunk(selector, null));
            }
            else {
                // check to see whether any of the requested intervals lie on this contig
                final List<GenomeLoc> contigIntervals = allIntervals.getOverlapping(contig);

                // if none of them do, don't bother with it at all
                if ( contigIntervals.isEmpty() )
                    continue;

                // if there are a reasonable number of intervals, let the selector do the work in the query
                if ( ! SiteSelector.hasTooManyIntervals(contigIntervals) ) {
                    selector.addIntervals(contigIntervals);
                    chunks.add(new SmartIteratorChunk(selector, null));
                }
                // otherwise, the manager will have to keep track of the intervals
                else {
                    selector.addInterval(contig);
                    chunks.add(new SmartIteratorChunk(selector, contigIntervals));
                }
            }
        }

        if ( chunks.isEmpty() )
            throw new IllegalArgumentException("No valid intervals were requested for the manager");
    }

    /**
     * Return an iterator over the requested sites
     *
     * @param sites      the database collection to use
     * @param sortOrder  the requested sort ordering
     * @return non-null SiteIterator
     */
    public SiteIterator<MongoVariantContext> getIterator(final DBCollection sites, final DBObject sortOrder) {
        return new SmartSiteIterator<>(chunks, sites, sortOrder);
    }

    /**
     * Return an iterator over the requested sites.
     * Throws an exception if the manager was not initialized with a DBCursor.
     *
     * @return non-null SiteIterator
     */
    public SiteIterator<MongoVariantContext> getIterator() {
        for ( final SmartIteratorChunk chunk : chunks ) {
            if ( chunk.cursor == null )
                throw new IllegalStateException("Trying to iterate over managed chunks that haven't been initialized with a cursor");
        }
        return new SmartSiteIterator<>(chunks);
    }

    /**
     * Tell the SiteManager to use only reviewed sites.
     *
     * @return non-null SiteManager for only reviewed sites
     */
    public SiteManager onlyReviewed() {
        for ( final SmartIteratorChunk chunk : chunks )
            chunk.selector.onlyReviewed();
        return this;
    }

    private final class SmartSiteIterator<T extends MongoVariantContext> implements Iterable<T>, CloseableIterator<T>, SiteIterator<T>  {

        final DBCollection sites;
        final DBObject sortOrder;
        InvalidRecordHandler<T> errorHandler = null;

        // the chunks in our work queue
        final LinkedList<SmartIteratorChunk> myChunks;

        // the chunk we are currently working on
        OneChunkIterator<T> currentChunkIterator = null;

        public SmartSiteIterator(final List<SmartIteratorChunk> chunks, final DBCollection sites, final DBObject sortOrder) {
            this.sites = sites;
            this.sortOrder = sortOrder;
            this.myChunks = new LinkedList<>(chunks);
            currentChunkIterator = makeIterator(sites, sortOrder, true);
        }

        public SmartSiteIterator(final List<SmartIteratorChunk> chunks) {
            this.sites = null;
            this.sortOrder = null;
            this.myChunks = new LinkedList<>(chunks);
            currentChunkIterator = makeIterator(sites, sortOrder, true);
        }

        private OneChunkIterator<T> makeIterator(final DBCollection sites, final DBObject sortOrder, final boolean errorIfDone) {
            if ( myChunks.isEmpty() ) {
                if ( errorIfDone )
                    throw new IllegalStateException("Trying to access data from a closed stream or one with no valid intervals");
                return null;
            }

            final SmartIteratorChunk chunk = myChunks.removeFirst();
            final DBCursor cursor = chunk.cursor != null ? chunk.cursor : sites.find(chunk.selector.toQuery()).sort(sortOrder);
            return new OneChunkIterator<>(parser, cursor, chunk.locs, errorHandler);
        }

        @Override
        public void setErrorHandler(final InvalidRecordHandler<T> errorHandler) {
            this.errorHandler = errorHandler;
            if ( currentChunkIterator != null )
                currentChunkIterator.setErrorHandler(errorHandler);
        }

        private void assureActiveStream() {
            if ( currentChunkIterator == null )
                throw new IllegalStateException("Trying to access data from a closed stream");

            // even if there's no data left don't remove the last iterator (since it's legal to have a used up stream)
            while ( ! currentChunkIterator.hasNext() && ! myChunks.isEmpty() ) {
                currentChunkIterator.close();
                currentChunkIterator = makeIterator(sites, sortOrder, false);
            }
        }

        @Override
        public List<T> getSitesBefore(final GenomeLoc loc) {
            assureActiveStream();
            return currentChunkIterator.getSitesBefore(loc);
        }

        @Override
        public List<T> getSitesAtLocation(final GenomeLoc loc) {
            assureActiveStream();
            return currentChunkIterator.getSitesAtLocation(loc);
        }

        @Override
        public List<T> getNextEquivalents() {
            assureActiveStream();
            return currentChunkIterator.getNextEquivalents();
        }

        @Override
        public List<T> toList() {
            final List<T> l = new LinkedList<>();

            while ( hasNext() )
                l.addAll(currentChunkIterator.toList());

            return l;
        }

        @Override
        public boolean hasNext() {
            if ( currentChunkIterator == null )
                return false;

            if ( ! currentChunkIterator.hasNext() ) {
                currentChunkIterator.close();
                currentChunkIterator = makeIterator(sites, sortOrder, false);
                return hasNext(); // recursively call on the next chunk
            }

            return true;
        }

        @Override
        public T next() {
            if ( currentChunkIterator == null )
                throw new IllegalStateException("Trying to access data from a closed stream");
            return currentChunkIterator.next();
        }

        @Override
        public Iterator<T> iterator() {
            return this;
        }

        @Override
        public void close() {
            if ( currentChunkIterator != null )
                currentChunkIterator.close();
        }

        /** Unsupported Operation. */
        public void remove() {
            throw new UnsupportedOperationException("Not supported: remove");
        }
    }
}
