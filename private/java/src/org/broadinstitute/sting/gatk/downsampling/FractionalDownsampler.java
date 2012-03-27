/*
 * Copyright (c) 2012, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.downsampling;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Fractional Downsampler: selects a specified fraction of the reads for inclusion
 *
 * @author David Roazen
 */
public class FractionalDownsampler implements ReadsDownsampler {

    private ArrayList<SAMRecord> selectedReads;

    private int cutoffForInclusion;

    private static final int RANDOM_POOL_SIZE = 10000;

    public FractionalDownsampler( double fraction ) {
        if ( fraction < 0.0 || fraction > 1.0 ) {
            throw new ReviewedStingException("Fraction of reads to include must be between 0.0 and 1.0, inclusive");
        }

        cutoffForInclusion = (int)(fraction * RANDOM_POOL_SIZE);
        clear();
    }

    public void submit( SAMRecord newRead ) {
        if ( GenomeAnalysisEngine.getRandomGenerator().nextInt(10000) < cutoffForInclusion ) {
            selectedReads.add(newRead);
        }
    }

    public void submit( Collection<? extends SAMRecord> newReads ) {
        for ( SAMRecord read : newReads ) {
            submit(read);
        }
    }

    public boolean hasDownsampledItems() {
        return selectedReads.size() > 0;
    }

    public Collection<SAMRecord> consumeDownsampledItems() {
        Collection<SAMRecord> downsampledItems = selectedReads;
        clear();
        return downsampledItems;
    }

    public boolean hasPendingItems() {
        return false;
    }

    public void signalEndOfInput() {
        // NO-OP
    }

    public void clear() {
        selectedReads = new ArrayList<SAMRecord>();
    }

    public boolean requiresCoordinateSortOrder() {
        return false;
    }
}
