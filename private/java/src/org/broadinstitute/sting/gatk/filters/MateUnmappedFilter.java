package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.SAMRecord;

/**
 * Filter out reads that have unmapped mates
 *
 * @author chartl
 * @since 10/17/11
 */
public class MateUnmappedFilter extends ReadFilter {

    public boolean filterOut(SAMRecord read) {
        return read.getReadPairedFlag() && read.getMateUnmappedFlag();
    }
}
