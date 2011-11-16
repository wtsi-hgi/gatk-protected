package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 10/17/11
 * Time: 11:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class MateUnmappedFilter extends ReadFilter {

    public boolean filterOut(SAMRecord read) {
        return read.getReadPairedFlag() && read.getMateUnmappedFlag();
    }
}
