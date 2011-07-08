package org.broadinstitute.sting.gatk.filters;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 6/27/11
 * Time: 2:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class AddAberrantInsertTagFilter extends ReadFilter {
    @Argument(shortName="aif",fullName = "aberrantInsertPercentileFile",required=false, doc="File containing aberrant insert size cutoffs per read group")
    File aif = null;

    HashMap<String,Pair<Integer,Integer>> rgCuts;
    public void initialize(GenomeAnalysisEngine engine) {
        rgCuts = new HashMap<String,Pair<Integer,Integer>>();
        XReadLines xrl;
        try {
            if ( aif != null  ) {
            xrl = new XReadLines(aif);
            } else {
                xrl = null;
            }
        } catch (FileNotFoundException e) {
            throw new ReviewedStingException("File not found: "+aif.getName(),e);
        }

        if ( xrl == null ) {
            rgCuts = null;
        } else {
            for ( String line : xrl ) {
                String[] args = line.split("\\s+");
                int lower = Integer.parseInt(args[1]);
                int upper = Integer.parseInt(args[2]);
                String rg = args[0];

                rgCuts.put(rg,new Pair<Integer,Integer>(lower,upper));
            }
        }
    }

    public boolean filterOut(SAMRecord read) {
        if ( rgCuts != null && ! read.getMateUnmappedFlag() ) {
            int iSize = Math.abs(read.getInferredInsertSize());
            Pair<Integer,Integer> threshs = rgCuts.get(read.getReadGroup().getId());
            boolean isAb = iSize < threshs.first || iSize > threshs.second;
            read.setAttribute("AI",isAb ? 1 : 0);
        }

        return false;
    }
}
