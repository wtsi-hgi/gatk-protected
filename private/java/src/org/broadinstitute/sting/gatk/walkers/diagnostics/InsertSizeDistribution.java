package org.broadinstitute.sting.gatk.walkers.diagnostics;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.*;

/**
 * This tool computes the insert size distributions (from the ISIZE field of the GATKSAMRecord) for each sample and read group in a BAM
 *
 * <h2>Input</h2>
 * <p>
 *     Any number of BAM files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 *     Emits a standard <a href="https://www.broadinstitute.org/gsa/wiki/index.php/GATKReport">GATKReport</a>
 *     with two tables:
 *
 *     <dl>
 *         <dt>InsertSizeDistributionBySample</dt>
 *         <dd>Table of read counts with each insert size, for each sample</dd>

 *         <dt>InsertSizeDistributionByReadGroup</dt>
 *         <dd>Table of read counts with each insert size, for each read group</dd>
 *     </dl>
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *     -T InsertSizeDistribution -I my.bam -o insert_sizes.gatkreport.txt
 * </pre>

 * @author Kiran Garimella and Mark DePristo
 * @since 2010-2011
 */
public class InsertSizeDistribution extends ReadWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    private GATKReport report;

    @Argument(doc="Maximum insert size (can't be bigger than an integer)",required=false,shortName="mis",fullName="maxInsertSize")
    public int maxInsert = Integer.MAX_VALUE;

    Map<String,Map<Integer,Integer>> rgCounts;
    Map<String,Map<Integer,Integer>> samCounts;

    public void initialize() {
        HashSet<String> samplesSeen = new HashSet<String>();
        HashSet<String> rgsSeen = new HashSet<String>();

        for (SAMReadGroupRecord rg : this.getToolkit().getSAMFileHeader().getReadGroups()) {
            samplesSeen.add(rg.getSample());
            rgsSeen.add(rg.getReadGroupId());
        }

        rgCounts = new HashMap<String,Map<Integer,Integer>>(rgsSeen.size());
        samCounts = new HashMap<String,Map<Integer,Integer>>(samplesSeen.size());

        for ( String s : samplesSeen ) {
            samCounts.put(s,new HashMap<Integer,Integer>(1200));
        }

        for ( String s : rgsSeen ) {
            rgCounts.put(s, new HashMap<Integer,Integer>(1200));
        }
    }

    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
        return (read.getReadPairedFlag() && read.getFirstOfPairFlag());
    }

    @Override
    public Integer map(ReferenceContext referenceContext, GATKSAMRecord samRecord, RefMetaDataTracker RefMetaDataTracker) {

        int insert = Math.abs(samRecord.getInferredInsertSize());
        if ( insert > maxInsert )
            insert = maxInsert;
        final String rg = samRecord.getReadGroup().getReadGroupId();
        final String sample = samRecord.getReadGroup().getSample();

        Map<Integer,Integer> sMap = samCounts.get(sample);
        Map<Integer,Integer> rgMap = rgCounts.get(rg);

        if ( ! sMap.containsKey(insert) ) {
            sMap.put(insert,0);
        }

        if ( ! rgMap.containsKey(insert) ) {
            rgMap.put(insert,0);
        }

        sMap.put(insert,sMap.get(insert)+1);
        rgMap.put(insert,rgMap.get(insert)+1);

        return null;
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer integer, Integer integer1) {
        return null;
    }

    public void onTraversalDone(Integer sum) {
        report = new GATKReport();
        report.addTable("InsertSizeDistributionBySample", "Table of insert size distributions", 1 + samCounts.size());
        report.addTable("InsertSizeDistributionByReadGroup", "Table of insert size distributions", 1 + rgCounts.size());

        GATKReportTable samTable = report.getTable("InsertSizeDistributionBySample");
        samTable.addColumn("Sample","%s");
        samTable.addColumn("InsertSize","%d");
        samTable.addColumn("Count","%d");
        GATKReportTable rgTable = report.getTable("InsertSizeDistributionByReadGroup");
        rgTable.addColumn("ReadGroup","%s");
        rgTable.addColumn("InsertSize","%d");
        rgTable.addColumn("Count","%d");

        // be explicit about inserts with 0 observations
        for ( Map.Entry<String,Map<Integer,Integer>> sampleInsertSize : samCounts.entrySet() ) {
            if ( sampleInsertSize.getValue().size() == 0 ) {
                logger.warn("Sample "+sampleInsertSize.getKey()+" present in header, but no reads counted");
                continue;
            }
            Integer maxInsert = Collections.max(new ArrayList<Integer>(sampleInsertSize.getValue().keySet()));
            for ( Integer insert = 0; insert < maxInsert; insert++ ) {
                String rid = sampleInsertSize.getKey()+"."+insert.toString();
                samTable.addRowID(rid);
                samTable.set(rid,"Sample",sampleInsertSize.getKey());
                samTable.set(rid,"InsertSize",insert);
                Integer count = sampleInsertSize.getValue().get(insert);
                if ( count == null )
                    count = 0;
                samTable.set(rid,"Count",count);
            }
        }

        for ( Map.Entry<String,Map<Integer,Integer>> rgInsertSize : rgCounts.entrySet() ) {
            if ( rgInsertSize.getValue().size() == 0 ) {
                logger.warn("Read group "+rgInsertSize.getKey()+" present in header, but no reads counted");
                continue;
            }
            Integer maxInsert = Collections.max(new ArrayList<Integer>(rgInsertSize.getValue().keySet()));
            for ( Integer insert = 0; insert < maxInsert; insert++ ) {
                String rid = rgInsertSize.getKey()+"."+insert.toString();
                rgTable.addRowID(rid);
                rgTable.set(rid,"ReadGroup",rgInsertSize.getKey());
                rgTable.set(rid,"InsertSize",insert);
                Integer count = rgInsertSize    .getValue().get(insert);
                if ( count == null )
                    count = 0;
                rgTable.set(rid,"Count",count);
            }
        }
        // Write report.
        report.print(out);
    }
}
