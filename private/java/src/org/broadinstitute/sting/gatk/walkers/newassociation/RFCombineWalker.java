package org.broadinstitute.sting.gatk.walkers.newassociation;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.table.TableFeature;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/23/11
 * Time: 5:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class RFCombineWalker extends RodWalker<Object,Object> {

    private static final String FIRST_COL = "chrm:start-stop";

    @Output
    PrintStream out;

    @Argument(fullName = "RFAOutput", shortName = "r", doc="Outputs from RFA walker")
    public List<RodBinding<TableFeature>> rfaOutputs = Collections.emptyList();

    private List<String> order;

    private GenomeLoc prevLoc;
    private boolean printHeader;

    public void initialize() {
        printHeader = true;
        prevLoc = null;
    }

    public Object reduceInit() { return null; }

    public Object map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || ref == null ) { return null;}

        GenomeLoc loc = null;
        boolean needPrint = false;
        List<String> eventBySample = new ArrayList<String>();
        if ( printHeader ) {
            for ( RodBinding<TableFeature> bound : rfaOutputs ) {
                StringBuffer outputHeader = new StringBuffer(FIRST_COL);
                for ( String entry : tracker.getFirstValue(bound).getHeader() ) {
                    if ( ! entry.startsWith("HEADER") ) {
                        outputHeader.append("\t");
                        outputHeader.append(entry);
                    }
                }
                out.printf("%s%n",outputHeader);
            }
        }

        for ( RodBinding<TableFeature> bound : rfaOutputs ) {
            TableFeature feature = tracker.getFirstValue(bound);
            loc = feature.getLocation();
            if ( prevLoc != null && loc.equals(prevLoc) ) {
                break;
            }

            for ( String s : feature.getAllValues().subList(1,feature.getAllValues().size()) ) {
                boolean has = ! ( s.charAt(0) == '0' );
                eventBySample.add(s);
                needPrint |= has;
            }
        }

        if ( needPrint && (loc != null)) {
            out.printf("%s",loc.toString());
            for ( String s : eventBySample ) {
                out.printf("\t%s", s);
            }
            out.printf("%n");
        }

        prevLoc = loc;

        return null;
    }

    public Object reduce(Object map, Object reduce) { return null; }
}
