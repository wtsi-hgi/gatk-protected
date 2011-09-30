package org.broadinstitute.sting.gatk.walkers.newassociation;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 6/14/11
 * Time: 10:14 AM
 * To change this template use File | Settings | File Templates.
 */
public class LocusDepthProportionWalker extends LocusWalker<double[],Boolean> implements TreeReducible<Boolean> {

    @Output
    PrintStream out;

    private Map<String,Integer> samOrder;

    public void initialize() {
        final Set<String> samples = SampleUtils.getSAMFileSamples(getToolkit());
        samOrder = new HashMap<String,Integer>(samples.size());
        int idx = 0;
        out.printf("pos");
        for ( final String s : samples ) {
            out.printf("\t");
            out.printf(s);
            samOrder.put(s,idx++);
        }
        out.printf("\t%s%n","total");
    }

    public Boolean reduceInit() { return null; }

    public double[] map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if ( ref == null || ! context.hasBasePileup() || ! context.hasReads() ) { return null; }

        out.print(ref.getLocus());
        double[] props = new double[1+samOrder.size()];

        // one pass this
        int nReads = context.size();

        for ( PileupElement e : context.getBasePileup() ) {
            props[samOrder.get(e.getRead().getReadGroup().getSample())] += 1;
        }

        for ( int idx = 0; idx < props.length -1 ; idx ++ ) {
            props[idx] /= nReads;
        }

        props[props.length-1] = nReads;

        return props;
    }

    public Boolean reduce(double[] map, Boolean pr) {
        if ( map == null ) { return null; }

        StringBuffer buf = new StringBuffer();
        for ( double d : map ) {
            buf.append("\t");
            buf.append(String.format("%.4f",d));
        }

        out.printf("%s%n",buf.toString());

        return null;
    }

    public Boolean treeReduce(Boolean a, Boolean b) { return null; }
}
