package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.variantcontext.GenotypeType;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: carneiro
 * Date: 11/29/12
 * Time: 9:59 AM
 * To change this template use File | Settings | File Templates.
 */
@Requires(value={DataSource.READS, DataSource.REFERENCE})
@Allows(value={DataSource.READS, DataSource.REFERENCE})
@By(DataSource.REFERENCE)
@Reference(window=@Window(start=-200,stop=200))
public class FindMinimumCallableCoverage extends RodWalker<Integer, Integer> {

    @Output
    public PrintStream out;

    /**
     * The input callset to evaluate
     */
    @Input(fullName="alleles", shortName = "alleles", doc="The set of alleles at which to genotype", required=true)
    public RodBinding<VariantContext> alleles;


    private GATKReport report;
    private UnifiedGenotyperEngine genotyperEngine;
    private double callConf = 0.0;

    @Override
    public void initialize() {
        super.initialize();
        final UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        genotyperEngine = new UnifiedGenotyperEngine(getToolkit(), uac);
        callConf = uac.STANDARD_CONFIDENCE_FOR_CALLING;
        report = GATKReport.newSimpleReport("Position", "MinimumCallableCoverage", "EventComplexity", "VariantType", "GenotypeType");
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return 0;

        VariantContext vcComp = tracker.getFirstValue(alleles);
        if( vcComp == null )
            return 0;

        final EventInfo eventInfo = new EventInfo();
        eventInfo.position = context.getLocation();
        eventInfo.eventLength = vcComp.getAltAlleleWithHighestAlleleCount().length();
        eventInfo.variantType = vcComp.isSNP() ? VariantContext.Type.SNP : VariantContext.Type.INDEL;
        eventInfo.genotypeType = vcComp.getGenotype(0).getType();   // assumes single sample analysis

        int coverage = context.getBasePileup().depthOfCoverage();
        boolean wasCalledOnce = false;
        while (coverage > 0) {
            List<VariantCallContext> callList = genotyperEngine.calculateLikelihoodsAndGenotypes(tracker, ref, context);
            VariantCallContext call = callList.get(0) == null && callList.size() > 1 ? callList.get(1) : callList.get(0);    // if we have more than one calls and the first one is false (snp) get the second.

            if (call != null && call.isCalledAlt(callConf)) {
                context.downsampleToCoverage(--coverage);
                wasCalledOnce = true;
            }
            else {
                if (wasCalledOnce) {
                    eventInfo.coverageNeededToCall = coverage + 1;
                    eventInfo.emit();
                }
                break;
            }
        }

        return wasCalledOnce ? 0 : 1;
    }

    @Override
    public Integer reduceInit() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @java.lang.Override
    public Integer reduce(Integer value, Integer sum) {
        return sum+value;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void onTraversalDone(Integer result) {
        super.onTraversalDone(result);    //To change body of overridden methods use File | Settings | File Templates.
        report.print(out);
    }

    private class EventInfo {
        public GenomeLoc position;
        public int coverageNeededToCall;
        public int eventLength;
        public VariantContext.Type variantType;
        public GenotypeType genotypeType;

        public EventComplexity getEventComplexity() {
            if (variantType == VariantContext.Type.SNP && eventLength == 1) return EventComplexity.SIMPLE_SNP;
            if (variantType == VariantContext.Type.SNP && eventLength < 10) return EventComplexity.MEDIUM_SNP;
            if (variantType == VariantContext.Type.SNP && eventLength >= 10) return EventComplexity.LONG_SNP;
            if (variantType == VariantContext.Type.INDEL && eventLength == 1) return EventComplexity.SIMPLE_INDEL;
            if (variantType == VariantContext.Type.INDEL && eventLength < 10) return EventComplexity.MEDIUM_INDEL;
            if (variantType == VariantContext.Type.INDEL && eventLength >= 10) return EventComplexity.LONG_INDEL;

            return EventComplexity.UNKNOWN;
        }

        public void emit() {
            report.addRow(position, coverageNeededToCall, getEventComplexity(), variantType, genotypeType);
        }
    }

    private enum EventComplexity {
        SIMPLE_SNP,
        SIMPLE_INDEL,
        MEDIUM_SNP,
        MEDIUM_INDEL,
        LONG_SNP,
        LONG_INDEL,
        UNKNOWN
    }
}

