package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.variant.variantcontext.GenotypeType;
import org.broadinstitute.variant.variantcontext.VariantContext;

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

    @Argument(fullName = "bootstrap", shortName = "boot", doc = "Number of bootstrap interations", required = false)
    public int bootstrapIterations = 100;

    @Hidden
    @Argument(fullName = "debugLevel", shortName = "dl", doc = "output debug information", required = false)
    public int debugLevel = 0;

    private GATKReport report;
    private UnifiedGenotyperEngine snpEngine;
    private UnifiedGenotyperEngine indelEngine;
    private double callConf = 0.0;
    private long mapCounter = 1;

    public int binarySearch(int left, int right, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariantContext.Type variantType) {
        int result;
        if (left < right) {
            int coverage = (int) Math.floor((left+right)/2);
            if (canCall(tracker, ref, context, variantType, coverage)) {
                result = binarySearch(left, coverage-1, tracker, ref, context, variantType);
            }
            else {
                result = binarySearch(coverage+1, right, tracker, ref, context, variantType);
            }
        }
        else if (left > right) {
            result = left;
        }
        else {
            result =  canCall(tracker, ref, context, variantType, left) ? left : left+1;
        }
        return result;
    }

    public boolean canCall(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariantContext.Type variantType, int downsampleTo) {
        AlignmentContext dsContext;
        if (downsampleTo < context.getBasePileup().depthOfCoverage()) {
            dsContext = new AlignmentContext(context.getLocation(), context.getBasePileup().copy(), context.getSkippedBases(), context.hasPileupBeenDownsampled());
            dsContext.downsampleToCoverage(downsampleTo);
        }
        else {
            dsContext = context;             // avoids unnecessary downsampling
        }

        final List<VariantCallContext> callList = variantType == VariantContext.Type.SNP ?
                snpEngine.calculateLikelihoodsAndGenotypes(tracker, ref, dsContext) :
                indelEngine.calculateLikelihoodsAndGenotypes(tracker, ref, dsContext);
        final VariantCallContext call = callList.isEmpty() ? null : callList.get(0);

        return (call != null && call.isCalledAlt(callConf) && call.getType() == variantType);
    }

    @Override
    public void initialize() {
        super.initialize();
        final UnifiedArgumentCollection snpUAC = new UnifiedArgumentCollection();
        final UnifiedArgumentCollection indelUAC = new UnifiedArgumentCollection();
        snpUAC.GLmodel = GenotypeLikelihoodsCalculationModel.Model.SNP;
        indelUAC.GLmodel = GenotypeLikelihoodsCalculationModel.Model.INDEL;
        snpEngine = new UnifiedGenotyperEngine(getToolkit(), snpUAC);
        indelEngine = new UnifiedGenotyperEngine(getToolkit(), indelUAC);
        callConf = snpUAC.STANDARD_CONFIDENCE_FOR_CALLING;
        report = GATKReport.newSimpleReport("MinCov", "Position", "MinimumCallableCoverage", "EventComplexity", "VariantType", "GenotypeType");
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

        int originalCoverage = context.getBasePileup().depthOfCoverage();
        boolean wasCalledOnce = false;
        int sum = 0;
        int result = 0;

        if (debugLevel > 0) System.out.print(mapCounter++ + "(" + context.getLocation() + "): ");

        if (canCall(tracker, ref, context, eventInfo.variantType, originalCoverage)) {


            for (int i=0; i<bootstrapIterations; i++) {
                int minCoverage = binarySearch(0, originalCoverage, tracker, ref, context, eventInfo.variantType);
                sum += minCoverage;
                if (debugLevel > 0) System.out.print(minCoverage + ", ");
            }

            eventInfo.coverageNeededToCall = sum/bootstrapIterations;
            eventInfo.emit();
            if (debugLevel > 0) System.out.print("[" + eventInfo.coverageNeededToCall + "]");
            result = 1;

            if (debugLevel > 0) System.out.println();
        }
        return result;

    }

    @Override
    public Integer reduceInit() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
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

