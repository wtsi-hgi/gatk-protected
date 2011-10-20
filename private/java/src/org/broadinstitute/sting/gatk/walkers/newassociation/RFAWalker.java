package org.broadinstitute.sting.gatk.walkers.newassociation;

import cern.jet.math.Arithmetic;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.ExactAFCalculationModel;
import org.broadinstitute.sting.gatk.walkers.newassociation.features.old.BinaryFeatureAggregator;
import org.broadinstitute.sting.gatk.walkers.newassociation.features.old.ReadFeatureAggregator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.MathUtils;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Read feature association walker -- associates read features between dichotomized, or multi-group cohorts
 */

@ReadFilters({MaxInsertSizeFilter.class,MappingQualityFilter.class,DuplicateReadFilter.class,
        FailsVendorQualityCheckFilter.class,NotPrimaryAlignmentFilter.class,UnmappedReadFilter.class,
        AddAberrantInsertTagFilter.class})
public abstract class RFAWalker extends ReadWalker<SAMRecord,RFWindow> {
    /*    // todo -- marginalization over genotypes and tables
    // todo -- segmenting the window to find max concentration of reads
    // todo -- estimation of purity of the actual window selection for the event
    @ArgumentCollection
    private RFAArgumentCollection collection = new RFAArgumentCollection();

    @Output
    public PrintStream out;

    Map<String,Boolean> caseStatus;

    protected List<BinaryFeatureAggregator> aggregators; // no re-instantiation, use a list to ensure ordering

    protected Iterator<GenomeLoc> locusIterator;
    protected GenomeLoc iteratorLoc;
    protected GenomeLoc loc;
    protected String sample;
    private List<String> EMPTY_LIST = new ArrayList<String>(0);

    private short nCase;
    private short nControl;
    private int[] caseR;
    private int[] controlR;
    private double[] caseA;
    private double[] controlA;

    private RFAGenotypeLikelihoodsCalculationModel rfaglcm = new RFAGenotypeLikelihoodsCalculationModel();

    public void initialize() {
        if ( collection.windowSize % collection.windowJump != 0 ) {
            throw new UserException("Window size is not divisible by window jump.");
        }

        if ( collection.caseFile == null || collection.controlFile == null ) {
            throw new UserException("You must provide both a case file (-case) and a control file (-control) each listing those samples belonging to the cohort");
        }

        caseStatus = new HashMap<String,Boolean>();
        nCase = 0;
        nControl = 0;
        try {
            for ( String sample : new XReadLines(collection.caseFile) ) {
                caseStatus.put(sample,true);
                ++nCase;
            }
            for ( String sample : new XReadLines(collection.controlFile)) {
                caseStatus.put(sample,false);
                ++nControl;
            }

            for ( final String sample : SampleUtils.getSAMFileSamples(getToolkit()) ) {
                if ( ! caseStatus.containsKey(sample)) {
                    throw new UserException("No case/control status for sample "+sample);
                }
            }

        } catch ( FileNotFoundException e ) {
            throw new UserException("Unable to open a case/control file",e);
        }

        Set<Class<? extends BinaryFeatureAggregator>> aggregatorSet = getFeatureAggregators(collection.inputFeatures);
        Set<BinaryFeatureAggregator> rfHolder1 = new HashSet<BinaryFeatureAggregator>(aggregatorSet.size());
        try {
            for ( Class<? extends BinaryFeatureAggregator> featureClass : aggregatorSet ) {
                BinaryFeatureAggregator readFeature = featureClass.getConstructor(RFAArgumentCollection.class).newInstance(collection);
                rfHolder1.add(readFeature);
            }
        } catch ( Exception e ) {
            throw new StingException("A read feature instantiation error occurred during initialization",e);
        }

        BinaryFeatureAggregator[] rfHolder2 = new BinaryFeatureAggregator[rfHolder1.size()];
        int idx = 0;
        for ( BinaryFeatureAggregator f : rfHolder1 ) {
            rfHolder2[idx++] = f;
        }
        Arrays.sort(rfHolder2, new Comparator<BinaryFeatureAggregator>() {
            public int compare(BinaryFeatureAggregator a, BinaryFeatureAggregator b) {
                return a.getClass().getSimpleName().compareTo(b.getClass().getSimpleName());
            }
        });
        aggregators = Arrays.asList(rfHolder2);

        writeHeader();

        locusIterator = getToolkit().getIntervals().iterator();
        iteratorLoc = locusIterator.hasNext() ? locusIterator.next() : null;

        caseR = new int[nCase];
        controlR = new int[nControl];
        caseA = new double[nCase];
        controlA = new double[nControl];
    }

    public RFWindow reduceInit() {
        Set<String> samples = new HashSet<String>(getSampleDB().getSamples().size());
        for ( Sample s : getSampleDB().getSamples() ) {
            samples.add(s.getID());
        }
        return new RFWindow(aggregators,collection,caseStatus,getToolkit().getGenomeLocParser());
    }

    public SAMRecord map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        if ( ref == null ) { return null; } // unmapped reads have null ref contexts
        //loc = getToolkit().getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(),read.getAlignmentStart());
        GenomeLoc newLoc = ref.getLocus().getStartLocation(); // can be problematic if read aligns prior to start of contig -- should never happen
        if ( newLoc.isPast(iteratorLoc.getStartLocation()) ) {
            loc = newLoc;
        } else {
            loc = iteratorLoc.getStartLocation();
        }
        if ( read == null ) { return null; }
        sample = read.getReadGroup().getSample();
        return read;
    }

    public RFWindow reduce(SAMRecord read, RFWindow prevReduce) {
        if ( iteratorLoc != null && iteratorLoc.isBefore(loc) ) {// test if read is past end of the user interval
            //logger.info(String.format("iteratorLoc: %s    loc: %s",iteratorLoc.toString(),loc.toString()));
            onIntervalDone(prevReduce);
            iteratorLoc = locusIterator.hasNext() ? locusIterator.next() : null;
            if ( loc.startsBefore(iteratorLoc) ) {
                loc = iteratorLoc.getStartLocation();
            }
            reduce(read,prevReduce);
        } else if ( read != null ) {
            // todo -- what happens if first read of an interval is not before or at the start of the interval?\
            List<Pair<GenomeLoc,Map<String,List<BinaryFeatureAggregator>>>> completed = prevReduce.inc(read, loc, sample,iteratorLoc);
            // todo -- run tests here; for now just log that a window/multiple windows are complete
            if ( completed.size() > 0 ) {
                // System.out.printf("At %s we have seen %d completed windows%n",loc,completed.size())
                // bed format
                int locShift = 0;
                for ( Pair<GenomeLoc,Map<String,List<BinaryFeatureAggregator>>> samWindow : completed ) {
                    GenomeLoc window = samWindow.first;
                    runWindowTests(samWindow.second, window);
                    locShift += collection.windowJump;
                }
            }
        }
        return prevReduce;
    }


    public RFWindow onIntervalDone(RFWindow rWindow) {
        //logger.info("In onIntervalDone at genome loc "+iteratorLoc.toString()+" with read loc "+loc.toString());
        List<Pair<GenomeLoc,Map<String,List<BinaryFeatureAggregator>>>> completed = rWindow.flush(iteratorLoc);
        int locShift = 0;
        for ( Pair<GenomeLoc,Map<String,List<BinaryFeatureAggregator>>> samWindow : completed ) {
            GenomeLoc window = samWindow.first;
            runWindowTests(samWindow.second, window);
            locShift += collection.windowJump;
        }

        return rWindow;
    }

    public Set<Class<? extends BinaryFeatureAggregator>> getFeatureAggregators(List<String> requestedFeatures) {
        HashSet<Class<? extends BinaryFeatureAggregator>> newFeatureSet = new HashSet<Class<? extends BinaryFeatureAggregator>>();
        List<Class<? extends BinaryFeatureAggregator>> availableFeatures = new PluginManager<BinaryFeatureAggregator>(BinaryFeatureAggregator.class).getPlugins();

        if ( collection.inputFeatures == null ) {
            newFeatureSet.addAll(availableFeatures);
            return newFeatureSet;
        }


        Map<String,Class<? extends BinaryFeatureAggregator>> classNameToClass = new HashMap<String,Class<? extends BinaryFeatureAggregator>>(collection.inputFeatures.size());
        for ( Class<? extends BinaryFeatureAggregator> clazz : availableFeatures ) {
            classNameToClass.put(clazz.getSimpleName(),clazz);
        }

        for ( String s : requestedFeatures) {
            if ( classNameToClass.containsKey(s) ) {
                newFeatureSet.add(classNameToClass.get(s));
            } else {
                throw new UserException("The name "+s+" does not correspond to an available read feature class.");
            }
        }

        return newFeatureSet;
    }

    public void runWindowTests(Map<String,List<BinaryFeatureAggregator>> window, GenomeLoc loc) {
        // two main tests: fixed-significance shift, and genotype-free skew
        // todo -- really the aggregators should be iterated over directly (rather than indirectly through the index)
        out.printf("%s\t%d\t%d",loc.getContig(),loc.getStart(),loc.getStop());
        for ( int agIdx = 0; agIdx < aggregators.size(); agIdx ++ ) {
            double fixedDelta = fixedSignificance(window.get("case").get(agIdx),window.get("control").get(agIdx));
            double genotypeFreePVal = fixedDelta > 0 ? calculateGenotypeFreeSkew(window,agIdx) : 1.0 ;
            out.printf("\t%.2e\t%.2e",fixedDelta,genotypeFreePVal);

        }
        out.printf("%n");
    }

    public double fixedSignificance(BinaryFeatureAggregator caseAg, BinaryFeatureAggregator controlAg) {
        if ( caseAg.getnReads() == 0 || controlAg.getnReads() == 0 ) {
            return 0.0;
        }
        double stat_num = caseAg.getMean() - controlAg.getMean();
        double stat_denom = Math.sqrt(caseAg.getUnbiasedVar()/caseAg.getnReads() + controlAg.getUnbiasedVar()/controlAg.getnReads());
        double stat = stat_num/stat_denom;
        if ( ! Double.isNaN(stat) && stat*stat < collection.fixedZ*collection.fixedZ ) {
            return 0.0;
        } else {
	    return stat_num;
        }
    }

    private double calculateGenotypeFreeSkew(Map<String,List<BinaryFeatureAggregator>> window, int offset) {
        Map<String,Genotype> caseGenotypeLikelihoods = rfaglcm.getLikelihoods(null,null,unwrap(window,offset,true));
        Map<String,Genotype> controlGenotypeLikelihoos = rfaglcm.getLikelihoods(null,null,unwrap(window,offset,false));
        double[] caseSpectrum = new double[nCase];
        double[] controlSpectrum = new double[nControl];
        // todo -- generalize to floating ploidy
        ExactAFCalculationModel.linearExact(caseGenotypeLikelihoods,caseSpectrum.clone(),caseSpectrum,0,1,2);
        ExactAFCalculationModel.linearExact(controlGenotypeLikelihoos,controlSpectrum.clone(),controlSpectrum,0,1,2);
        return MathUtils.marginalizedFisherExact(caseSpectrum,controlSpectrum,nCase,nControl);
    }

    private Map<String,BinaryFeatureAggregator> unwrap(Map<String,List<BinaryFeatureAggregator>> map, int o, boolean getCase) {
        Map<String,BinaryFeatureAggregator> toRet = new HashMap<String,BinaryFeatureAggregator>(map.size());
        for ( Map.Entry<String,List<BinaryFeatureAggregator>> entry : map.entrySet() ) {
            if ( entry.getKey().equals("case") || entry.getKey().equals("control") ) {
                continue;
            }

            if ( caseStatus.get(entry.getKey()).equals(getCase) ) {
                toRet.put(entry.getKey(),entry.getValue().get(o));
            }
        }

        return toRet;
    }

    public void writeHeader() {
        // "%.2e\t%d:%d\t%s,%s\t%.2e"
        StringBuffer buf = new StringBuffer();
        buf.append("description=chr,start,stop");
        for ( BinaryFeatureAggregator f : aggregators ) {
            buf.append(",");
            buf.append(f.getClass().getSimpleName());
            buf.append("-d,");
            buf.append(f.getClass().getSimpleName());
            buf.append("-p");
        }
        out.printf("track type=bedTable %s%n",buf);
    }
    */
}

/*    DEAD CODE BELOW HERE

    private double calculateGenotypeFreeSkew(Map<String,List<BinaryFeatureAggregator>> window, int rfaIndex) {
        // two-step process. Step 1: identify the best simple state space (e.g. integer N) that can represent discrete
        // genotypes within the data. N can vary because of subclonal populations or mixtures of CNVs in the data. For now
        // assume N is the same between cases and controls.

        // lots of iteration, so cast these to arrays first
        int caseidx = 0;
        int controlidx = 0;
        for ( Map.Entry<String,Boolean> sample : caseStatus.entrySet() ) {
            if ( sample.getValue() ) {
                caseR[caseidx] = window.get(sample.getKey()).get(rfaIndex).getnReads();
                caseA[caseidx] = window.get(sample.getKey()).get(rfaIndex).getMean();
                caseidx++;
            } else {
                controlR[controlidx] = window.get(sample.getKey()).get(rfaIndex).getnReads();
                controlA[controlidx] = window.get(sample.getKey()).get(rfaIndex).getMean();
                controlidx++;
            }
        }

        // find the maximum possible likelihood of the data under the binomial hypothesis by setting each samples 'frequency'
        // to the observed frequency

        double L_Star = 0;
        for ( int i = 0; i < nCase; i ++ ) {
            L_Star += MathUtils.log10BinomialProbability(caseR[i],(int)(0.5+caseR[i]*caseA[i]),caseA[i]);
        }
        for ( int j = 0; j < nControl; j++ ) {
            L_Star += MathUtils.log10BinomialProbability(controlR[j],(int)(0.5+controlR[j]*controlA[j]),controlA[j]);
        }

        // calculate the base likelihood - binomial
        // note that this is a heuristic proxy: likelihood is only calculated at the closest points (and not summed over all)
        // this heuristic is used for speed.
        double bestlikelihood = Double.NEGATIVE_INFINITY;
        int bestN = 2;
        int N = 1;
        do {
            double likelihood = 0;
            logger.debug(String.format("likelihood:%.2f%n",likelihood));
            N++;
            for ( int i = 0; i < nCase; i++ ) {
                // find the best k
                int k = (int) (0.5+N*caseA[i]);
                likelihood += MathUtils.log10BinomialProbability(caseR[i],(int)(0.5+caseR[i]*caseA[i]),Math.max(((double) k)/N,0.01));
            }

            for ( int j = 0; j < nControl; j++ ) {
                int k = (int) (0.5+N*controlA[j]);
                likelihood += MathUtils.log10BinomialProbability(controlR[j],(int)(0.5+controlR[j]*controlA[j]),Math.max(((double) k)/N,0.01));
            }

            if ( 2*(likelihood - bestlikelihood) > 3.8414 ) {
                bestlikelihood = likelihood;
                bestN = N;
            }

        } while ( N <= 10 && 2*(L_Star - bestlikelihood) > 3.8414);

        return 1.0;
    }

    private double calculateGL(int indexToSwap, int newK, int bestN) {
        double logLik = 0;
        int idx = 0;
        for ( int i = 0; i < nCase; i++ ) {
            int k = idx == indexToSwap ? newK : (int) (0.5+bestN*caseA[i]);
            double bprob = MathUtils.log10BinomialProbability(caseR[i],(int)(0.5+caseR[i]*caseA[i]),Math.log10(rectify(((double) k)/bestN)));
            //logger.debug(String.format("N: %d K: %d P: %.2f Prob: %.2f",caseR[i],(int)(0.5+caseR[i]*caseA[i]),rectify(((double) k)/bestN),bprob));
            logLik += bprob;
            idx++;
        }

        for ( int i = 0; i < nControl; i++ ) {
            int k = idx == indexToSwap ? newK : (int) (0.5+bestN*controlA[i]);
            logLik += MathUtils.log10BinomialProbability(controlR[i],(int)(0.5+controlR[i]*controlA[i]),Math.log10(rectify(((double) k)/bestN)));
            idx++;
        }

        return logLik;
    }


    private double rectify(double p) {
        return Math.min(Math.max(p,0.005),0.995);
    }

 */
