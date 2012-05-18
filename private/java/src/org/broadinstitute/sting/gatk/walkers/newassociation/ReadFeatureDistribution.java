package org.broadinstitute.sting.gatk.walkers.newassociation;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.newassociation.features.ReadFeature;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.*;

public class ReadFeatureDistribution extends ReadWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @ArgumentCollection
    RFAArgumentCollection collection = new RFAArgumentCollection();

    private GATKReport report;

    private Set<ReadFeature> featuresToExtract;

    public void initialize() {
        report = new GATKReport();
        Set<Class<? extends ReadFeature>> aggregatorSet = getFeatureAggregators(collection.inputFeatures);
        featuresToExtract = new HashSet<ReadFeature>();
        try {
            List<SAMReadGroupRecord> readGroups = this.getToolkit().getSAMFileHeader().getReadGroups();
            for ( Class<? extends ReadFeature> featureClass : aggregatorSet ) {
                ReadFeature readFeature = featureClass.getConstructor(RFAArgumentCollection.class).newInstance(collection);
                report.addTable(readFeature.getName(),readFeature.getDescription(), 1 + readGroups.size());
                GATKReportTable table = report.getTable(readFeature.getName());
                table.addColumn(readFeature.getKey());

                for (SAMReadGroupRecord rg : readGroups) {
                    table.addColumn(rg.getId());
                }

                featuresToExtract.add(featureClass.getConstructor(RFAArgumentCollection.class).newInstance(collection));
            }
        } catch ( Exception e ) {
            throw new StingException("A read feature instantiation error occurred during initialization",e);
        }
    }

    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
        return (read.getReadPairedFlag() && read.getFirstOfPairFlag() && ! read.getMateUnmappedFlag());
    }

    @Override
    public Integer map(ReferenceContext referenceContext, GATKSAMRecord samRecord, ReadMetaDataTracker readMetaDataTracker) {
        for ( ReadFeature feature : featuresToExtract ) {
            if ( feature.isDefinedFor(samRecord) ) {
                GATKReportTable table = report.getTable(feature.getName());
                Object value = feature.getFeature(samRecord);
                String rgid = samRecord.getReadGroup().getReadGroupId();
                table.set(value, feature.getKey(), value);
                final Object currentVal = table.get(value,rgid);
                if (  currentVal == null )
                    table.set(value,rgid,1);
                else
                    table.set(value,rgid,(Integer)currentVal+1);
            }
        }

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
        // Write report.
        report.print(out);
    }

    public Set<Class<? extends ReadFeature>> getFeatureAggregators(List<String> requestedFeatures) {
        HashSet<Class<? extends ReadFeature>> newFeatureSet = new HashSet<Class<? extends ReadFeature>>();
        List<Class<? extends ReadFeature>> availableFeatures = new PluginManager<ReadFeature>(ReadFeature.class).getPlugins();

        if ( collection.inputFeatures == null ) {
            newFeatureSet.addAll(availableFeatures);
            return newFeatureSet;
        }


        Map<String,Class<? extends ReadFeature>> classNameToClass = new HashMap<String,Class<? extends ReadFeature>>(collection.inputFeatures.size());
        for ( Class<? extends ReadFeature> clazz : availableFeatures ) {
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
}
