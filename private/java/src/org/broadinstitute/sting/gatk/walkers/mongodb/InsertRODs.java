package org.broadinstitute.sting.gatk.walkers.mongodb;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 3/30/12
 * Time: 4:47 PM
 */

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Inserts all of the RODs in the input data set into MongoDB.
 *
 * Sites are grouped into blocks for storage into MongoDB.
 * Sites and samples data are aggregated by block and inserted once a block is complete.
 *
 * TODO: In the current implementation, a new sample at a location not currently referenced in the site table
 * will NOT insert the location into the sites table (collection in MongoDB-speak)
 */
@PartitionBy(PartitionType.NONE)
public class InsertRODs extends RodWalker<Integer, Integer> {
    @Input(fullName="input", shortName = "input", doc="The input ROD which should be inserted into the DB.", required=true)
    public RodBinding<Feature> input;

    private String RODFileName;
    private MongoBlockKey activeBlockKey = null;

    private List<MongoSiteData> activeBlockSiteData = new ArrayList<MongoSiteData>();
    private Map<String,List<MongoSampleData>> activeBlockSampleMap = new HashMap<String,List<MongoSampleData>>();

    @Override
    public void initialize() {
        RODFileName = input.getSource();
        int lastSep = RODFileName.lastIndexOf(File.separator);
        RODFileName = RODFileName.substring(lastSep + 1);

        MongoSiteData.ensurePrimaryKey();
        MongoSampleData.ensurePrimaryKey();
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    @Override
    public Integer reduceInit() { return 0; }

    /**
     * Insert the site and sample data from the ROD
     * TODO? Currently only operates on VariantContext
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        for ( Feature feature : tracker.getValues(Feature.class, context.getLocation()) ) {
            if ( feature instanceof VariantContext ) {
                VariantContext vc = (VariantContext) feature;
                addToActiveBlock(vc, RODFileName);
            }
        }

        return 1;
    }

    /**
     * Adds the VariantContext to the active block and inserts the previous block if complete
     *
     * @param vc            the VariantContext
     * @param sourceROD     the source filename
     */
    public void addToActiveBlock(VariantContext vc, String sourceROD) {
        MongoBlockKey block_id = new MongoBlockKey(vc);

        // Blocks must be inserted in order.  If the previous block is complete, insert it.

        if (!block_id.equals(activeBlockKey)) {
            if (activeBlockKey != null) {
                insertActiveBlockIntoMongo();
            }
            activeBlockKey = block_id;
            activeBlockSiteData = new ArrayList<MongoSiteData>();
            activeBlockSampleMap = new HashMap<String, List<MongoSampleData>>();
        }

        activeBlockSiteData.add(new MongoSiteData(vc, sourceROD));

        for (Genotype genotype : vc.getGenotypes()) {
            String sampleName = genotype.getSampleName();
            if (!activeBlockSampleMap.containsKey(sampleName)) {
                activeBlockSampleMap.put(sampleName, new ArrayList<MongoSampleData>());
            }

            activeBlockSampleMap.get(sampleName).add(new MongoSampleData(vc, genotype, sourceROD));
        }
    }

    /**
     * Increment the number of rods processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of rods processed.
     */
    @Override
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    /**
     * Insert the final block into MongoDB and close the connection
     *
     * @param result    the number of rods processed
     */
    @Override
    public void onTraversalDone(Integer result) {
        insertActiveBlockIntoMongo();
        MongoDB.close();
    }

    /**
     * Insert the current block into MongoDB
     */
    private void insertActiveBlockIntoMongo() {
        MongoSiteData.insertIntoMongoDb(activeBlockKey, activeBlockSiteData);
        for (String sampleName : activeBlockSampleMap.keySet()) {
            MongoSampleData.insertIntoMongoDb(activeBlockKey, sampleName, activeBlockSampleMap.get(sampleName));
        }
    }
}