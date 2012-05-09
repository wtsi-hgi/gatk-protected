package org.broadinstitute.sting.gatk.walkers.mongodb;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 3/30/12
 * Time: 4:47 PM
 * To change this template use File | Settings | File Templates.
 */

import com.mongodb.BasicDBObject;
import com.mongodb.DBCollection;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Inserts all of the RODs in the input data set. Data is inserted using toMongoDB().
 */
public class InsertRODsWalker extends RodWalker<Integer, Integer> {
    @Input(fullName="input", shortName = "input", doc="The input ROD which should be inserted into the DB.", required=true)
    public RodBinding<Feature> input;

    @Output
    PrintStream out;


    private String RODFileName;

    @Override
    public void initialize() {
        DBCollection mongoAttributes = MongoDB.getAttributesCollection();
        DBCollection mongoSamples = MongoDB.getSamplesCollection();

        RODFileName = input.getSource();
        int lastSep = RODFileName.lastIndexOf(File.separator);
        RODFileName = RODFileName.substring(lastSep + 1);

        // set up primary keys

        mongoAttributes.ensureIndex(new BasicDBObject("location", 1).append("sourceROD", 1).append("alleles", 1), new BasicDBObject("unique", 1));
        mongoSamples.ensureIndex(new BasicDBObject("location", 1).append("sourceROD", 1).append("alleles", 1).append("sample", 1), new BasicDBObject("unique", 1));
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }

    /**
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        DBCollection mongoAttributes = MongoDB.getAttributesCollection();
        DBCollection mongoSamples = MongoDB.getSamplesCollection();

        for ( Feature feature : tracker.getValues(Feature.class, context.getLocation()) ) {
            if ( feature instanceof VariantContext ) {
                VariantContext vc = (VariantContext) feature;

                Pair<BasicDBObject,List<BasicDBObject>> mongoCollections = toMongoDB(vc, RODFileName);
                mongoAttributes.insert(mongoCollections.first);
                for (BasicDBObject sampleForMongo : mongoCollections.second) {
                    mongoSamples.insert(sampleForMongo);
                }
            }
        }

        return 1;
    }

    /**
     * Generate a Mongo DB attributes collection element and a set of samples collection elements
     * @param sourceROD
     * @return
     */
    public Pair<BasicDBObject,List<BasicDBObject>> toMongoDB(VariantContext vc, String sourceROD) {
        // fields common to both attributes and samples collections
        BasicDBObject siteDoc = new BasicDBObject();

        String contig = vc.getChr();
        int start = vc.getStart();
        int stop = vc.getEnd();

        siteDoc.put("location", contig + ":" + (start - stop == 0 ? start : start + "-" + stop));
        siteDoc.put("contig", contig);
        siteDoc.put("start", start);
        siteDoc.put("stop", stop);
        siteDoc.put("id", vc.getID());
        siteDoc.put("error", vc.getLog10PError());
        siteDoc.put("source", vc.getSource());
        siteDoc.put("sourceROD", sourceROD);
        siteDoc.put("type", vc.getType().toString());

        Integer alleleIndex = 0;
        BasicDBObject allelesDoc = new BasicDBObject();
        for (Allele allele : vc.getAlleles())
        {
            String index = alleleIndex.toString();
            allelesDoc.put(index, allele.toString());
            alleleIndex++;
        }
        siteDoc.put("alleles", allelesDoc);

        Integer filterIndex = 0;
        BasicDBObject filtersDoc = new BasicDBObject();
        for (String filter : vc.getFilters())
        {
            String index = filterIndex.toString();
            filtersDoc.put(index, filter.toString());
            filterIndex++;
        }
        if (filterIndex > 0) {
            siteDoc.put("filters", filtersDoc);
        }

        // attributes collection

        BasicDBObject attributesDoc = new BasicDBObject(siteDoc);
        List<BasicDBObject> attributeKVPs = new ArrayList<BasicDBObject>();
        for (Map.Entry<String, Object> attribute : vc.getAttributes().entrySet() )
        {
            String key = attribute.getKey();
            Object value = attribute.getValue();
            BasicDBObject attributeKVP = new BasicDBObject();
            attributeKVP.put("key", key);
            attributeKVP.put("value", value);
            attributeKVPs.add(attributeKVP);
        }
        attributesDoc.put("attributes", attributeKVPs);

        // samples collection

        List<BasicDBObject> samplesDocs = new ArrayList<BasicDBObject>();
        for (Genotype genotype : vc.getGenotypes()) {
            BasicDBObject sampleDoc = new BasicDBObject(siteDoc);
            sampleDoc.put("sample", genotype.getSampleName());

            BasicDBObject genotypesDoc = new BasicDBObject();
            Integer genotypeAlleleIndex = 0;
            BasicDBObject genotypeAllelesDoc = new BasicDBObject();
            for (Allele allele : genotype.getAlleles())
            {
                String index = genotypeAlleleIndex.toString();
                genotypeAllelesDoc.put(index, allele.toString());
                genotypeAlleleIndex++;
            }
            genotypesDoc.put("alleles", genotypeAllelesDoc);

            List<BasicDBObject> genotypesAttributesDocs = new ArrayList<BasicDBObject>();
            for (Map.Entry<String, Object> attribute : genotype.getAttributes().entrySet() )
            {
                String key = attribute.getKey();
                Object value = attribute.getValue();
                BasicDBObject genotypesAttributesDoc = new BasicDBObject();
                genotypesAttributesDoc.put("key", key);
                genotypesAttributesDoc.put("value", value);
                genotypesAttributesDocs.add(genotypesAttributesDoc);
            }
            genotypesDoc.put("attributes", genotypesAttributesDocs);
            genotypesDoc.put("error", genotype.getLog10PError());

            sampleDoc.put("genotype", genotypesDoc);

            samplesDocs.add(sampleDoc);
        }

        return new Pair<BasicDBObject,List<BasicDBObject>>(attributesDoc, samplesDocs);
    }

    /**
     * Increment the number of rods processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of rods processed.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        MongoDB.close();
    }
}