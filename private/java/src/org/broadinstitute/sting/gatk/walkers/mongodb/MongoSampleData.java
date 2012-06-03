package org.broadinstitute.sting.gatk.walkers.mongodb;

import com.mongodb.*;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 5/11/12
 * Time: 3:23 PM
 *
 * Handles the MongoDB samples collection data processing
 */
public class MongoSampleData {
    private final String contig;
    private final Integer start;
    private final Integer stop;

    private final List<Allele> siteAlleles;
    private final String sourceROD;

    private final Genotype genotype;

    protected Genotype getGenotype() {
        return genotype;
    }

    protected static void ensurePrimaryKey() {
        MongoDB.getSamplesCollection().ensureIndex(new BasicDBObject("sample", 1).append("block", 1), new BasicDBObject("unique", 1));
    }

    protected MongoSampleData(VariantContext vc, Genotype pGenotype, String pSourceROD) {
        contig = vc.getChr();
        start = vc.getStart();
        stop = vc.getEnd();

        siteAlleles = vc.getAlleles();
        sourceROD = pSourceROD;

        genotype = pGenotype;
    }

    /**
     * Private constructor for internal use
     *
     * @param contig
     * @param start
     * @param stop
     * @param siteAlleles
     * @param sourceROD
     * @param genotype
     */
    private MongoSampleData(String contig, Integer start, Integer stop, List<Allele> siteAlleles, String sourceROD, Genotype genotype) {
        this.contig = contig;
        this.start = start;
        this.stop = stop;
        this.siteAlleles = siteAlleles;
        this.sourceROD = sourceROD;
        this.genotype = genotype;
    }

    protected static void insertIntoMongoDb(MongoBlockKey block_id, String sampleName, List<MongoSampleData> samplesPerBlock) {
        DBCollection collection = MongoDB.getSamplesCollection();

        BasicDBObject blockDoc = new BasicDBObject();
        blockDoc.put("block", block_id.getKey());
        blockDoc.put("sample", sampleName);

        BasicDBList sitesList = new BasicDBList();
        for (MongoSampleData sample: samplesPerBlock) {
            BasicDBObject siteDoc = new BasicDBObject();
            siteDoc.put("contig", sample.contig);
            siteDoc.put("start", sample.start);
            siteDoc.put("stop", sample.stop);

            Integer alleleIndex = 0;
            BasicDBObject allelesDoc = new BasicDBObject();
            for (Allele allele : sample.siteAlleles)
            {
                String index = alleleIndex.toString();
                allelesDoc.put(index, allele.toString());
                alleleIndex++;
            }
            siteDoc.put("alleles", allelesDoc);

            siteDoc.put("sourceROD", sample.sourceROD);

            BasicDBObject genotypesDoc = new BasicDBObject();
            Integer genotypeAlleleIndex = 0;
            BasicDBObject genotypeAllelesDoc = new BasicDBObject();
            for (Allele allele : sample.genotype.getAlleles())
            {
                String index = genotypeAlleleIndex.toString();
                genotypeAllelesDoc.put(index, allele.toString());
                genotypeAlleleIndex++;
            }
            genotypesDoc.put("alleles", genotypeAllelesDoc);

            BasicDBObject attributesList = new BasicDBObject();
            // TODO -- add code to handle GQ, AD, PL
            for (Map.Entry<String, Object> attribute : sample.genotype.getExtendedAttributes().entrySet())
            {
                String key = attribute.getKey();
                Object value = attribute.getValue();
                attributesList.put(key, value);
            }
            genotypesDoc.put("attributes", attributesList);

            genotypesDoc.put("error", sample.genotype.getLog10PError());
            siteDoc.put("genotype", genotypesDoc);

            sitesList.add(siteDoc);
        }

        blockDoc.put("sites", sitesList);
        collection.insert(blockDoc);
    }

    protected static Map<Integer, List<MongoSampleData>> retrieveFromMongo(MongoBlockKey block_id, String sample) {
        Map<Integer, List<MongoSampleData>> returnMap = new HashMap<Integer, List<MongoSampleData>>();

        BasicDBObject query = new BasicDBObject();
        query.put("sample", sample);
        query.put("block", block_id.getKey());

        DBCursor cursor = MongoDB.getSamplesCollection().find(query);
        while(cursor.hasNext()) {
            DBObject blockResult = cursor.next();

            for (BasicDBObject siteDoc : (List<BasicDBObject>)blockResult.get("sites")) {
                String pContig = (String)siteDoc.get("contig");
                Integer pStart = (Integer)siteDoc.get("start");
                Integer pStop = (Integer)siteDoc.get("stop");

                ArrayList<Allele> pAlleles = new ArrayList<Allele>();
                BasicDBObject allelesInDb = (BasicDBObject)siteDoc.get("alleles");
                for (Object alleleInDb : allelesInDb.values()) {
                    String rawAllele = (String)alleleInDb;
                    boolean isRef = rawAllele.contains("*");
                    String allele = rawAllele.replace("*", "");
                    pAlleles.add(Allele.create(allele, isRef));
                }

                String pSourceROD = (String)siteDoc.get("sourceROD");

                BasicDBObject genotypeDoc = (BasicDBObject)siteDoc.get("genotype");

                ArrayList<Allele> genotypeAlleles = new ArrayList<Allele>();
                BasicDBObject genotypeAllelesInDb = (BasicDBObject)genotypeDoc.get("alleles");
                for (Object alleleInDb : genotypeAllelesInDb.values()) {
                    String rawAllele = (String)alleleInDb;
                    boolean isRef = rawAllele.contains("*");
                    String allele = rawAllele.replace("*", "");
                    genotypeAlleles.add(Allele.create(allele, isRef));
                }
                Double genotypeError = (Double)genotypeDoc.get("error");

                Map<String, Object> genotypeAttributes = new HashMap<String, Object>();
                BasicDBObject genotypeAttrsInDb = (BasicDBObject)genotypeDoc.get("attributes");
                for (String key : genotypeAttrsInDb.keySet()) {
                    Object value = genotypeAttrsInDb.get(key);
                    genotypeAttributes.put(key, value);
                }

                Genotype pGenotype = new GenotypeBuilder(sample, genotypeAlleles)
                        .log10PError(genotypeError).attributes(genotypeAttributes).make();

                if (!returnMap.containsKey(pStart)) {
                    returnMap.put(pStart, new ArrayList<MongoSampleData>());
                }
                returnMap.get(pStart).add(new MongoSampleData(pContig, pStart, pStop, pAlleles, pSourceROD, pGenotype));
            }
        }

        return returnMap;
    }

    protected boolean matches(MongoSiteData site) {
        return contig.equals(site.getContig())
                && start.equals(site.getStart())
                && siteAlleles.equals(site.getAlleles())
                && sourceROD.equals(site.getSourceROD());
    }
}
