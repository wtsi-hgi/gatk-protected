package org.broadinstitute.sting.gatk.walkers.mongodb;

import com.mongodb.*;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.List;
import java.util.Map;

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

    private final List<Allele> alleles;
    private final String sourceROD;

    private final Genotype genotype;

    protected static void ensurePrimaryKey() {
        MongoDB.getSamplesCollection().ensureIndex(new BasicDBObject("sample", 1).append("block", 1), new BasicDBObject("unique", 1));
    }

    protected MongoSampleData(VariantContext vc, Genotype pGenotype, String pSourceROD) {
        contig = vc.getChr();
        start = vc.getStart();
        stop = vc.getEnd();

        alleles = vc.getAlleles();
        sourceROD = pSourceROD;

        genotype = pGenotype;
    }

    /**
     * Private constructor for internal use
     *
     * @param contig
     * @param start
     * @param stop
     * @param alleles
     * @param sourceROD
     * @param genotype
     */
    private MongoSampleData(String contig, Integer start, Integer stop, List<Allele> alleles, String sourceROD, Genotype genotype) {
        this.contig = contig;
        this.start = start;
        this.stop = stop;
        this.alleles = alleles;
        this.sourceROD = sourceROD;
        this.genotype = genotype;
    }

    protected static void insertIntoMongoDb(MongoBlockKey block_id, String sampleName, List<MongoSampleData> samplesPerBlock) {
        DBCollection collection = MongoDB.getSamplesCollection();

        BasicDBObject blockDoc = new BasicDBObject();
        blockDoc.put("block", block_id.key);
        blockDoc.put("sample", sampleName);

        BasicDBList sitesList = new BasicDBList();
        for (MongoSampleData sample: samplesPerBlock) {
            BasicDBObject siteDoc = new BasicDBObject();
            siteDoc.put("contig", sample.contig);
            siteDoc.put("start", sample.start);
            siteDoc.put("stop", sample.stop);

            Integer alleleIndex = 0;
            BasicDBObject allelesDoc = new BasicDBObject();
            for (Allele allele : sample.alleles)
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
            for (Map.Entry<String, Object> attribute : sample.genotype.getAttributes().entrySet())
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
}
