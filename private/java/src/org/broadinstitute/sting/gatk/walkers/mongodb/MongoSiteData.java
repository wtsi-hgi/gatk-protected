package org.broadinstitute.sting.gatk.walkers.mongodb;

import com.mongodb.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.List;
import java.util.Map;
import java.util.Set;

/**                   
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 5/11/12
 * Time: 3:23 PM
 *
 * Handles the MongoDB sites collection data processing
 */
public class MongoSiteData {
    private final String contig;
    private final Integer start;
    private final Integer stop;

    private final List<Allele> alleles;
    private final String sourceROD;

    private final String id;
    private final Double error;
    private final String source;
    private final VariantContext.Type type;
    private final Set<String> filters;
    private final Map<String, Object> attributes;

    protected static void ensurePrimaryKey() {
        MongoDB.getSitesCollection().ensureIndex(new BasicDBObject("block", 1), new BasicDBObject("unique", 1));
    }

    protected MongoSiteData(VariantContext vc, String pSourceROD) {
        contig = vc.getChr();
        start = vc.getStart();
        stop = vc.getEnd();

        alleles = vc.getAlleles();
        sourceROD = pSourceROD;

        id = vc.getID();
        error = vc.getLog10PError();
        source = vc.getSource();
        type = vc.getType();
        filters = vc.getFilters();
        attributes = vc.getAttributes();
    }

    /**
     * Private constructor for internal use
     *
     * @param contig
     * @param start
     * @param stop
     * @param alleles
     * @param sourceROD
     * @param id
     * @param error
     * @param source
     * @param type
     * @param filters
     * @param attributes
     */
    private MongoSiteData(String contig, Integer start, Integer stop, List<Allele> alleles, String sourceROD, String id, Double error, String source, VariantContext.Type type, Set<String> filters, Map<String, Object> attributes) {
        this.contig = contig;
        this.start = start;
        this.stop = stop;
        this.alleles = alleles;
        this.sourceROD = sourceROD;
        this.id = id;
        this.error = error;
        this.source = source;
        this.type = type;
        this.filters = filters;
        this.attributes = attributes;
    }

    // Important: does NOT insert/update if the DB already has data for this block
    // The implication here is that new samples are very likely to include sites which are not in the sites collection
    // Let's call that a TODO.

    protected static void insertIntoMongoDb(MongoBlockKey block_id, List<MongoSiteData> sites) {
        DBCollection collection = MongoDB.getSitesCollection();

        BasicDBObject blockDoc = new BasicDBObject();
        blockDoc.put("block", block_id.key);

        BasicDBList sitesList = new BasicDBList();
        for (MongoSiteData site: sites) {
            BasicDBObject siteDoc = new BasicDBObject();
            siteDoc.put("contig", site.contig);
            siteDoc.put("start", site.start);
            siteDoc.put("stop", site.stop);

            Integer alleleIndex = 0;
            BasicDBObject allelesDoc = new BasicDBObject();
            for (Allele allele : site.alleles)
            {
                String index = alleleIndex.toString();
                allelesDoc.put(index, allele.toString());
                alleleIndex++;
            }
            siteDoc.put("alleles", allelesDoc);

            siteDoc.put("sourceROD", site.sourceROD);
            siteDoc.put("id", site.id);
            siteDoc.put("error", site.error);
            siteDoc.put("source", site.source);
            siteDoc.put("type", site.type.toString());

            Integer filterIndex = 0;
            BasicDBList filtersList = new BasicDBList();
            for (String filter : site.filters) {
                String index = filterIndex.toString();
                filtersList.put(index, filter);
                filterIndex++;
            }
            if (filterIndex > 0) {
                siteDoc.put("filters", filtersList);
            }

            BasicDBObject attributes = new BasicDBObject();
            for (Map.Entry<String, Object> attribute : site.attributes.entrySet() )
            {
                String key = attribute.getKey();
                Object value = attribute.getValue();
                attributes.put(key, value);
            }
            siteDoc.put("attributes", attributes);
            sitesList.add(siteDoc);
        }

        blockDoc.put("sites", sitesList);
        collection.insert(blockDoc);
    }
}
