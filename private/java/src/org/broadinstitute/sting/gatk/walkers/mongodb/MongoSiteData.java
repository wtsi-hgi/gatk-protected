package org.broadinstitute.sting.gatk.walkers.mongodb;

import com.mongodb.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.util.*;

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

    protected String getContig() {
        return contig;
    }

    protected Integer getStart() {
        return start;
    }

    protected List<Allele> getAlleles() {
        return alleles;
    }

    protected String getSourceROD() {
        return sourceROD;
    }

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
        blockDoc.put("block", block_id.getKey());

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

    protected static Map<Integer, List<MongoSiteData>> retrieveFromMongo(MongoBlockKey block_id) {
        Map<Integer, List<MongoSiteData>> returnMap = new HashMap<Integer, List<MongoSiteData>>();

        BasicDBObject query = new BasicDBObject();
        query.put("block", block_id.getKey());

        DBCursor cursor = MongoDB.getSitesCollection().find(query);
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
                String pId = (String)siteDoc.get("id");
                Double pError = (Double)siteDoc.get("error");
                String pSource = (String)siteDoc.get("source");
                VariantContext.Type pType = VariantContext.Type.valueOf((String)siteDoc.get("type"));

                Set<String> pFilters = new HashSet<String>();
                BasicDBList filtersInDb = (BasicDBList)siteDoc.get("filters");
                if (filtersInDb != null) {
                    for (Object filterInDb : filtersInDb) {
                        pFilters.add((String)filterInDb);
                    }
                }

                Map<String, Object> pAttributes = new HashMap<String, Object>();
                BasicDBObject attrsInDb = (BasicDBObject)siteDoc.get("attributes");
                for (String key : attrsInDb.keySet()) {
                    Object value = attrsInDb.get(key);
                    pAttributes.put(key, value);
                }

                if (!returnMap.containsKey(pStart)) {
                    returnMap.put(pStart, new ArrayList<MongoSiteData>());
                }
                returnMap.get(pStart).add(new MongoSiteData(pContig, pStart, pStop, pAlleles, pSourceROD, pId, pError, pSource, pType, pFilters, pAttributes));
            }
        }

        return returnMap;
    }

    protected VariantContextBuilder builder(ReferenceContext ref) {
        VariantContextBuilder builder = new VariantContextBuilder(source, contig, start, stop, alleles);

        builder.id(id);
        builder.log10PError(error);
        builder.attributes(attributes);
        builder.filters(filters);

        return builder;
    }
}

