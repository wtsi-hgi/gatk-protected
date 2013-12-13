/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.mongodb;

import com.mongodb.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;

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

