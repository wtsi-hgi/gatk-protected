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
import org.broadinstitute.variant.variantcontext.*;

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

            // populate inline attributes
            if (sample.genotype.hasPL()) {
                attributesList.put("PL", sample.genotype.getPL());
            }
            if (sample.genotype.hasDP()) {
                attributesList.put("DP", sample.genotype.getDP());
            }
            if (sample.genotype.hasAD()) {
                attributesList.put("AD", sample.genotype.getAD());
            }
            if (sample.genotype.hasGQ()) {
                attributesList.put("GQ", sample.genotype.getGQ());
            }

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

                Map<String, Object> genotypeExtendedAttributes = new HashMap<String, Object>();

                // defaults from GenotypeBuilder.  TODO: change to constants?
                int GQ = -1;
                int DP = -1;
                int[] AD = null;
                int[] PL = null;

                BasicDBObject genotypeAttrsInDb = (BasicDBObject)genotypeDoc.get("attributes");
                for (String key : genotypeAttrsInDb.keySet()) {
                    Object value = genotypeAttrsInDb.get(key);

                    if (key.equals("GQ")) {
                        GQ = (Integer)value;
                    }
                    else if (key.equals("DP")) {
                        DP = (Integer)value;
                    }
                    else if (key.equals("AD")) {
                        BasicDBList ADInDb = (BasicDBList)value;
                        AD = new int[ADInDb.size()];
                        for (int counter = 0; counter < ADInDb.size(); counter++) {
                            AD[counter] = (Integer)ADInDb.get(counter);
                        }
                    }
                    else if (key.equals("PL")) {
                        BasicDBList PLInDb = (BasicDBList)value;
                        PL = new int[PLInDb.size()];
                        for (int counter = 0; counter < PLInDb.size(); counter++) {
                            PL[counter] = (Integer)PLInDb.get(counter);
                        }
                    }
                    else {
                        genotypeExtendedAttributes.put(key, value);
                    }
                }

                Genotype pGenotype = new GenotypeBuilder(sample, genotypeAlleles)
                        .log10PError(genotypeError).GQ(GQ).DP(DP).AD(AD).PL(PL)
                        .attributes(genotypeExtendedAttributes)
                        .make();

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
