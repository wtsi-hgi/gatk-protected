package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import com.google.java.contract.Requires;
import com.mongodb.*;
import org.apache.log4j.Logger;
import org.apache.log4j.Priority;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.errors.InvalidRecordsRemove;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFStandardHeaderLines;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;

import java.util.*;

public class NA12878KnowledgeBase {
    private final static Logger logger = Logger.getLogger(NA12878KnowledgeBase.class);

    protected DBCollection sites;
    protected DBCollection callSets;
    protected DBCollection consensusSites;
    protected GenomeLocParser parser;

    private final MongoDBManager.Locator dblocator;

    public NA12878KnowledgeBase(final GenomeLocParser parser, final NA12878DBArgumentCollection args) {
        this.parser = parser;
        this.dblocator = args.getLocator();

        MongoDBManager.DBWrapper dbWrapper = MongoDBManager.getDB(this.dblocator);

        sites = dbWrapper.getSites();
        sites.setObjectClass(MongoVariantContext.class);
        sites.ensureIndex(sitesOrder());
        sites.ensureIndex(essentialIndex());

        callSets = dbWrapper.getCallsets();
        callSets.setObjectClass(CallSet.class);

        consensusSites = dbWrapper.getConsensus();
        consensusSites.setObjectClass(MongoVariantContext.class);
        consensusSites.ensureIndex(sitesOrder());
        consensusSites.ensureIndex(essentialIndex());

        dbWrapper.getMongo().setWriteConcern(WriteConcern.SAFE);

        if ( args.resetDB ) reset();
    }

    /**
     * Close down the connections for this KnowledgeBase
     */
    public void close() {
        MongoDBManager.getDB(dblocator).close();
    }

    public void delete() {
        MongoDBManager.getDB(dblocator).delete();
    }

    protected void printStatus() {
        printCollectionStatus("sites", sites);
        printCollectionStatus("callsets", callSets);
        printCollectionStatus("consensus", consensusSites);
    }

    private void printCollectionStatus(final String name, final DBCollection collection) {
        logger.info("Number of " + name + ":" + collection.getCount());

        for ( final DBObject object : collection.find() ) {
            logger.info("        " + object);
        }
    }

    public void reset() {
        logger.warn("Resetting all documents in " + this);
        sites.drop();
        callSets.drop();
        clearConsensus();
    }

    public void clearConsensus() {
        consensusSites.drop();
    }

    // ---------------------------------------------------------------------
    // working with CallSets
    // ---------------------------------------------------------------------

    public WriteResult addCallset(final CallSet callSet) {
        return callSets.insert(callSet);
    }

    public List<CallSet> getCallSets() {
        final List<CallSet> callsets = new LinkedList<CallSet>();

        for ( final DBObject obj : callSets.find() ) {
            callsets.add((CallSet)obj);
        }

        return callsets;
    }

    public CallSet getCallSet(final String name) {
        final BasicDBObject query = new BasicDBObject("Name", name);
        return (CallSet)callSets.findOne(query);
    }

    // ---------------------------------------------------------------------
    // Working with sites
    // ---------------------------------------------------------------------

    public WriteResult addCall(final MongoVariantContext mvc) {
        return sites.insert(mvc);
    }

    public List<WriteResult> removeCall(final MongoVariantContext mvc){
        Set<String> matchKeys = new HashSet<String>(mvc.keySet());
        matchKeys.remove("Date");
        matchKeys.remove("_id");
        DBObject matchObject = new BasicDBObject();
        for(String key: matchKeys){
            matchObject.put(key, mvc.get(key));
        }
        DBCursor cursor = sites.find(matchObject);
        List<WriteResult> results = new ArrayList<WriteResult>(cursor.size());
        for(DBObject next: cursor){
            results.add(sites.remove(next));
        }
        return results;
    }

    public void addCalls(final Collection<MongoVariantContext> mvcs) {
        for ( final MongoVariantContext mvc : mvcs )
            addCall(mvc);
    }

    private DBObject essentialIndex() {
        final DBObject sortOrder = new BasicDBObject();
        sortOrder.put("Chr", 1);
        sortOrder.put("Start", 1);
        return sortOrder;
    }

    protected DBCollection getSitesCollection() {
        return sites;
    }

    public DBObject sitesOrder() {
        final DBObject sortOrder = new BasicDBObject();
        sortOrder.put("Chr", 1);
        sortOrder.put("Start", 1);
        sortOrder.put("Stop", 1);
        sortOrder.put("Ref", 1);
        sortOrder.put("Alt", 1);
        return sortOrder;
    }

    /**
     * Get all of the calls in the db
     * @return
     */
    public SiteIterator<MongoVariantContext> getCalls() {
        return getCalls(new SiteSelector(parser));
    }

    /**
     * Get the sites, if any, that overlap the interval specified as the loc
     *
     * @param criteria
     * @return
     */
    public SiteIterator<MongoVariantContext> getCalls(final SiteSelector criteria) {
        final DBObject sortOrder = sitesOrder();
        final DBCursor cursor = sites.find(criteria.toQuery()).sort(sortOrder);
        final SiteIterator<MongoVariantContext> it = new SiteIterator<MongoVariantContext>(parser, cursor);
        it.setErrorHandler(new InvalidRecordsRemove<MongoVariantContext>(sites));
        return it;
    }

    /**
     * Get the consensus sites, if any, that overlap the interval specified as the loc
     *
     * @param criteria a SiteSelector to select the consensus sites we are interested in
     * @return
     */
    public SiteIterator<MongoVariantContext> getConsensusSites(final SiteSelector criteria) {
        final DBCursor cursor = consensusSites.find(criteria.toQuery()).sort(sitesOrder());
        return new SiteIterator<MongoVariantContext>(parser, cursor);
    }

    public ConsensusSummarizer updateConsensus(final SiteSelector selector) {
        return updateConsensus(selector, Priority.DEBUG);
    }

    public ConsensusSummarizer updateConsensus(final SiteSelector selector, final Priority logPriority) {
        final ConsensusSummarizer summary = new ConsensusSummarizer();

        final SiteIterator<MongoVariantContext> siteIterator = getCalls(selector);
        while ( siteIterator.hasNext() ) {
            final Collection<MongoVariantContext> callsAtSite = siteIterator.getNextEquivalents();
            final MongoVariantContext consensus = makeConsensus(callsAtSite);
            addConsensus(consensus);
            logger.log(logPriority, "Updating consensus at site " + consensus);
            summary.add(consensus);
        }

        return summary;
    }

    private MongoVariantContext makeConsensus(final Collection<MongoVariantContext> allSupportingCalls) {
        final List<String> supportingCallSets = new LinkedList<String>();
        final Collection<MongoVariantContext> callsForConsensus = selectCallsForConsensus(allSupportingCalls);
        final boolean isReviewed = isReviewed(callsForConsensus);
        final VariantContext first = callsForConsensus.iterator().next().getVariantContext();

        final VariantContextBuilder builder = new VariantContextBuilder();
        builder.chr(first.getChr()).start(first.getStart()).stop(first.getEnd());

        final Set<Allele> alleles = new LinkedHashSet<Allele>();
        for ( final MongoVariantContext vc : allSupportingCalls ) {
            alleles.addAll(vc.getVariantContext().getAlleles());
            supportingCallSets.addAll(vc.getSupportingCallSets());
        }
        builder.alleles(alleles);

        final TruthStatus type = determineTruth(callsForConsensus);
        final PolymorphicStatus status = determinePolymorphicStatus(callsForConsensus);
        final Genotype gt = consensusGT(type, status, new LinkedList<Allele>(alleles), callsForConsensus);

        return MongoVariantContext.create(supportingCallSets, builder.make(), type, new Date(), gt, isReviewed);
    }

    /**
     * Is at least one call in individualCalls a reviewed call?
     * @param individualCalls a collection of calls to consider
     * @return true if at least one call in individualCalls is a review, false otherwise
     */
    @Requires("individualCalls != null")
    private boolean isReviewed(final Collection<MongoVariantContext> individualCalls) {
        for ( final MongoVariantContext vc : individualCalls )
            if ( vc.isReviewed() )
                return true;
        return false;
    }

    private Collection<MongoVariantContext> selectCallsForConsensus(final Collection<MongoVariantContext> individualCalls) {
        final List<MongoVariantContext> reviewed = new LinkedList<MongoVariantContext>();

        for ( final MongoVariantContext vc : individualCalls )
            if ( vc.isReviewed() ) reviewed.add(vc);
        return reviewed.isEmpty() ? individualCalls : reviewed;
    }

    /**
     * Create a consensus genotype appropriate for a site backed by individualCalls with given
     * truthStatus and polyStatus
     *
     * @param truthStatus the determined truth status for this site
     * @param polyStatus the determined polymorphic status of this site // TODO -- why is this necessary?
     * @param alleles the list of alleles segregating at this site
     * @param individualCalls the individual call sets we should use to make the consensus gt
     * @return a Genotype appropriate for this consensus site
     */
    protected Genotype consensusGT(final TruthStatus truthStatus,
                                   final PolymorphicStatus polyStatus,
                                   final List<Allele> alleles,
                                   final Collection<MongoVariantContext> individualCalls) {
        if ( ! truthStatus.isTruePositive() ) {
            return MongoGenotype.NO_CALL;
        } else if ( polyStatus.isDiscordant() || polyStatus.isUnknown() ) {
            return MongoGenotype.NO_CALL;
        } else {
            Genotype g = MongoGenotype.NO_CALL;

            // we are a TP, we need to compute the consensus genotype
            for ( final MongoVariantContext mvc : individualCalls ) {
                final Genotype mvcG = mvc.getGt().toGenotype(alleles);
                if ( g.isNoCall() )
                    g = mvcG;
                else if ( mvcG.isNoCall() )
                    ; // keep g
                else if ( g.isMixed() || ! g.isAvailable() )
                    throw new IllegalStateException("Unexpected genotype in mongo db " + g + " at " + individualCalls);
                else if ( g.getType() != mvcG.getType() )
                    return MongoGenotype.createDiscordant(mvcG);
                else
                    ; // TODO -- should try to capture more DP and GQ
            }

            return g;
        }
    }

    private TruthStatus determineTruth(final Collection<MongoVariantContext> individualCalls) {
        final boolean hasReview = isReviewed(individualCalls);
        TruthStatus type = TruthStatus.UNKNOWN;
        for ( final MongoVariantContext vc : individualCalls ) {
            if ( ! hasReview || vc.isReviewed() )
                // if we have some reviews, only include those, otherwise use everything
                type = type.makeConsensus(vc.getType());
        }

        return type;
    }

    private PolymorphicStatus determinePolymorphicStatus(final Collection<MongoVariantContext> individualCalls) {
        final boolean hasReview = isReviewed(individualCalls);
        PolymorphicStatus status = PolymorphicStatus.UNKNOWN;
        for ( final MongoVariantContext vc : individualCalls ) {
            if ( ! hasReview || vc.isReviewed() )
                // if we have some reviews, only include those, otherwise use everything
                status = status.makeConsensus(vc.getPolymorphicStatus());
        }

        return status;
    }

    private WriteResult addConsensus(final MongoVariantContext site) {
        return consensusSites.insert(site);
    }

    /**
     * Convenience function to write all reviewed MongoVariantContexts as VariantContext
     * from the sites collection to the VariantContextWriter outputStream
     *
     * @param writer where the VariantContexts should be written
     * @param selector the root SiteSelector that lets us know where we should be getting our records from
     * @return the number of reviews written to disk
     */
    public int writeReviews(final VariantContextWriter writer, final SiteSelector selector) {
        return writeSelectedSites(writer, selector.onlyReviewed());
    }

    /**
     * Write all selected sites from the knowledge base to the writer
     *
     * @param writer
     * @param selector
     * @return
     */
    public int writeSelectedSites(final VariantContextWriter writer, final SiteSelector selector) {
        writer.writeHeader(makeStandardVCFHeader());

        int counter = 0;
        for ( final MongoVariantContext vc : getCalls(selector)) {
            writer.add(vc.getVariantContext());
            counter++;
        }

        return counter;
    }

    /**
     * @return a VCF header containing all of the standard NA12878 knowledge base metadata (INFO and FORMAT) fields
     */
    public VCFHeader makeStandardVCFHeader() {
        final Set<VCFHeaderLine> metadata = new HashSet<VCFHeaderLine>();

        for ( final VCFHeaderLine line : MongoVariantContext.reviewHeaderLines() )
            metadata.add(line);

        VCFStandardHeaderLines.addStandardFormatLines(metadata, true,
                VCFConstants.GENOTYPE_KEY, VCFConstants.DEPTH_KEY, VCFConstants.GENOTYPE_QUALITY_KEY);

        return new VCFHeader(metadata, Collections.singleton("NA12878"));
    }

    @Override
    public String toString() {
        String msg = String.format("NA12878KnowledgeBase{%s}", dblocator);
        return msg;
    }
}
