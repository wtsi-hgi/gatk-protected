package org.broadinstitute.sting.gatk.walkers.na12878kb;

import com.mongodb.*;
import org.apache.log4j.Logger;
import org.apache.log4j.Priority;
import org.broadinstitute.sting.gatk.walkers.na12878kb.errors.InvalidRecordsRemove;
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

    private MongoVariantContext makeConsensus(final Collection<MongoVariantContext> individualCalls) {
        final List<String> supportingCallSets = new LinkedList<String>();
        final VariantContext first = individualCalls.iterator().next().getVariantContext();

        final VariantContextBuilder builder = new VariantContextBuilder();
        builder.chr(first.getChr()).start(first.getStart()).stop(first.getEnd());

        final boolean isReviewed = isReviewed(individualCalls);
        final Set<Allele> alleles = new LinkedHashSet<Allele>();
        for ( final MongoVariantContext vc : individualCalls ) {
            alleles.addAll(vc.getVariantContext().getAlleles());
            supportingCallSets.addAll(vc.getSupportingCallSets());
        }
        builder.alleles(alleles);

        final TruthStatus type = determineTruth(individualCalls);
        final PolymorphicStatus status = determinePolymorphicStatus(individualCalls);
        final Genotype gt = consensusGT(type, status, new LinkedList<Allele>(alleles), individualCalls);

        return MongoVariantContext.create(supportingCallSets, builder.make(), type, new Date(), gt, isReviewed);
    }

    private boolean isReviewed(final Collection<MongoVariantContext> individualCalls) {
        boolean hasReview = false;
        for ( final MongoVariantContext vc : individualCalls )
            hasReview = hasReview || vc.isReviewed();
        return hasReview;
    }

    private Genotype consensusGT(final TruthStatus truthStatus,
                                 final PolymorphicStatus polyStatus,
                                 final List<Allele> alleles,
                                 final Collection<MongoVariantContext> individualCalls) {
        if ( ! truthStatus.isTruePositive() ) {
            return MongoGenotype.NO_CALL;
        } else if ( polyStatus.isDiscordant() || polyStatus.isUnknown() ) {
            return MongoGenotype.NO_CALL;
        } else {
            // we are a TP, we need to compute the consensus genotype
            for ( final MongoVariantContext mvc : individualCalls ) {
                return mvc.getGt().toGenotype(alleles); // TODO -- fix me
            }

            return MongoGenotype.NO_CALL;
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
        String msg = String.format("NA12878KnowledgeBase{db=%s, locator=%s}", dblocator.name, dblocator);
        return msg;
    }
}
