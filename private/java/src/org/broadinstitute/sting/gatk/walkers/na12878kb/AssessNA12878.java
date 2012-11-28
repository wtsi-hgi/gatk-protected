package org.broadinstitute.sting.gatk.walkers.na12878kb;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;

import java.io.PrintStream;
import java.util.*;

/**
 * Assess the quality of an NA12878 callset against the NA12878 knowledge base
 *
 * <p>
 *     This walker takes a single VCF file contains calls of any type (SNPs, Indels)
 *     using at least the sample NA12878 (i.e., the VCF must contain the sample NA12878)
 *     and provides an itemized summary of the different types of true positives, false positives
 *     and false negatives in the input callset relative to the NA12878 knowledge base.
 * </p>Additionally, writes out a bad sites VCF that contains the following data by default: calls at known
 *     false positives, calls in NA12878 known to be monomorphic, sites that are TP but are
 *     either called but filtered out or not called at all, and finally calls in the input
 *     VCF not in the DB at all or in the DB with unknown status.  This VCF contains INFO field
 *     key/value pairs describing why the site was included (i.e., it is a false negative).
 * </p>
 *
 * See http://gatkforums.broadinstitute.org/discussion/1848/using-the-na12878-knowledge-base for more information.
 *
 * @author depristo
 * @since 11/2012
 * @version 0.1
 */
public class AssessNA12878 extends NA12878DBWalker {
    /**
     * Variants from this VCF file are used by this tool as input.
     * The file must at least contain the standard VCF header lines, but
     * can be empty (i.e., no variants are contained in the file).
     */
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=false)
    public RodBinding<VariantContext> variants;

    @Output(doc="Summary GATKReport will be written here")
    public PrintStream out;

    /**
     * A VCF file containing bad sites (FN/FP) in the input callset w.r.t. the current NA12878 knowledge base
     */
    @Output(fullName = "badSites", shortName = "badSites", doc="VCF file containing information on FP/FNs in the input callset")
    public VariantContextWriter badSites;

    @Input(fullName="maxToWrite", shortName = "maxToWrite", doc="Max. number of bad sites to write out", required=false)
    public int maxToWrite = 10000;

    /**
     * Useful when some state isn't interesting as a bad site but is extremely prevalent in the input callset
     */
    @Input(fullName="AssessmentsToExclude", shortName = "AssessmentsToExclude", doc="If provided, we will prevent any of these states from being written out to the badSites VCF.", required=false)
    public Set<AssessmentType> AssessmentsToExclude = EnumSet.noneOf(AssessmentType.class);

    @Hidden
    @Argument(shortName = "debug", required=false)
    protected boolean debug = false;

    SiteIterator<MongoVariantContext> consensusSiteIterator;
    int nWritten = 0;

    private enum AssessmentType {
        TRUE_POSITIVE(false),
        CORRECTLY_FILTERED(false),
        REASONABLE_FILTERS_WOULD_FILTER_FP_SITE(false),
        FALSE_POSITIVE_SITE_IS_FP(true),
        FALSE_POSITIVE_MONO_IN_NA12878(true),
        FALSE_NEGATIVE_CALLED_BUT_FILTERED(true),
        FALSE_NEGATIVE_NOT_CALLED_AT_ALL(true),
        CALLED_IN_DB_UNKNOWN_STATUS(true),
        CALLED_NOT_IN_DB_AT_ALL(true),
        NOT_RELEVANT(false);

        private final boolean interesting;

        private AssessmentType(boolean interesting) {
            this.interesting = interesting;
        }
    }

    private class Assessment {
        final EnumMap<AssessmentType, Integer> counts;

        public Assessment() {
            counts = new EnumMap<AssessmentType, Integer>(AssessmentType.class);
            for ( final AssessmentType type : AssessmentType.values() )
                counts.put(type, 0);
        }

        public final void inc(final AssessmentType type) {
            counts.put(type, counts.get(type) + 1);
        }

        public final int get(final AssessmentType type) {
            return counts.get(type);
        }
    }

    final Assessment SNPAssessments = new Assessment();
    final Assessment IndelAssessments = new Assessment();

    @Override
    public NA12878DBArgumentCollection.DBType getDefaultDB() {
        return NA12878DBArgumentCollection.DBType.PRODUCTION;
    }

    public void initialize() {
        super.initialize();
        consensusSiteIterator = db.getConsensusSites(makeSiteSelector());

        final Set<VCFHeaderLine> lines = VCFUtils.getHeaderFields(getToolkit());
        lines.add(new VCFInfoHeaderLine("WHY", 1, VCFHeaderLineType.String, "Why was the site considered bad"));
        lines.add(new VCFInfoHeaderLine("SupportingCallsets", 1, VCFHeaderLineType.String, "Callsets supporting the consensus, where available"));
        lines.add(new VCFHeaderLine("CallSetBeingEvaluated", variants.getSource()));
        lines.addAll(MongoVariantContext.reviewHeaderLines());
        badSites.writeHeader(new VCFHeader(lines, Collections.singleton("NA12878")));
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) return 0;

        if ( debug ) logger.info("Processing " + context.getLocation());
        includeMissingCalls(consensusSiteIterator, context.getLocation());
        final List<MongoVariantContext> consensusSites = consensusSiteIterator.getSitesAtLocation(context.getLocation());

        for ( VariantContext vcRaw : tracker.getValues(variants, ref.getLocus()) ) {
            final VariantContext vc = vcRaw.subContextFromSample("NA12878");

            if ( ! vc.isBiallelic() ) {
                logger.info("Skipping unsupported multi-allelic variant " + vc);
                continue;
            }

            final MongoVariantContext matchedConsensus = findMatching(vc, consensusSites);
            accessSite(vc, matchedConsensus);
        }

        return 1;
    }

    final MongoVariantContext findMatching(final VariantContext vc, final Collection<MongoVariantContext> consensusSites ) {
        for ( final MongoVariantContext site : consensusSites )
            if ( site.matches(vc) )
                return site;
        return null;
    }

    private void includeMissingCalls(final SiteIterator<MongoVariantContext> siteIterator, final GenomeLoc loc) {
        includeMissingCalls(siteIterator.getSitesBefore(loc));
    }

    private void includeMissingCalls(final List<MongoVariantContext> missedSites) {
        for ( final MongoVariantContext missedSite : missedSites ) {
            logger.info("Missed site " + missedSite);
            accessSite(null, missedSite);
        }
    }

    private void accessSite(final VariantContext call, final MongoVariantContext consensusSite) {
        final VariantContext vc = call != null ? call : consensusSite.getVariantContext();
        final AssessmentType type = figureOutAssessmentType(call, consensusSite);
        final Assessment assessment = vc.isSNP() ? SNPAssessments : IndelAssessments;
        assessment.inc(type);

        if ( type.interesting && ! AssessmentsToExclude.contains(type) && nWritten++ < maxToWrite) {
            final VariantContextBuilder builder = new VariantContextBuilder(vc);
            if ( consensusSite != null )
                builder.attribute("SupportingCallsets", Utils.join(",", consensusSite.getSupportingCallSets()));
            builder.attribute("WHY", type.toString());
            badSites.add(builder.make());
        }

        logger.info("Accessed site " + call + " consensus " + consensusSite);
    }

    private AssessmentType figureOutAssessmentType(final VariantContext call, final MongoVariantContext consensusSite) {
        final boolean consensusTP = consensusSite != null && consensusSite.getType().isTruePositive() && consensusSite.getPolymorphicStatus().isPolymorphic();
        final boolean consensusFP = consensusSite != null && consensusSite.getType().isFalsePositive();

        if ( call != null ) {
            if ( consensusTP ) {
                return call.isFiltered()
                        ? AssessmentType.FALSE_NEGATIVE_CALLED_BUT_FILTERED
                        : AssessmentType.TRUE_POSITIVE;
            } else if (consensusSite != null && consensusSite.getType().isTruePositive() && consensusSite.getPolymorphicStatus().isMonomorphic() ) {
                return AssessmentType.FALSE_POSITIVE_MONO_IN_NA12878;
            } else if ( consensusFP ) {
                if ( call.isFiltered() )
                    return AssessmentType.CORRECTLY_FILTERED;
                else if ( likelyWouldBeFiltered(call) )
                    return AssessmentType.REASONABLE_FILTERS_WOULD_FILTER_FP_SITE;
                else
                    return AssessmentType.FALSE_POSITIVE_SITE_IS_FP;
            } else if ( consensusSite != null && consensusSite.getType().isUnknown() ) {
                return AssessmentType.CALLED_IN_DB_UNKNOWN_STATUS;
            } else if ( consensusSite == null && call.isNotFiltered() ) {
                return AssessmentType.CALLED_NOT_IN_DB_AT_ALL;
            } else {
                return AssessmentType.NOT_RELEVANT;
            }
        } else {
            // maybe false negative
            return consensusTP ? AssessmentType.FALSE_NEGATIVE_NOT_CALLED_AT_ALL : AssessmentType.NOT_RELEVANT;
        }
    }

    /**
     * Returns true is a simple set of reasonable filters would likely remove it
     *
     * For SNPs:
     * QD < 2.0
     * MQ < 40.0
     * FS > 60.0
     * HaplotypeScore > 13.0
     * MQRankSum < -12.5
     * ReadPosRankSum < -8.0
     *
     * For indels:
     *
     * QD < 2.0
     * ReadPosRankSum < -20.0
     * InbreedingCoeff < -0.8
     * FS > 200.0
     *
     * @param vc
     * @return
     */
    public boolean likelyWouldBeFiltered(final VariantContext vc) {
        final double FS = vc.getAttributeAsDouble("FS", 0.0);
        final double QD = vc.getAttributeAsDouble("QD", 20.0);
        final double MQ = vc.getAttributeAsDouble("MQ", 50.0);

        if ( vc.isSNP() ) {
            return FS > 60 || QD < 2 || MQ < 40;
        } else {
            return FS > 200 || QD < 2;
        }
    }

    public void onTraversalDone(Integer result) {
        super.onTraversalDone(result);
        includeMissingCalls(consensusSiteIterator.toList());

        final GATKReport report = GATKReport.newSimpleReportWithDescription("NA12878Assessment", "Evaluation of " + variants.getSource(),
                "VariantType", "AssessmentType", "Count");

        // snps
        for ( final AssessmentType type : AssessmentType.values() ) {
            report.addRow("SNP", type, SNPAssessments.get(type));
        }

        // indels
        for ( final AssessmentType type : AssessmentType.values() ) {
            report.addRow("Indel", type, IndelAssessments.get(type));
        }

        report.print(out);
    }
}