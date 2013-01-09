package org.broadinstitute.sting.gatk.walkers.na12878kb;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.NA12878DBArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.SiteIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
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
     * Variants from these VCF files are used by this tool as input.
     * The files must at least contain the standard VCF header lines, but
     * can be empty (i.e., no variants are contained in the file).
     */
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public List<RodBinding<VariantContext>> variants;

    @Input(fullName="BAM", shortName = "BAM", doc="Input BAM file.  If provided, we will differentiate false negative sites into those truly missed and those without coverage", required=false)
    public File BAM = null;

    @Output(doc="Summary GATKReport will be written here")
    public PrintStream out;

    @Argument(fullName="excludeCallset", shortName = "excludeCallset", doc="Don't count calls that come from only these excluded callsets", required=false)
    public Set<String> excludeCallset = null;

    /**
     * An output VCF file containing the bad sites (FN/FP) that were found in the input callset w.r.t. the current NA12878 knowledge base
     */
    @Output(fullName = "badSites", shortName = "badSites", doc="VCF file containing information on FP/FNs in the input callset", required=false)
    public VariantContextWriter badSites = null;

    @Argument(fullName="maxToWrite", shortName = "maxToWrite", doc="Max. number of bad sites to write out", required=false)
    public int maxToWrite = 10000;

    @Argument(fullName="minDepthForLowCoverage", shortName = "minDepthForLowCoverage", doc="A false negative will be flagged as due to low coverage if the (optional) BAM is provided and the coverage overlapping the site is less than this value", required=false)
    public int minDepthForLowCoverage = 5;

    @Argument(fullName="typesToInclude", shortName = "typesToInclude", doc="Should we analyze SNPs, INDELs, or both?", required=false)
    public TypesToInclude typesToInclude = TypesToInclude.BOTH;

    public enum TypesToInclude {
        SNPS,
        INDELS,
        BOTH
    }

    /**
     * Useful when some state isn't interesting as a bad site but is extremely prevalent in the input callset
     */
    @Argument(fullName="AssessmentsToExclude", shortName = "AssessmentsToExclude", doc="If provided, we will prevent any of these states from being written out to the badSites VCF.", required=false)
    public Set<AssessmentType> AssessmentsToExclude = EnumSet.noneOf(AssessmentType.class);

    @Hidden
    @Argument(shortName = "debug", required=false)
    protected boolean debug = false;

    SiteIterator<MongoVariantContext> consensusSiteIterator;
    boolean captureBadSites = true;
    int nWritten = 0;

    public enum AssessmentType {
        TRUE_POSITIVE(false),
        CORRECTLY_FILTERED(false),
        CORRECTLY_UNCALLED(false),
        REASONABLE_FILTERS_WOULD_FILTER_FP_SITE(false),
        FALSE_POSITIVE_SITE_IS_FP(true),
        FALSE_POSITIVE_MONO_IN_NA12878(true),
        FALSE_NEGATIVE_CALLED_BUT_FILTERED(true),
        FALSE_NEGATIVE_NOT_CALLED_BUT_LOW_COVERAGE(false),
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

        public List<Integer> getCounts() {
            final List<Integer> returnList = new ArrayList<Integer>();
            for( final AssessmentType type : AssessmentType.values() ) {
                returnList.add(get(type));
            }
            return returnList;
        }
    }

    final Map<String,Assessment> SNPAssessments = new HashMap<String,Assessment>();
    final Map<String,Assessment> IndelAssessments = new HashMap<String,Assessment>();
    SAMFileReader bamReader = null;

    @Override
    public NA12878DBArgumentCollection.DBType getDefaultDB() {
        return NA12878DBArgumentCollection.DBType.PRODUCTION;
    }

    public void initialize() {
        super.initialize();
        consensusSiteIterator = db.getConsensusSites(makeSiteSelector());
        captureBadSites = badSites != null;

        if( captureBadSites ) {
            final Set<VCFHeaderLine> lines = GATKVCFUtils.getHeaderFields(getToolkit());
            lines.add(new VCFInfoHeaderLine("WHY", 1, VCFHeaderLineType.String, "Why was the site considered bad"));
            lines.add(new VCFInfoHeaderLine("SupportingCallsets", 1, VCFHeaderLineType.String, "Callsets supporting the consensus, where available"));
            lines.addAll(MongoVariantContext.reviewHeaderLines());
            badSites.writeHeader(new VCFHeader(lines, Collections.singleton("NA12878")));
        }

        if ( BAM != null ) {
            bamReader = new SAMFileReader(BAM);
        }
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) return 0;

        if ( debug ) logger.info("Processing " + context.getLocation());
        includeMissingCalls(consensusSiteIterator, context.getLocation());
        final List<MongoVariantContext> consensusSites = consensusSiteIterator.getSitesAtLocation(context.getLocation());

        for ( final RodBinding<VariantContext> rod : variants ) {
            if( tracker.getValues(rod, ref.getLocus()).size() == 0 ) {
                // missed consensus site(s)
                for ( final MongoVariantContext site : consensusSites ) {
                    accessSite(rod.getName(), null, site);
                }
            } else {
                for( final VariantContext vcRaw : tracker.getValues(rod, ref.getLocus()) ) {
                    final VariantContext vc = vcRaw.subContextFromSample("NA12878");

                    if ( ! vc.isBiallelic() ) {
                        logger.info("Skipping unsupported multi-allelic variant " + vc);
                        continue;
                    }

                    accessSite(rod.getName(), vc, findMatching(vc, consensusSites)); // BUGBUG: Should this be called for not only the matching consensusSite but for all consensusSites?
                }
            }
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
            for ( final RodBinding<VariantContext> rod : variants ) {
                accessSite(rod.getName(), null, missedSite);
            }
        }
    }

    /**
     * Should we include the variant vc in our assessment?
     *
     * @param vc a VariantContext to potentially include
     * @return true if VC should be included, false otherwise
     */
    private boolean includeVariant(final VariantContext vc) {
        switch ( typesToInclude ) {
            case BOTH: return true;
            case SNPS: return vc.isSNP();
            case INDELS: return ! vc.isSNP();
            default:
                throw new IllegalStateException("Unexpected enum " + typesToInclude);
        }
    }

    private void accessSite(final String rodName, final VariantContext call, final MongoVariantContext consensusSite) {

        final VariantContext vc = call != null ? call : consensusSite.getVariantContext();

        if ( ! includeVariant(vc) )
            return;

        final AssessmentType type = figureOutAssessmentType(call, consensusSite);
        final Map<String,Assessment> assessmentMap = vc.isSNP() ? SNPAssessments : IndelAssessments;
        Assessment assessment = assessmentMap.get(rodName);
        if(assessment == null) {
            assessment = new Assessment();
            assessmentMap.put(rodName, assessment);
        }
        assessment.inc(type);

        if ( captureBadSites && type.interesting && ! AssessmentsToExclude.contains(type) && nWritten++ < maxToWrite) {
            final VariantContextBuilder builder = new VariantContextBuilder(vc);
            if ( consensusSite != null )
                builder.attribute("SupportingCallsets", Utils.join(",", consensusSite.getSupportingCallSets()));
            builder.attribute("WHY", type.toString());
            badSites.add(builder.make());
            logger.info("Accessed site " + call + " consensus " + consensusSite);
        }
    }

    private AssessmentType figureOutAssessmentType(final VariantContext call, final MongoVariantContext consensusSite) {
        final boolean consensusTP = consensusSite != null && !isExcluded(consensusSite) && consensusSite.getType().isTruePositive() && consensusSite.getPolymorphicStatus().isPolymorphic();
        final boolean consensusFP = consensusSite != null && !isExcluded(consensusSite) && consensusSite.getType().isFalsePositive();

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
        } else if ( consensusTP ) { // call == null
            if ( BAM != null ) {
                return sufficientDepthToCall(consensusSite)
                        ? AssessmentType.FALSE_NEGATIVE_NOT_CALLED_AT_ALL
                        : AssessmentType.FALSE_NEGATIVE_NOT_CALLED_BUT_LOW_COVERAGE;
            } else {
                return AssessmentType.FALSE_NEGATIVE_NOT_CALLED_AT_ALL;
            }
        } else if (consensusFP ) { // call == null
            return AssessmentType.CORRECTLY_UNCALLED;
        } else {
            return AssessmentType.NOT_RELEVANT;
        }
    }

    private boolean isExcluded( final MongoVariantContext consensusSite ) {
        return excludeCallset != null && excludeCallset.containsAll(consensusSite.getSupportingCallSets());
    }

    /**
     * If a BAM is provided, assess whether there's sufficient coverage to call the site
     *
     * @param falseNegative a site that was missed in the call set
     * @return true if there's enough coverage at the site in the BAM to likely make a call
     */
    private boolean sufficientDepthToCall(final MongoVariantContext falseNegative) {
        final SAMRecordIterator it = bamReader.queryOverlapping(falseNegative.getChr(), falseNegative.getStart(), falseNegative.getStart());

        int depth = 0;
        while ( it.hasNext() && depth < minDepthForLowCoverage ) {
            final SAMRecord read = it.next();
            // TODO -- filter for MAPQ?
            depth++;
        }
        it.close();

        return depth >= minDepthForLowCoverage;
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

        if( variants.size() == 1 ) {
            final GATKReport report = GATKReport.newSimpleReportWithDescription("NA12878Assessment", "Evaluation of input variant callsets",
                    "Name", "VariantType", "AssessmentType", "Count");
            for( final RodBinding rod : variants ) {
                // snps
                for ( final AssessmentType type : AssessmentType.values() ) {
                    report.addRow(rod.getName(), "SNP", type, SNPAssessments.get(rod.getName()).get(type));
                }

                // indels
                for ( final AssessmentType type : AssessmentType.values() ) {
                    report.addRow(rod.getName(), "Indel", type, IndelAssessments.get(rod.getName()).get(type));
                }
            }
            report.print(out);
        } else {
            final List<String> columns = new ArrayList<String>();
            columns.add("Name");
            columns.add("VariantType");
            for( final AssessmentType type : AssessmentType.values() ) {
                columns.add(type.toString());
            }
            final GATKReport report = GATKReport.newSimpleReport("NA12878Assessment", columns);
            for( final RodBinding rod : variants ) {
                // snps
                final List<Object> row = new ArrayList<Object>();
                row.add(rod.getName());
                row.add("SNP");
                row.addAll(SNPAssessments.get(rod.getName()).getCounts());
                report.addRowList(row);

                // indels
                row.clear();
                row.add(rod.getName());
                row.add("indel");
                row.addAll(IndelAssessments.get(rod.getName()).getCounts());
                report.addRowList(row);
            }
            report.print(out);
        }
    }
}