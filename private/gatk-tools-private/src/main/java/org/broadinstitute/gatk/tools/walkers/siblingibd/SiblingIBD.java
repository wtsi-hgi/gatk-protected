package org.broadinstitute.gatk.tools.walkers.siblingibd;

import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.engine.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.utils.SampleUtils;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.variant.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;

import java.io.PrintStream;
import java.util.*;

/**
 * Determines Identical-By-Descent regions in two siblings
 *
 * IBD0 regions - Siblings share neither parental haplotype
 * IBD1 regions - Siblings share one parental haplotype
 * IBD2 regions - Siblings share two parental haplotypes
 *
 * This tool implements two different methods for IBD region detection, which
 * should be used depending on whether the parental genotypes are available
 * (ie in a quartet sequencing project). Users can specify which type of analysis
 * to perform using the useParentalGenotypes parameter. If set to true, the tool
 * uses an HMM to label IBD regions based on the informativeness of each site (see
 * IDBObservation for details). If useParentalGenotypes is set to false, the tool
 * looks for sites at which minor alleles are shared between siblings and,
 * based on the level of sharing in sliding windows along the genome executes
 * a k-means model to determine the most likely IBD state combined
 * with median filters for smoothing. The latter method (based only on rare allele sharing)
 * is experimental and not recommended for general use.
 *
 * The tool scans the input pedigree file to detect sibling pairs, and computes
 * IBD regions for each pair of siblings detected in the pedigree file. We recommend
 * processing only one or two families at a time.
 *
 * IBD states are written to the output text file with the following columns:
 *
 * Sibling Pair ID: [Sibling Sample ID #1]-[Sibling Sample ID #2]
 * Chromosome
 * Region Start
 * Region End
 * IBD State: One of ZERO, ONE, or TWO
 */
public class SiblingIBD extends RodWalker<List<VariantIBDState>, Map<SiblingPair,IBDRegionAccumulator>> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Input(fullName="mask", shortName = "mask", doc="BED file of regions to mask variants in", required=true)
    public RodBinding<BEDFeature> mask = null;

    @Argument(shortName = "doNotUseParentalGenotypes",required = true,fullName = "doNotUseParentalGenotypes", doc="Don't use parental genotypes, just count minor allele sharing between two siblings")
    private boolean doNotUseParentalGenotypes = false;

    @Argument(shortName = "gqThreshold",required = false,fullName = "genotypeQualityThreshold", doc="Threshold on minimum GQ to include a site in the IBD calculation")
    private Integer gqThreshold = 50;

    @Argument(shortName = "slidingWindowSize",required = false,fullName = "slidingWindowSize", doc="Size of the sliding window to examine for shared alleles (Rare Allele Sharing Method Only)")
    private Integer slidingWindowSize = 1000;

    @Argument(shortName = "medianFilterSize",required = false,fullName = "medianFilterSize", doc="Size of median filter to use to clean up IBD class predictions (Rare Allele Sharing Method Only)")
    private Integer medianFilterSize = 500;

    @Argument(shortName = "ibdRegionsFile", required = false, doc = "File to which IBD regions should be written")
    protected PrintStream ibdRegionsFile = null;

    @Output(doc = "Output VCF file")
    protected VariantContextWriter output;

    @Argument(shortName = "useOriginalAF",required = false,fullName = "useOriginalAF", doc="Use the AF_Orig field rather than the AF field for allele frequency (Rare Allele Sharing Method Only)")
    private boolean useOriginalAF = false;

    @Argument(shortName = "spFile",required = false,fullName = "spFile", doc="File to output IBD observations at every unfiltered site for debugging purposes (Quartet Method Only)")
    private PrintStream spFile = null;

    @Argument(shortName = "countsFile",required = false,fullName = "countsFile", doc="File to output sliding window counts for debugging purposes (Rare Allele Sharing Method Only)")
    private PrintStream countsFile = null;

    @Argument(shortName = "unfilteredIbdClassFile",required = false,fullName = "unfilteredIbdClassFile", doc="File to output unfiltered IBD classes for debugging purposes (Rare Allele Sharing Method Only)")
    private PrintStream unfilteredIbdClassFile = null;

    @Argument(shortName = "filteredIbdClassFile",required = false,fullName = "filteredIbdClassFile", doc="File to output filtered IBD classes for debugging purposes (Rare Allele Sharing Method Only)")
    private PrintStream filteredIbdClassFile = null;

    private List<SiblingPair> siblingPairs;

    private IBDStateModel model = null;

    @Override
    public void initialize() {
        super.initialize();
        final ArrayList<String> rodNames = new ArrayList<>();
        rodNames.add(variantCollection.variants.getName());
        final Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        final Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        siblingPairs = initSiblingPairs(vcfSamples);

        if (doNotUseParentalGenotypes) {
            model = new SiblingRareAlleleSharingModel(siblingPairs, gqThreshold, countsFile, unfilteredIbdClassFile, filteredIbdClassFile, useOriginalAF, slidingWindowSize, medianFilterSize);
        } else {
            model = new QuartetIBDStateModel(siblingPairs, gqThreshold, spFile);
        }

        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        headerLines.add(new VCFFormatHeaderLine("IBDS", VCFHeaderLineCount.G, VCFHeaderLineType.String, "IBD Max Likelihood State"));
        headerLines.add(new VCFFormatHeaderLine("IBDQ", VCFHeaderLineCount.G, VCFHeaderLineType.Float, "IBD State Posteriors"));
        headerLines.add(new VCFHeaderLine("source", "SiblingIBD"));

        output.writeHeader(new VCFHeader(headerLines, vcfSamples));

    }

    private List<SiblingPair> initSiblingPairs(final Set<String> vcfSamples) {
        final List<SiblingPair> siblingPairs = new ArrayList<>();
        final Map<String,Set<Sample>> families = this.getSampleDB().getFamilies(vcfSamples);
        for(final String familyName : families.keySet()){
            final Set<Sample> family = families.get(familyName);
            final Map<String, List<Sample>> familyMembersByParents = new HashMap<>();

            // first build up the list of family members by parents
            for(final Sample familyMember : family) {
                if (vcfSamples.contains(familyMember.getID())) {
                    final List<Sample> parents = familyMember.getParents();
                    if (parents != null && parents.size() == 2) {
                        final String parentPairId = parents.get(0).getID() + "-" + parents.get(1).getID();
                        if (!familyMembersByParents.containsKey(parentPairId)) {
                            familyMembersByParents.put(parentPairId, new ArrayList<Sample>());
                        }
                        familyMembersByParents.get(parentPairId).add(familyMember);
                    }
                }
            }

            // now find all sibling pairs in the family
            for (final String parentPairId : familyMembersByParents.keySet()) {
                final List<Sample> siblings = familyMembersByParents.get(parentPairId);
                for (int i  = 0; i < siblings.size(); i++) {
                    for (int j = i + 1; j < siblings.size(); j++) {
                        final SiblingPair siblingPair = new SiblingPair(siblings.get(i), siblings.get(j));
                        siblingPairs.add(siblingPair);
                    }
                }
            }
        }
        return siblingPairs;
    }


    @Override
    public List<VariantIBDState> map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if (tracker == null)
            return null;

        final VariantContext vc = tracker.getFirstValue(variantCollection.variants, context.getLocation());
        if ( vc == null )
            return null;

        if ( !vc.isBiallelic()) {
            //vcsToWrite.add(vc);
            return null;
        }

        if (! vc.isSNP()) {
            //vcsToWrite.add(vc);
            return null;
        }

        if (vc.isFiltered()) {
            //vcsToWrite.add(vc);
            return null;
        }

        if (tracker.getFirstValue(mask) != null) {
            //vcsToWrite.add(vc);
            return null;
        }

        final double af = SiblingRareAlleleSharingModel.getAF(useOriginalAF, vc);
        if (doNotUseParentalGenotypes && af > .4 && af < .6) {
            //vcsToWrite.add(vc);
            return null;
        }

        final List<VariantIBDState> result = new ArrayList<>();
        final Map<Integer, List<VariantIBDState>> ibdStateModifications = new HashMap<>();
        for (final SiblingPair siblingPair : siblingPairs) {
            final Genotype sib1Gt = vc.getGenotype(siblingPair.sib1.getID());
            final Genotype sib2Gt = vc.getGenotype(siblingPair.sib2.getID());

            if (sib1Gt.isCalled() && sib2Gt.isCalled()) {
                final List<VariantIBDState> modelResults = model.addSite(vc, siblingPair, sib1Gt, sib2Gt);
                if (modelResults.size() == 0) {
                    continue;
                }
                for (final VariantIBDState variantIBDState : modelResults) {
                    if (! ibdStateModifications.containsKey(variantIBDState.vc.getStart())) {
                        ibdStateModifications.put(variantIBDState.vc.getStart(), new ArrayList<VariantIBDState>());
                    }
                    ibdStateModifications.get(variantIBDState.vc.getStart()).add(variantIBDState);
                }
                result.addAll(modelResults);
            }
        }

        if (ibdStateModifications.keySet().size() == 0) {
            return null;
        }

        for (final List<VariantIBDState> variantIDBStates : ibdStateModifications.values()) {
            if (variantIDBStates.size() == 0) {
                continue;
            }
            final VariantContext modifiedVC = variantIDBStates.get(0).vc;
            final VariantContextBuilder builder = new VariantContextBuilder(modifiedVC);
            final VariantContext newVC = addIBDStatesToVC(builder, modifiedVC, variantIDBStates);
            output.add(newVC);
        }

        if (logger.isDebugEnabled() && !result.isEmpty()) {
            logger.debug("emitting " + result);
        }
        return result;
    }

    private VariantContext addIBDStatesToVC(final VariantContextBuilder builder, final VariantContext vc, final List<VariantIBDState> variantIBDStates) {
        final Map<String, String> ibdStateAttributes = new HashMap<>();
        final Map<String, String> ibdPosteriorAttributes = new HashMap<>();

        for (final VariantIBDState variantIBDState : variantIBDStates) {
            final SiblingPair siblingPair = variantIBDState.siblingPair;
            final Sample s1 = siblingPair.sib1;
            final Sample s2 = siblingPair.sib2;

            addIBDAttributesToMaps(s1, ibdStateAttributes, ibdPosteriorAttributes, variantIBDState);
            addIBDAttributesToMaps(s2, ibdStateAttributes, ibdPosteriorAttributes, variantIBDState);
        }

        final GenotypesContext newContext = GenotypesContext.create();
        for (final Genotype gt : vc.getGenotypes()) {
            GenotypeBuilder gtBuilder = new GenotypeBuilder(gt);
            if (ibdStateAttributes.containsKey(gt.getSampleName())) {
                gtBuilder = gtBuilder.attribute("IBDS", ibdStateAttributes.get(gt.getSampleName()));
            }
            if (ibdStateAttributes.containsKey(gt.getSampleName())) {
                gtBuilder = gtBuilder.attribute("IBDQ", ibdPosteriorAttributes.get(gt.getSampleName()));
            }
            newContext.add(gtBuilder.make());
        }
        return builder.genotypes(newContext).make();
    }

    private void addIBDAttributesToMaps(final Sample s1, final Map<String, String> ibdStateAttributes, final Map<String, String> ibdPosteriorAttributes, final VariantIBDState variantIBDState) {
        final String ibdPairStateString = variantIBDState.siblingPair.sib1.getID() + "-" + variantIBDState.siblingPair.sib2.getID() + "=" + variantIBDState.ibdState.ordinal();
        if (ibdStateAttributes.containsKey(s1.getID())) {
            ibdStateAttributes.put(s1.getID(), ibdStateAttributes.get(s1.getID()) + ";" + ibdPairStateString);
        } else {
            ibdStateAttributes.put(s1.getID(), ibdPairStateString);
        }
        if (variantIBDState.statePosteriors != null) {
            final String ibdPosteriorString = variantIBDState.siblingPair.sib1.getID() + "-" + variantIBDState.siblingPair.sib2.getID() + "=" + formatPosteriors(variantIBDState);
            if (ibdPosteriorAttributes.containsKey(s1.getID())) {
                ibdPosteriorAttributes.put(s1.getID(), ibdPosteriorAttributes.get(s1.getID()) + ";" + ibdPosteriorString  );
            } else {
                ibdPosteriorAttributes.put(s1.getID(), ibdPosteriorString);
            }
        }
    }

    private String formatPosteriors(final VariantIBDState variantIBDState) {
        final StringBuilder builder = new StringBuilder();
        for (int i = 0; i < variantIBDState.statePosteriors.length; i++) {
            final double p = variantIBDState.statePosteriors[i];
            builder.append(new Double(-1 * p).intValue());
            if (i < variantIBDState.statePosteriors.length - 1) {
                builder.append(",");
            }
        }
        return builder.toString();
    }

    @Override
    public Map<SiblingPair, IBDRegionAccumulator> reduceInit() {
        final Map<SiblingPair, IBDRegionAccumulator> ibdRegionAccumulators = new HashMap<>();
        for (final SiblingPair siblingPair : siblingPairs) {
            ibdRegionAccumulators.put(siblingPair, new IBDRegionAccumulator());
        }
        return ibdRegionAccumulators;
    }

    @Override
    public Map<SiblingPair, IBDRegionAccumulator> reduce(final List<VariantIBDState> valueList, final Map<SiblingPair, IBDRegionAccumulator> accumulators) {
        if (valueList != null) {
            accumulateValues(valueList, accumulators);
        }
        return accumulators;
    }

    private void accumulateValues(final List<VariantIBDState> valueList, final Map<SiblingPair, IBDRegionAccumulator> accumulators) {
        if (ibdRegionsFile != null) {
            for (final VariantIBDState value : valueList) {
                if (value != null) {
                    final SiblingPair siblingPair = value.siblingPair;
                    final IBDRegionAccumulator.IBDRegion region = accumulators.get(siblingPair).regionChange(value.vc.getChr(), value.vc.getStart(), value.ibdState);
                    if (region != null) {
                        ibdRegionsFile.print(siblingPair.getName() + "\t" + region.chr + "\t" + region.start + "\t" + region.end + "\t" + region.state + "\n");
                    }
                }
            }
        }
    }

    @Override
    public void onTraversalDone(final Map<SiblingPair, IBDRegionAccumulator> accumulators) {
        final Map<Integer, List<VariantIBDState>> ibdStateModifications = new HashMap<>();
        final List<VariantIBDState> finalValues = model.finalizeModel();
        for (final VariantIBDState variantIBDState : finalValues) {
            if (! ibdStateModifications.containsKey(variantIBDState.vc.getStart())) {
                ibdStateModifications.put(variantIBDState.vc.getStart(), new ArrayList<VariantIBDState>());
            }
            ibdStateModifications.get(variantIBDState.vc.getStart()).add(variantIBDState);
        }

        for (final List<VariantIBDState> variantIDBStates : ibdStateModifications.values()) {
            if (variantIDBStates.size() == 0) {
                continue;
            }
            final VariantContext modifiedVC = variantIDBStates.get(0).vc;
            final VariantContextBuilder builder = new VariantContextBuilder(modifiedVC);
            final VariantContext newVC = addIBDStatesToVC(builder, modifiedVC, variantIDBStates);
            output.add(newVC);
        }

        accumulateValues(finalValues, accumulators);
        if (ibdRegionsFile != null) {
            for (final SiblingPair siblingPair : accumulators.keySet()) {
                final IBDRegionAccumulator accumulator = accumulators.get(siblingPair);
                final IBDRegionAccumulator.IBDRegion region = accumulator.getFinalRegion();
                ibdRegionsFile.print(siblingPair.getName() + "\t" + region.chr + "\t" + region.start + "\t" + region.end + "\t" + region.state + "\n");
            }
        }
        output.close();
    }


}
