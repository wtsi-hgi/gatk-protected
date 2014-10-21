/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.siblingibd;

import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
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
 *
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
