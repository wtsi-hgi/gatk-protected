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
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.na12878kb.assess;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.TruthStatus;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.io.Resource;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.tools.walkers.na12878kb.NA12878DBWalker;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.NA12878DBArgumentCollection;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.SiteIterator;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.R.RScriptExecutor;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.io.*;
import java.util.*;

/**
 * Create a ROC curve by walking over the NA12878 KB and ranking the input variants by their VQSLOD score
 * User: rpoplin
 * Date: 12/10/13
 */

public class ROCCurveNA12878 extends NA12878DBWalker {

    final static int SNP_INDEX = 0, INDEL_INDEX = 1;
    final static int TP_INDEX = 0, FP_INDEX = 1;
    final static int CALLED_INDEX = 0, TOTAL_INDEX = 1;
    final static AlignmentContext lastLoc = new AlignmentContext(GenomeLoc.END_OF_GENOME, new ReadBackedPileupImpl(GenomeLoc.END_OF_GENOME));

    /**
     * Variants from these VCF files are used by this tool as input.
     * The files must at least contain the standard VCF header lines, but
     * can be empty (i.e., no variants are contained in the file).
     */
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public RodBinding<VariantContext> variants;

    @Argument(fullName="N_SNP_Bins", shortName = "n_snp_bins", doc="number of bins to use for SNPs in making the ROC curve", required=false)
    private int n_bin_snp = 100;

    @Argument(fullName="N_Indel_Bins", shortName = "n_indel_bins", doc="number of bins to use for Indels in making the ROC curve", required=false)
    private int n_bin_indel = 50;

    @Output(doc="Summary GATKReport will be written here", required=false)
    public static File report_file;

    @Argument(fullName="project", shortName = "project", doc="String project tag", required=true)
    public String project = null;

    @Argument(fullName="requireReviewed", shortName = "requireReviewed", doc="If specified will only use reviewed sites in the knowledgebase for the assessment", required=false)
    public boolean REQUIRE_REVIEWED = false;

    /**
     * Unknown sites in the knowledgebase (like unreviewed Mills sites) are likely to be true positives and will be counted as such if this argument is used.
     * The UNKNOWN status indicates that the variant has been seen in someone, but it is unknown if this variant is in NA12878.
     * This argument must be specified in order to achieve concordance of TPs and FPs with AssessNA12878
     */
    @Argument(fullName="includeUnknownTruthSites", shortName = "includeUnknowns", doc="If specified, sites that are listed as UNKNOWN in the knowledgebase will be counted as TPs")
    public boolean includeUnknowns = false;

    @Argument(fullName="highConfidenceMode", shortName = "hiConf", doc="The name of the database that should be used if running within a high confidence region given by -L", required=false)
    public String hiConf = null;

    @Argument(fullName="sampleNameToCompare", shortName = "sample", doc="Specifies the sample name in the VCF that will be compared to the knowledgebase.", required=false)
    public String sampleNameToCompare = "NA12878";

    @Output(fullName="plotFiles", shortName="plotFiles", doc="The base name of two plots files (TP vs FP counts and ratio) output by the R-script", required=false, defaultToStdout=false)
    private static String PLOTS_BASENAME = null;

    @Argument(fullName="maxToWrite", shortName = "maxToWrite", doc="Max. number of bad sites to write out", required=false)
    public int maxToWrite = 100_000_000;

    /**
     * An output VCF file containing the used sites (TP/FP) that were found in the input callset w.r.t. the current NA12878 knowledge base
     */
    @Output(fullName = "rocSites", shortName = "rocSites", doc="VCF file containing information on TP/FPss in the input callset", required=false, defaultToStdout=false)
    public VariantContextWriter rocSites = null;

    private SiteIterator<MongoVariantContext> siteIterator;
    private List<ROCDatum> data = new ArrayList<>();
    private int[][] uncalledSites = {{0,0},{0,0}};
    private static final String PLOT_RSCRIPT = "roc_Plot.R";
    public NA12878DBArgumentCollection.DBType getDefaultDB() {
        return NA12878DBArgumentCollection.DBType.PRODUCTION;
    }

    private SitesWriter sitesWriter;
    private Set<AssessmentType> AssessmentsToExclude;


    @Override
    public void initialize() {
        super.initialize();

        AssessmentsToExclude = new HashSet<>();
        AssessmentsToExclude.add(AssessmentType.FALSE_NEGATIVE);
        AssessmentsToExclude.add(AssessmentType.TRUE_NEGATIVE);

        if(hiConf!=null){
            if(db.getCallSet(hiConf)==null){
                throw new UserException(hiConf + "callset name is not included in the database.");
            }

            if(this.getToolkit().getIntervals() == null){
                throw new UserException("hiConfMode set but no interval list given. An interval list must be provided with -L of the high confidence region.");
            }

            // We need to get calls here instead of consensus sites because even if a site has a supporting callset of
            // hiConf, its result might disagree with the consensus.
            siteIterator = db.getCalls(makeSiteManager(false));
        }else{
            siteIterator = db.getConsensusSites(makeSiteManager(false));
        }

        if (PLOTS_BASENAME != null && !RScriptExecutor.RSCRIPT_EXISTS)
            Utils.warnUser(logger, String.format(
                    "Rscript %s not found in environment path. Plots %s will not be generated.",
                    PLOT_RSCRIPT, PLOTS_BASENAME));

        if (rocSites == null)
            sitesWriter = SitesWriter.NOOP_WRITER;
        else
            sitesWriter = new AllSitesWriter(maxToWrite, AssessmentsToExclude, rocSites);
        sitesWriter.initialize(GATKVCFUtils.getHeaderFields(getToolkit()));
    }

    @Override
    public Integer map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if ( tracker == null ) return 0;

        countUncalledBefore(context);

        List<MongoVariantContext> kbSitesAtThisLoc = siteIterator.getSitesAtLocation(context.getLocation());
        // Does this input site overlap a KB site
        for( final MongoVariantContext cs : kbSitesAtThisLoc ) {
            TruthStatus csType = cs.getType();
            //Check if we are in hiConf mode that this MongoVariantContext is from the correct callSet
            if(hiConf!=null && !cs.getSupportingCallSets().contains(hiConf)){
                continue;
            }

            // Is it either a SNP or indel and is it marked as either a true positive or false positive
            if( (cs.getVariantContext().isSNP() || cs.getVariantContext().isIndel()) && (cs.getType().isTruePositive() || cs.getType().isFalsePositive()) || (includeUnknowns && cs.getType().isUnknown()) ) {
                if( !REQUIRE_REVIEWED || cs.isReviewed() ) {
                    boolean foundMatchingAllele = false;
                    for ( final VariantContext vcRaw : tracker.getValues(variants, ref.getLocus()) ) {
                        //for either true positives for false positives, this call isn't a "positive" at all if it's homRef in NA12878
                        if (!vcRaw.hasGenotypes() || !vcRaw.hasGenotype(sampleNameToCompare) || vcRaw.getGenotype(sampleNameToCompare).isHomRef() || vcRaw.getGenotype(sampleNameToCompare).isNoCall())
                            continue;
                        VariantContext vcTrimmed = GATKVariantContextUtils.trimAlleles(vcRaw.subContextFromSample(sampleNameToCompare), false, true);
                        Set<VariantContext> biallelics = new HashSet<>();
                        for (final VariantContext biallelic : GATKVariantContextUtils.splitVariantContextToBiallelics(vcTrimmed, true, GATKVariantContextUtils.GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL)) {
                            biallelics.add(biallelic);
                        }
                        for (final VariantContext vc : biallelics) {
                            // Do the alleles match between the input site and the KB site
                            if (cs.getVariantContext().hasSameAllelesAs(vc)) {
                                foundMatchingAllele = true;
                                if (csType.isTruePositive() && !cs.isMonomorphic() || (includeUnknowns && csType.isUnknown())) {
                                    data.add(new ROCDatum(true, cs.getVariantContext().isSNP(), (vc.hasAttribute("VQSLOD") ? vc.getAttributeAsDouble("VQSLOD", Double.NaN) : vc.getPhredScaledQual()), vc.getFiltersMaybeNull()));
                                    sitesWriter.notifyOfSite(AssessmentType.TRUE_POSITIVE, vc, cs);
                                } else if (csType.isFalsePositive() || (csType.isTruePositive() && cs.isMonomorphic())) {
                                    data.add(new ROCDatum(false, cs.getVariantContext().isSNP(), (vc.hasAttribute("VQSLOD") ? vc.getAttributeAsDouble("VQSLOD", Double.NaN) : vc.getPhredScaledQual()), vc.getFiltersMaybeNull()));
                                    sitesWriter.notifyOfSite(AssessmentType.FALSE_POSITIVE, vc, cs);
                                }
                                //unknown, discordant, and suspect consensus sites won't be counted in the ROC curve data (unknowns are counted if the appropriate arg is specified)
                            }
                            // Don't add this unmatched allele to the uncalled Sites because it is actually a false positive in hiConf mode since it exists in
                            // the input file but not in the KB
                            if (hiConf != null) {
                                foundMatchingAllele = true;
                            }
                        }
                    }
                    //TODO: uncalled sites may be a little off when there are het-non-ref genotypes (i.e. 1/2s)
                    if (!foundMatchingAllele) {
                        uncalledSites[cs.getVariantContext().isSNP() ? SNP_INDEX : INDEL_INDEX][cs.getType().isTruePositive() ? TP_INDEX : FP_INDEX]++;
                    }
                }
            }
        }

        // Find sites that are not in the KB and make them false positives if we are in a high confidence region
        if(hiConf!=null){
            for( final VariantContext vc : tracker.getValues(variants, ref.getLocus())){
                boolean foundMatchingAllele = false;
                for( final MongoVariantContext cs : kbSitesAtThisLoc ) {
                    //If this MongoVariantContext is not from the correct callset ignore it
                    if(!cs.getSupportingCallSets().contains(hiConf)){
                        continue;
                    }
                    if ((cs.getVariantContext().isSNP() || cs.getVariantContext().isIndel()) && (cs.getType().isTruePositive() || cs.getType().isFalsePositive())) {
                        if (!REQUIRE_REVIEWED || cs.isReviewed() ) {
                            if( cs.getVariantContext().hasSameAllelesAs(vc) ) {
                                foundMatchingAllele = true;
                            }
                        }
                    }
                }
                if(!foundMatchingAllele){
                    data.add( new ROCDatum( false, vc.isSNP(), (vc.hasAttribute("VQSLOD") ? vc.getAttributeAsDouble("VQSLOD", Double.NaN) : vc.getPhredScaledQual()), vc.getFiltersMaybeNull()));
                }
            }
        }


        return 1;
    }

    private void countUncalledBefore(final AlignmentContext context) {
        for( final MongoVariantContext cs : siteIterator.getSitesBefore(context.getLocation()) ) {
            //Check if we are in hiConf mode that this MongoVariantContext is from the correct callSet
            if(hiConf!=null && !cs.getSupportingCallSets().contains(hiConf)){
                continue;
            }

            if( (cs.getVariantContext().isSNP() || cs.getVariantContext().isIndel()) && (cs.getType().isTruePositive() || cs.getType().isFalsePositive()) ) {
                if( !REQUIRE_REVIEWED || cs.isReviewed() ) {
                    uncalledSites[cs.getVariantContext().isSNP() ? SNP_INDEX : INDEL_INDEX][cs.getType().isTruePositive() ? TP_INDEX : FP_INDEX]++;
                }
            }
        }
    }

    @Override
    public void onTraversalDone(final Integer result) {
        countUncalledBefore(lastLoc);
        super.onTraversalDone(result);
        final GATKReport report = calculateROCCurve(data, n_bin_snp, n_bin_indel, project, variants.getSource());
        try {
            PrintStream out_stream = new PrintStream(report_file);
            report.print(out_stream);
            out_stream.close();
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(report_file, e);
        }
        if (PLOTS_BASENAME != null) {
            plotROCCurves(report_file);
        }
    }

    /**
     * Create the GATK report which holds the ROC-like curve information for this data
     * @param data      the data to use
     * @param n_bin_snp   the number of partitions for the snps in the ROC curve
     * @param n_bin_indel   the number of partitions for the indels in the ROC curve
     * @param project   the project ID
     * @param vcf_name      the name of the input VCF
     * @return          the GATK report written to report_file
     */
    protected GATKReport calculateROCCurve(final List<ROCDatum> data, final int n_bin_snp, final int n_bin_indel, final String project, final String vcf_name) {
        final GATKReport report = GATKReport.newSimpleReportWithDescription("NA12878Assessment", "Evaluation of input variant callsets", "project", "vcf_name", "variation", "vqslod", "N_TP", "N_FP", "filter");
        Collections.sort(data); // sort by the LOD score
        final int[][][] rocData = new int[2][2][2]; //[SNP/INDEL][TP/FP][called/total]

        // loop over the sorted data and compute a ROC curve by looking at the fraction of called TPs and FPs
        for( final ROCDatum datum : data ) {
            rocData[datum.isSNP ? SNP_INDEX : INDEL_INDEX][datum.isTP ? TP_INDEX : FP_INDEX][TOTAL_INDEX]++;
        }
        report.addRow(project, vcf_name, "uncalled_SNPs", "NA", uncalledSites[SNP_INDEX][TP_INDEX], uncalledSites[SNP_INDEX][FP_INDEX],"NA");
        report.addRow(project, vcf_name, "uncalled_Indels", "NA", uncalledSites[INDEL_INDEX][TP_INDEX], uncalledSites[INDEL_INDEX][FP_INDEX],"NA");
        for( final boolean calcSNP : new boolean[]{true, false} ) {
            int numVariants = 0;
            final int varIndex = calcSNP? SNP_INDEX : INDEL_INDEX;
            final int total = rocData[varIndex][TP_INDEX][TOTAL_INDEX] + rocData[varIndex][FP_INDEX][TOTAL_INDEX];
            final int stepSize = Math.max(total / (calcSNP? n_bin_snp : n_bin_indel), 1);
            for( final ROCDatum datum : data ) {
                if( datum.isSNP != calcSNP ) { continue; }
                rocData[datum.isSNP ? SNP_INDEX : INDEL_INDEX][datum.isTP ? TP_INDEX : FP_INDEX][CALLED_INDEX]++;
                if( (numVariants+1) % stepSize == 0 ) {
                    report.addRow(project, vcf_name, calcSNP ? "SNPs" : "Indels", datum.lod,
                            (double)rocData[varIndex][TP_INDEX][CALLED_INDEX],
                            (double)rocData[varIndex][FP_INDEX][CALLED_INDEX],
                            datum.filterField);
                }
                numVariants++;
            }
        }

        logger.info("SNPs");
        logger.info("\tTotal # of true positives: " + rocData[SNP_INDEX][TP_INDEX][TOTAL_INDEX]);
        logger.info("\tTotal # of false positives: " + rocData[SNP_INDEX][FP_INDEX][TOTAL_INDEX]);
        logger.info("Indels");
        logger.info("\tTotal # of true positives: " + rocData[INDEL_INDEX][TP_INDEX][TOTAL_INDEX]);
        logger.info("\tTotal # of false positives: " + rocData[INDEL_INDEX][FP_INDEX][TOTAL_INDEX]);

        return report;
    }

    protected void plotROCCurves(final File report_file) {
        String RATIO_PLOT_FNAME = new String(PLOTS_BASENAME + "_ratio.png");
        String COUNT_PLOT_FNAME = new String(PLOTS_BASENAME + "_count.png");
        String[] fileNames = new String[2];
        if (!RScriptExecutor.RSCRIPT_EXISTS) {
            logger.info("RScript exists but executor not found.");
            logger.info("Plots will  not be generated.");
        } else {
            RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(new Resource(PLOT_RSCRIPT, ROCCurveNA12878.class));
            executor.addArgs(report_file.getAbsoluteFile(), COUNT_PLOT_FNAME, RATIO_PLOT_FNAME);
            logger.info("Executing: " + executor.getApproximateCommandLine());
            executor.exec();
        }
    }

    // private class to hold the data for use when computing ROC curves
    protected static class ROCDatum implements Comparable<ROCDatum> {
        public final boolean isTP;
        public final boolean isSNP;
        public final double lod;
        public final String filterField;

        public ROCDatum( final boolean isTP, final boolean isSNP, final double lod, final Set<String> filters ) {
            this.isTP = isTP;
            this.isSNP = isSNP;
            this.lod = lod;
            filterField = ( filters == null || filters.isEmpty() ? "PASS" : filters.iterator().next() );
        }

        @Override
        public int compareTo(final ROCDatum datum) {
            return Double.compare(datum.lod, this.lod);
        }

    }

}
