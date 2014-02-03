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

package org.broadinstitute.sting.gatk.walkers.na12878kb.assess;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.picard.filter.FilteringIterator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.filters.BadCigarFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.sting.utils.sam.GATKSamRecordFactory;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeType;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.File;
import java.util.*;

/**
 * Assess the quality of a single callset
 *
 * Keeps track of the Assessment values for SNPs and indels for a single
 * callset.  Doesn't get data at all, just manages the assessment of a single input
 * callset compared to the KB.
 *
 * The basic interface here is through a single function:
 *
 * A function that assesses the calls at a single locus
 *
 * @see #assessSite(java.util.List, java.util.List, boolean)
 *
 * That finds equivalent calls at a single site in a single callset compared to the calls
 * in the KB at that location (list two)
 *
 * @author depristo
 * @since 11/2012
 * @version 0.1
 */
public class Assessor {
    private final Logger logger = Logger.getLogger(Assessor.class);

    private final AssessNA12878.TypesToInclude typesToInclude;
    protected final String name;
    private final SAMFileReader bamReader;
    private final int minDepthForLowCoverage;
    private final Set<String> excludeKBSitesSupportedByOnlyTheseCallset;
    private final SitesWriter badSitesWriter;
    private final boolean ignoreFilters;
    private final int minPNonRef;
    private final int minGQ;

    Assessment SNPAssessment = new Assessment(AssessmentType.DETAILED_ASSESSMENTS);
    Assessment IndelAssessment = new Assessment(AssessmentType.DETAILED_ASSESSMENTS);


    /**
     * Create a new simple Assessor suitable for testing only
     *
     * @param name the name of our assessor
     */
    protected Assessor(final String name) {
        this(name, AssessNA12878.TypesToInclude.BOTH, Collections.<String>emptySet(), BadSitesWriter.NOOP_WRITER, null, 0, -1, 20, false);
    }

    /**
     * Create a new assess for a single callset
     *
     * @param name a non-null string name for callset
     * @param typesToInclude the types of variants to keep track off
     * @param excludeKBSitesSupportedByOnlyTheseCallset a non-null set of strings.  We will ignore calls in the KB
     *                                                  that are only supported by callsets in this set
     * @param badSitesWriter a bad sites writer that will be called for each assessed pair of call / KB sites
     *                       that may choose to write out interesting sites
     * @param bamReader an optional (can be null) BAM reader.  If not null we'll look at depth of coverage
     *                  for FNs using this reader
     * @param minDepthForLowCoverage if bamReader is provided, coverage at a FN will be considered "low coverage"
     *                               if its below this value
     * @param minPNonRef  if PLs against 0/0 are below this number, do not consider the site called; use -1 to ignore
     * @param minGQ  if GQ is below this number, do not consider the genotype called; use -1 to ignore
     * @param ignoreFilters if true, ignore the filter status of calls and use them all
     */
    public Assessor(final String name,
                    final AssessNA12878.TypesToInclude typesToInclude,
                    final Set<String> excludeKBSitesSupportedByOnlyTheseCallset,
                    final SitesWriter badSitesWriter,
                    final SAMFileReader bamReader,
                    final int minDepthForLowCoverage,
                    final int minPNonRef,
                    final int minGQ,
                    final boolean ignoreFilters) {
        if ( name == null ) throw new IllegalArgumentException("ROD name cannot be null");
        if ( name.equals("") ) throw new IllegalArgumentException("ROD name cannot be the empty string");
        if ( typesToInclude == null ) throw new IllegalArgumentException("typesToInclude cannot be null");
        if ( excludeKBSitesSupportedByOnlyTheseCallset == null ) throw new IllegalArgumentException("excludeKBSitesSupportedByOnlyTheseCallset cannot be null");
        if ( badSitesWriter == null ) throw new IllegalArgumentException("badSitesWriter cannot be null");

        this.name = name;
        this.typesToInclude = typesToInclude;
        this.bamReader = bamReader;
        this.minDepthForLowCoverage = minDepthForLowCoverage;
        this.minPNonRef = minPNonRef;
        this.minGQ = minGQ;
        this.excludeKBSitesSupportedByOnlyTheseCallset = excludeKBSitesSupportedByOnlyTheseCallset;
        this.badSitesWriter = badSitesWriter;
        this.ignoreFilters = ignoreFilters;
    }

    public Assessment getSNPAssessment() { return SNPAssessment; }
    public Assessment getIndelAssessment() { return IndelAssessment; }

    /**
     * Should we include the variant vc in our assessment?
     *
     * @param vc a VariantContext to potentially include
     * @return true if VC should be included, false otherwise
     */
    @Requires("vc != null")
    private boolean includeVariant(final VariantContext vc) {
        switch ( typesToInclude ) {
            case BOTH: return true;
            case SNPS: return vc.isSNP();
            case INDELS: return ! vc.isSNP();
            default:
                throw new IllegalStateException("Unexpected enum " + typesToInclude);
        }
    }

    /**
     * @see #assessSite(java.util.List, java.util.List, boolean, boolean) with okayToMiss=false
     */
    public void assessSite(final List<VariantContext> vcs, final List<MongoVariantContext> consensusSites, final boolean onlyReviewed) {
        assessSite(vcs, consensusSites, onlyReviewed, false);
    }

    /**
     * Access a single locus that contains vcs calls from a single callset and consensusSites calls from the KB
     *
     * Finds equivalent pairs of calls between vcs and consensusSites, and updates the TP/FN/etc assessment
     * data based on these pairings.
     *
     * vcs can contains calls not solely based on NA12878, as the input VCs will be subsetted down to NA12878
     * on the fly
     *
     * @param vcs a non-null list of calls in a single callset.  If empty, means that no calls were present
     *            in the callset at consensusSites
     * @param consensusSites a non-null list of calls in the KB.  If empty, vcs are interpreted as having
     *                       no potential equivalents in the KB
     * @param onlyReviewed if true, only consider reviewed sites during the assessment
     * @param okayToMiss   if true, we will not penalize for any FALSE_NEGATIVES
     */
    public void assessSite(final List<VariantContext> vcs, final List<MongoVariantContext> consensusSites, final boolean onlyReviewed, final boolean okayToMiss) {
        if ( vcs == null ) throw new IllegalArgumentException("vcs cannot be null");
        if ( consensusSites == null ) throw new IllegalArgumentException("consensusSites cannot be null");

        // TODO -- check that these lists are in some sense meaningful to compare (locations, for instance)

        if ( vcs.isEmpty() ) {
            // missed consensus site(s)
            for ( final MongoVariantContext site : consensusSites ) {
                if ( logger.isDebugEnabled() ) logger.debug("Missed site in " + name + " site = " + site);
                assessMatchedCallWithKB(null, site, onlyReviewed, okayToMiss, true);
            }
        } else {
            final Set<VariantContext> biallelics = new HashSet<VariantContext>();

            for( final VariantContext vcRaw : vcs ) {
                if ( vcRaw.getAlternateAlleles().isEmpty() ) // skip sites without alt allele
                    continue;

                // allow sites only VCs to be evaluated as though they are just NA12878 calls
                final VariantContext na12878vc = subsetToNA12878(vcRaw);
                if ( na12878vc != null && na12878vc.isVariant() ) {
                    for ( final VariantContext biallelic : GATKVariantContextUtils.splitVariantContextToBiallelics(na12878vc, true, GATKVariantContextUtils.GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL)) {
                        if ( biallelic != null && biallelic.getStart() > na12878vc.getStart() ) {
                            // trimming the variant moved the variant forward on the genome, which can happen but
                            // means that the input VCF had a very strange multi-allelic structure
                            logger.warn("Biallelic split in " + name + " moved a variant into the future " + biallelic + " from " + na12878vc);
                        } else
                            biallelics.add(biallelic);
                    }
                }
            }

            for ( final Pair<VariantContext, MongoVariantContext> match : matchCallsWithKB(biallelics, consensusSites) ) {
                final VariantContext biallelic = match.getFirst();
                final MongoVariantContext consensusSite = match.getSecond();
                // because we split VariantContexts into component parts it is possible that records that were originally
                // of MIXED type can look 0/1 instead of 1/2 (in VCF parlance), so we don't want to assess genotype concordance
                // in such cases (hence the check for biallelics having at most 1 element).
                assessMatchedCallWithKB(biallelic, consensusSite, onlyReviewed, okayToMiss, biallelics.size() <= 1);
            }
        }
    }

    private VariantContext subsetToNA12878(final VariantContext vcRaw) {
        if ( vcRaw.hasGenotypes() ) {
            final VariantContext vcNA12878 = GATKVariantContextUtils.trimAlleles(vcRaw.subContextFromSample("NA12878"), false, true);
            if (!vcNA12878.hasGenotype("NA12878"))
                throw new UserException.BadInput("The input file does not contain NA12878. VCFs with genotypes must have NA12878 as one of the samples present, otherwise you must use a sites-only VCF. If your VCF is from NA12878 make sure that the sample name is actually 'NA12878'.");
            final Genotype na12878 = vcNA12878.getGenotype("NA12878");
            if (na12878.hasPL() && na12878.getPL()[0] < minPNonRef)
                return null;
            else {
                return vcNA12878;
            }
        } else {
            return vcRaw;
        }
    }

    protected List<VariantContext> biallelizeVarantContextList(final List<VariantContext> vcs) {
        final List<VariantContext> biallelics = new LinkedList<VariantContext>();

        for( final VariantContext vcRaw : vcs ) {
            if ( vcRaw.getAlternateAlleles().isEmpty() ) // skip sites without alt allele
                continue;

            // allow sites only VCs to be evaluated as though they are just NA12878 calls
            final VariantContext na12878vc = vcRaw.hasGenotypes() ? GATKVariantContextUtils.trimAlleles(vcRaw.subContextFromSample("NA12878"), false, true) : vcRaw;
            if ( na12878vc.isVariant() ) {
                for ( final VariantContext biallelic : GATKVariantContextUtils.splitVariantContextToBiallelics(na12878vc, true, GATKVariantContextUtils.GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL)) {
                    if ( biallelic != null && biallelic.getStart() > na12878vc.getStart() ) {
                        // trimming the variant moved the variant forward on the genome, which can happen but
                        // means that the input VCF had a very strange multi-allelic structure
                        logger.warn("Biallelic split in " + name + " moved a variant into the future " + biallelic + " from " + na12878vc);
                    } else
                        biallelics.add(biallelic);
                }
            }
        }
        return biallelics;
    }

    /**
     * Match called variants with consensus sites at this location
     *
     * Walks over the consensusSites, looking for a call in calls that matches it.  Consensus sites that match
     * are removed, so it's not possible for multiple equal calls (which shouldn't be passed in) to have matches
     * in the db.  Calls that don't have matches are returned as (call, null) and consensus sites that don't
     * have matches are returned as (null, consensusSite).  The returned list ensures that all elements of both
     * calls and consensusSites occur exactly once in the collection of pairs
     *
     * Uses MongoVariantContext.match to determine if a consensus site matches.
     *
     * @param calls a non-null collection of calls at this location
     * @param consensusSites a non-null collection of consensus sites at this location
     * @return a non-null collection of pairs of matched calls / consensus sites
     */
    @Requires({"calls != null", "consensusSites != null"})
    @Ensures("result != null")
    protected Collection<Pair<VariantContext, MongoVariantContext>> matchCallsWithKB(final Collection<VariantContext> calls,
                                                                                     final Collection<MongoVariantContext> consensusSites) {
        final Set<MongoVariantContext> unmatchedConsensusSites = new HashSet<MongoVariantContext>(consensusSites);
        final List<Pair<VariantContext, MongoVariantContext>> matches = new LinkedList<Pair<VariantContext, MongoVariantContext>>();

        // find matches for every call
        for ( final VariantContext call : calls ) {
            final MongoVariantContext match = findMatching(call, unmatchedConsensusSites);
            if ( match != null ) unmatchedConsensusSites.remove(match);
            matches.add(new Pair<VariantContext, MongoVariantContext>(call, match));
        }

        // add in the unmatched consensus sites
        for ( final MongoVariantContext unmatched : unmatchedConsensusSites ) {
            matches.add(new Pair<VariantContext, MongoVariantContext>(null, unmatched));
        }

        return matches;
    }

    /**
     * Find a matching consensus site for vc in the collection of sites
     *
     * @param vc a non-null VariantContext we want to match
     * @param consensusSites a collection of potential consensus sites to find our match in
     * @return a matching site in consensusSites, or null if none could be found
     */
    @Requires({"vc != null", "consensusSites != null"})
    protected MongoVariantContext findMatching(final VariantContext vc, final Collection<MongoVariantContext> consensusSites ) {
        for ( final MongoVariantContext site : consensusSites )
            if ( site.matches(vc) )
                return site;
        return null;
    }

    /**
     * Assess a call with its matched consensus site
     *
     * If call and consensusSite are both not null, then we interpret this as a equivalent call, and the
     * assessment type is computed as such.  If call is null, consensusSite must not be null and it's interpreted
     * as a KB site that wasn't present in the call set.  If call is not null and consensusSite is null, it's
     * intrepreted as a call without an equivalent record in the KB.  It's an error if both are null.
     *
     * @param call a potentially null VariantContext call
     * @param consensusSite a potentially null consensus Site
     * @param okayToMiss   if true, we will not penalize for any FALSE_NEGATIVES
     * @param okayToAssessGenotypes if false, we will not penalize for possible genotype discordance
     */
    protected void assessMatchedCallWithKB(final VariantContext call, final MongoVariantContext consensusSite,
                                           final boolean onlyReviewed, final boolean okayToMiss, final boolean okayToAssessGenotypes) {
        if ( call == null && consensusSite == null ) throw new IllegalArgumentException("both call and consensusSite cannot be null");
        if ( call != null && consensusSite != null && ( call.getStart() != consensusSite.getStart() ) )
            throw new IllegalArgumentException("Call and consensusSite don't start at the same position! " + call + " consensus " + consensusSite);

        final VariantContext vc = call != null ? call : consensusSite.getVariantContext();

        if ( ! includeVariant(vc) )
            return;

        if ( onlyReviewed && (consensusSite == null || ! consensusSite.isReviewed()) )
            return;

        final AssessmentType type = figureOutAssessmentType(call, consensusSite, okayToMiss);
        final Assessment assessment = vc.isSNP() ? SNPAssessment : IndelAssessment;
        assessment.inc(type);

        // determine if we have called the site correctly but failed to genotype it properly
        boolean genotypeDiscordance = false;
        if ( okayToAssessGenotypes && call != null && consensusSite != null && type == AssessmentType.TRUE_POSITIVE &&
                call.hasGenotype("NA12878") && ! isNotUsableCall(call) && consensusSite.isPolymorphic()) {
            final List<Allele> alleles = vc.getAlleles();
            final Genotype consensusGT = consensusSite.getGt().toGenotype(alleles);
            if ( consensusGT.getType() == GenotypeType.HET || consensusGT.getType() == GenotypeType.HOM_VAR ) {
                final Genotype callGT = call.getGenotype("NA12878");
                if ( !callGT.hasGQ() || callGT.getGQ() >= minGQ ) {
                    final boolean concordant = consensusGT.getType() == callGT.getType();
                    genotypeDiscordance = ! concordant;
                    assessment.incGenotypingAccuracy(consensusGT.getType(), concordant);
                }
            }
        }

        if ( logger.isDebugEnabled() ) logger.debug("Assessing site " + name + " call " + call + " against consensus " + consensusSite);

        badSitesWriter.notifyOfSite(genotypeDiscordance ? AssessmentType.GENOTYPE_DISCORDANCE : type, vc, consensusSite);
    }

    /**
     * Look at two equivalent calls (call from the input callset and consensusSite from the KB) and decide
     * the AssessmentType (TP, FN, etc) appropriate for their pairing
     *
     * Note that call and consensus cannot both be null at the same time, but either can be null individually
     * indicating that one call was present and the other was missing.
     *
     * @param call a single VariantContext representing a call in a input callset potentially equivalent to consensus site
     *             If null, then no equivalent call was found in the input callset for consensusSite
     * @param consensusSite a single MongoVariantContext representing a KB consensus call equivalent to call.  If
     *                      null, indicates that there was no equivalent call in the KB for call.
     * @param okayToMiss   if true, we will not penalize for any FALSE_NEGATIVES
     * @return a non-null AssessmentType
     */
    @Ensures("result != null")
    protected AssessmentType figureOutAssessmentType(final VariantContext call, final MongoVariantContext consensusSite, final boolean okayToMiss) {
        if ( call == null && consensusSite == null ) throw new IllegalArgumentException("Both call and consensusSite cannot be null");

        final boolean consensusTP = consensusSite != null && !excludeConsensusSite(consensusSite) && consensusSite.getType().isTruePositive()
                && !consensusSite.getPolymorphicStatus().isMonomorphic(); // discordant genotypes should be allowed through (since it means at least one call must be polymorphic)
        final boolean consensusFP = consensusSite != null && !excludeConsensusSite(consensusSite) && consensusSite.getType().isFalsePositive();

        if ( call != null ) {
            if ( consensusTP ) {
                if ( isNotUsableCall(call) )
                    return AssessmentType.FALSE_NEGATIVE_CALLED_BUT_FILTERED;
                else if ( likelyWouldBeFiltered(call) )
                    // note that we don't consider how we might potentially filter a site that at a TP
                    return AssessmentType.TRUE_POSITIVE;
                    //return AssessmentType.FALSE_NEGATIVE_CALLED_BUT_WOULD_BE_FILTERED;
                else
                    return AssessmentType.TRUE_POSITIVE;
            } else if (consensusSite != null && consensusSite.getType().isTruePositive() && consensusSite.getPolymorphicStatus().isMonomorphic() ) {
                return AssessmentType.FALSE_POSITIVE_MONO_IN_NA12878;
            } else if ( consensusFP ) {
                if ( isNotUsableCall(call) )
                    return AssessmentType.CORRECTLY_FILTERED;
                else if ( likelyWouldBeFiltered(call) )
                    return AssessmentType.REASONABLE_FILTERS_WOULD_FILTER_FP_SITE;
                else
                    return AssessmentType.FALSE_POSITIVE_SITE_IS_FP;
            } else if ( consensusSite != null && consensusSite.getType().isUnknown() ) {
                return AssessmentType.CALLED_IN_DB_UNKNOWN_STATUS;
            } else if ( consensusSite == null && ! isNotUsableCall(call) ) {
                return AssessmentType.CALLED_NOT_IN_DB_AT_ALL;
            } else {
                return AssessmentType.NOT_RELEVANT;
            }
        } else if ( consensusTP ) { // call == null
            // if it's a complex event, just ignore it (because we may have called it with a different representation in the VCF)
            if ( okayToMiss || consensusSite.isComplexEvent() )
                return AssessmentType.NOT_RELEVANT;

            if ( bamReader != null ) {
                return sufficientDepthToCall(consensusSite)
                        ? AssessmentType.FALSE_NEGATIVE_NOT_CALLED_AT_ALL
                        : AssessmentType.FALSE_NEGATIVE_NOT_CALLED_BUT_LOW_COVERAGE;
            } else {
                return AssessmentType.FALSE_NEGATIVE_NOT_CALLED_AT_ALL;
            }
        } else if ( consensusFP ) { // call == null
            return AssessmentType.CORRECTLY_UNCALLED;
        } else {
            return AssessmentType.NOT_RELEVANT;
        }
    }

    private boolean isNotUsableCall(final VariantContext vc) {
        return ! ignoreFilters && vc.isFiltered();
    }

    /**
     * Should consensusSite be excluded from the analysis?
     *
     * @param consensusSite the consensus site we may want to exclude
     * @return true if we shouldn't consider consensusSite during the analysis
     */
    @Requires("consensusSite != null")
    private boolean excludeConsensusSite(final MongoVariantContext consensusSite) {
        return excludeKBSitesSupportedByOnlyTheseCallset.containsAll(consensusSite.getSupportingCallSets());
    }

    /**
     * Create a SAMFileReader appropriately initialized for getting the DoC of sites
     *
     * @param bam a file pointing to a BAM file
     * @return a SAMFileReader that reads reads from bam
     */
    public static SAMFileReader makeSAMFileReaderForDoCInBAM(final File bam) {
        final SAMFileReader bamReader = new SAMFileReader(bam);
        bamReader.setSAMRecordFactory(new GATKSamRecordFactory());
        bamReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        return bamReader;
    }

    /**
     * If a BAM is provided, assess whether there's sufficient coverage to call the site
     *
     * @param falseNegative a site that was missed in the call set
     * @return true if there's enough coverage at the site in the BAM to likely make a call
     */
    @Requires({"bamReader != null", "falseNegative != null"})
    private boolean sufficientDepthToCall(final MongoVariantContext falseNegative) {
        final int depth = getDepthAtLocus(falseNegative.getChr(), falseNegative.getStart());
        return depth >= minDepthForLowCoverage;
    }

    /**
     * Get the depth of coverage at position on chr
     *
     * @param chr the chr we want to query
     * @param position the position we want to query
     * @return the depth at chr:position
     */
    protected int getDepthAtLocus(final String chr, final int position) {
        // set up the query and wrapping filtering iterators
        CloseableIterator<SAMRecord> it = bamReader.queryOverlapping(chr, position, position);
        it = new FilteringIterator(it, new BadCigarFilter());
        it = new FilteringIterator(it, new DuplicateReadFilter());

        final LocusIteratorByState libs = new LocusIteratorByState(bamReader, it);
        final AlignmentContext context = libs.advanceToLocus(position, false);
        int depth = 0;
        if ( context != null ) {
            // need to remove low quality reads/bases
            depth = context.getBasePileup().getBaseAndMappingFilteredPileup(20, 20).depthOfCoverage();
        }
        it.close();
        return depth;
    }

    /**
     * Transform the assessments to simpler forms
     */
    public void simplifyAssessments() {
        SNPAssessment = SNPAssessment.simplify();
        IndelAssessment = IndelAssessment.simplify();
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
    @Requires("vc != null")
    private boolean likelyWouldBeFiltered(final VariantContext vc) {
        final double FS = vc.getAttributeAsDouble("FS", 0.0);
        final double QD = vc.getAttributeAsDouble("QD", 20.0);
        final double MQ = vc.getAttributeAsDouble("MQ", 50.0);

        if ( vc.isSNP() ) {
            return FS > 40 || QD < 2 || MQ < 40;
        } else {
            return FS > 200 || QD < 2;
        }
    }
}