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

package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;

import java.util.*;

/**
 * Create a consensus MongoVariantContext from one or more supporting call sets
 *
 * User: depristo
 * Date: 2/18/13
 * Time: 2:55 PM
 */
public class ConsensusMaker {
    /**
     * Key interface: create a consensus based on joining the evidence from allSupportingCalls
     *
     * Looks at the supporting calls, and creates a consensus that reflects the TP/FP etc status
     * of each and the genotypes.  Weighs the result heavily (potentially entirely) on supporting
     * calls tagged as reviews.  The resulting consensus is at the same site as allSupportingCalls
     * and same alleles, and the supporting call set field will reflect the information in allSupportingCalls.
     *
     * Supports having multiple equivalent supporting calls (two depristo records, for instance); it
     * will use the most recent result.
     *
     * @param allSupportingCalls a non-null, non-empty collection of supporting calls
     * @return a non-null MongoVariantContext consensus
     */
    public MongoVariantContext makeConsensus(final Collection<MongoVariantContext> allSupportingCalls) {
        if ( allSupportingCalls == null || allSupportingCalls.isEmpty())
            throw new IllegalArgumentException("allSupportingCalls must be non-null, not empty collection");

        // make sure the input MVCs can be safely combined into a consensus
        final MongoVariantContext firstMVC = allSupportingCalls.iterator().next();
        for ( final MongoVariantContext mvc : allSupportingCalls )
            ensureSafeToCombineInConsensus(firstMVC, mvc);

        // only the most recent call from each call set is used
        final Collection<MongoVariantContext> callsForConsensus = mostRecentCalls(allSupportingCalls);

        final VariantContextBuilder builder = new VariantContextBuilder();
        final VariantContext first = firstMVC.getVariantContext();
        builder.chr(first.getChr()).start(first.getStart()).stop(first.getEnd());

        final List<Allele> alleles = first.getAlleles();
        builder.alleles(alleles);

        // iteration is over allSupportingCalls so we get the names all correct
        final LinkedHashSet<String> supportingCallSets = new LinkedHashSet<>();
        for ( final MongoVariantContext vc : allSupportingCalls ) {
            supportingCallSets.addAll(vc.getSupportingCallSets());
        }

        final TruthStatus type = determineBayesianTruthEstimate(callsForConsensus);
        final PolymorphicStatus status = determinePolymorphicStatus(callsForConsensus);
        final Genotype gt = consensusGT(type, status, new LinkedList<>(alleles), callsForConsensus);
        final boolean isReviewed = isReviewed(callsForConsensus);
        final boolean isComplexEvent = isComplexEvent(callsForConsensus);

        return MongoVariantContext.create(new LinkedList<>(supportingCallSets), builder.make(), type, new Date(), gt, isReviewed, isComplexEvent);
    }

    /**
     * Make sure the canon and test can be combined in a consensus
     * @param canon the first to test
     * @param test the second to test
     * @throws IllegalArgumentException if canon and test aren't safe to merge
     */
    private void ensureSafeToCombineInConsensus(final MongoVariantContext canon, final MongoVariantContext test) {
        if ( ! canon.getChr().equals(test.getChr()) )
            throw new IllegalArgumentException("Tried to make consensus from non-equivalent supporting calls (not equal chromosomes): " + canon + " vs " + test);

        if ( canon.getStart() != test.getStart() )
            throw new IllegalArgumentException("Tried to make consensus from non-equivalent supporting calls (not equal start): " + canon + " vs " + test);

        if ( canon.getStop() != test.getStop() )
            throw new IllegalArgumentException("Tried to make consensus from non-equivalent supporting calls (not equal end): " + canon + " vs " + test);

        if ( ! canon.getRefAllele().equals(test.getRefAllele()) || ! canon.getAltAllele().equals(test.getAltAllele()))
            throw new IllegalArgumentException("Tried to make consensus from non-equivalent supporting calls (not equal alleles): " + canon + " vs " + test);
    }

    /**
     * Is at least one call in individualCalls a reviewed call?
     * @param individualCalls a collection of calls to consider
     * @return true if at least one call in individualCalls is a review, false otherwise
     */
    @Requires("individualCalls != null")
    protected boolean isReviewed(final Collection<MongoVariantContext> individualCalls) {
        for ( final MongoVariantContext vc : individualCalls )
            if ( vc.isReviewed() )
                return true;
        return false;
    }

    /**
     * Is at least one call in individualCalls marked as complex?
     * @param individualCalls a collection of calls to consider
     * @return true if at least one call in individualCalls is marked as complex, false otherwise
     */
    @Requires("individualCalls != null")
    protected boolean isComplexEvent(final Collection<MongoVariantContext> individualCalls) {
        for ( final MongoVariantContext vc : individualCalls )
            if ( vc.isComplexEvent() )
                return true;
        return false;
    }

    /**
     * Get the calls from individualCalls that should be used to create the consensus
     *
     * Potentially selects a subset of the calls to actually build the consensus.  The subset
     * can be based on only those with reviews, only the most recent from each call sets, or
     * other criteria that increase the quality of the consensus
     *
     * @param individualCalls a non-empty collection of calls
     * @return a non-empty subset of individualCalls
     */
    @Requires("individualCalls != null && ! individualCalls.isEmpty()")
    @Ensures("! result.isEmpty()")
    protected Collection<MongoVariantContext> selectCallsForConsensus(final Collection<MongoVariantContext> individualCalls) {
        final List<MongoVariantContext> reviewed = new LinkedList<>();

        for ( final MongoVariantContext vc : mostRecentCalls(individualCalls) )
            if ( vc.isReviewed() ) reviewed.add(vc);

        return reviewed.isEmpty() ? individualCalls : reviewed;
    }

    /**
     * Get the most recent call for each call set in individualCalls
     *
     * @param individualCalls a collection of calls
     * @return a subset of individualCalls
     */
    @Requires("individualCalls != null")
    @Ensures("result != null")
    protected List<MongoVariantContext> mostRecentCalls(final Collection<MongoVariantContext> individualCalls) {
        final Map<String, List<MongoVariantContext>> byCallSet = new LinkedHashMap<>();

        for ( final MongoVariantContext vc : individualCalls ) {
            if ( vc.getSupportingCallSets().size() != 1 )
                throw new IllegalArgumentException("Expected exactly one supporting call set but got " + vc.getSupportingCallSets());

            List<MongoVariantContext> calls = byCallSet.get(vc.getCallSetName());
            if ( calls == null ) {
                calls = new ArrayList<>();
                byCallSet.put(vc.getCallSetName(), calls);
            }
            calls.add(vc);
        }

        final List<MongoVariantContext> uniques = new LinkedList<>();
        for ( final List<MongoVariantContext> callsFor1Callset : byCallSet.values() ) {
            Collections.sort(callsFor1Callset, new Comparator<MongoVariantContext>() {
                @Override
                public int compare(MongoVariantContext o1, MongoVariantContext o2) {
                    return -1 * o1.getDate().compareTo(o2.getDate());
                }
            });
            uniques.add(callsFor1Callset.get(0));
        }

        return uniques;
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

    /**
     * This is the previous version of the method to determine truth status, kept around for reference.
     * Uses only reviewed sites if they are present; otherwise uses all sites.  Calls into TruthStatus.makeConsensus()
     * to get a consensus truth type.
     *
     * @param individualCalls   the calls at the site
     * @return non-null truth status
     */
    @Deprecated
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

    private static final double CONFIDENCE_THRESHOLD = 0.8;

    /**
     * Determines the truth status of the site given the calls using a Bayesian estimate
     *
     * @param calls   the calls at the site
     * @return non-null truth status
     */
    protected TruthStatus determineBayesianTruthEstimate(final Collection<MongoVariantContext> calls) {

        final double LofTP = calculateLikelihood(calls, TruthStatus.TRUE_POSITIVE);
        final double LofFP = calculateLikelihood(calls, TruthStatus.FALSE_POSITIVE);
        final double LofSuspect = calculateLikelihood(calls, TruthStatus.SUSPECT);
        final double totalL = LofTP + LofFP + LofSuspect;

        final TruthStatus status;
        if ( LofTP / totalL >= CONFIDENCE_THRESHOLD )
            status = TruthStatus.TRUE_POSITIVE;
        else if ( LofFP / totalL >= CONFIDENCE_THRESHOLD )
            status = TruthStatus.FALSE_POSITIVE;
        else if ( LofSuspect / totalL >= CONFIDENCE_THRESHOLD )
            status = TruthStatus.SUSPECT;
        else if ( allStatusesAreUnknown(calls) )
            status = TruthStatus.UNKNOWN;
        else
            status = TruthStatus.DISCORDANT;

        return status;
    }

    /**
     * Determines the likelihood of the given truth status
     *
     * @param calls   the calls at the site
     * @param status  the truth status to assess
     * @return likelihood
     */
    private double calculateLikelihood(final Collection<MongoVariantContext> calls, final TruthStatus status) {

        // for now let's just multiply through the confidences

        double confidence = 1.0;
        for ( final MongoVariantContext call : calls ) {
            if ( call.getType() != TruthStatus.UNKNOWN )
                confidence *= (call.getType() == status ? call.getConfidence() : 1.0 - call.getConfidence());
        }
        return confidence;
    }

    /**
     * Determines whether are calls are UNKNOWN status
     *
     * @param calls   the calls to assess
     * @return true if all calls have UNKNOWN truth status
     */
    private boolean allStatusesAreUnknown(final Collection<MongoVariantContext> calls) {
        for ( final MongoVariantContext call : calls ) {
            if ( call.getType() != TruthStatus.UNKNOWN )
                return false;
        }
        return true;
    }
}
