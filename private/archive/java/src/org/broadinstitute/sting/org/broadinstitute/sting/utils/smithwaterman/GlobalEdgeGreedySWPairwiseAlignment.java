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

package org.broadinstitute.sting.utils.smithwaterman;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;

import java.util.*;

/**
 * Pairwise discrete Smith-Waterman alignment with an edge greedy implementation
 *
 * ************************************************************************
 * ****                    IMPORTANT NOTE:                             ****
 * ****  This class assumes that all bytes come from UPPERCASED chars! ****
 * ************************************************************************
 *
 * User: ebanks
 */
public final class GlobalEdgeGreedySWPairwiseAlignment extends SWPairwiseAlignment {

    private final static boolean DEBUG_MODE = false;

    /**
     * Create a new greedy SW pairwise aligner
     *
     * @param reference the reference sequence we want to align
     * @param alternate the alternate sequence we want to align
     * @param parameters the SW parameters to use
     */
    public GlobalEdgeGreedySWPairwiseAlignment(final byte[] reference, final byte[] alternate, final Parameters parameters) {
        super(reference, alternate, parameters);
    }

    /**
     * Create a new SW pairwise aligner
     *
     * After creating the object the two sequences are aligned with an internal call to align(seq1, seq2)
     *
     * @param reference the reference sequence we want to align
     * @param alternate the alternate sequence we want to align
     * @param namedParameters the named parameter set to get our parameters from
     */
    public GlobalEdgeGreedySWPairwiseAlignment(final byte[] reference, final byte[] alternate, final SWParameterSet namedParameters) {
        this(reference, alternate, namedParameters.parameters);
    }

    /**
     * @see #GlobalEdgeGreedySWPairwiseAlignment(byte[], byte[], SWParameterSet) with original default parameters
     */
    public GlobalEdgeGreedySWPairwiseAlignment(byte[] reference, byte[] alternate) {
        this(reference, alternate, SWParameterSet.ORIGINAL_DEFAULT);
    }

    /**
     * Aligns the alternate sequence to the reference sequence
     *
     * @param reference  ref sequence
     * @param alternate  alt sequence
     */
    @Override
    protected void align(final byte[] reference, final byte[] alternate) {
        if ( reference == null || reference.length == 0 )
            throw new IllegalArgumentException("Non-null, non-empty reference sequences are required for the Smith-Waterman calculation");
        if ( alternate == null || alternate.length == 0 )
            throw new IllegalArgumentException("Non-null, non-empty alternate sequences are required for the Smith-Waterman calculation");

        final int forwardEdgeMatch = Utils.longestCommonPrefix(reference, alternate, Integer.MAX_VALUE);

        // edge case: one sequence is a strict prefix of the other
        if ( forwardEdgeMatch == reference.length || forwardEdgeMatch == alternate.length ) {
            alignmentResult = new SWPairwiseAlignmentResult(makeCigarForStrictPrefixAndSuffix(reference, alternate, forwardEdgeMatch, 0), 0);
            return;
        }

        int reverseEdgeMatch = Utils.longestCommonSuffix(reference, alternate, Integer.MAX_VALUE);

        // edge case: one sequence is a strict suffix of the other
        if ( reverseEdgeMatch == reference.length || reverseEdgeMatch == alternate.length ) {
            alignmentResult = new SWPairwiseAlignmentResult(makeCigarForStrictPrefixAndSuffix(reference, alternate, 0, reverseEdgeMatch), 0);
            return;
        }

        final int sizeOfRefToAlign = reference.length - forwardEdgeMatch - reverseEdgeMatch;
        final int sizeOfAltToAlign = alternate.length - forwardEdgeMatch - reverseEdgeMatch;

        // edge case: one sequence is a strict subset of the other accounting for both prefix and suffix
        final int minSizeToAlign = Math.min(sizeOfRefToAlign, sizeOfAltToAlign);
        if ( minSizeToAlign < 0 )
            reverseEdgeMatch += minSizeToAlign;
        if ( sizeOfRefToAlign <= 0 || sizeOfAltToAlign <= 0 ) {
            alignmentResult = new SWPairwiseAlignmentResult(makeCigarForStrictPrefixAndSuffix(reference, alternate, forwardEdgeMatch, reverseEdgeMatch), 0);
            return;
        }

        final byte[] refToAlign = Utils.trimArray(reference, forwardEdgeMatch, reverseEdgeMatch);
        final byte[] altToAlign = Utils.trimArray(alternate, forwardEdgeMatch, reverseEdgeMatch);

        final double[] sw = new double[(sizeOfRefToAlign+1)*(sizeOfAltToAlign+1)];
        if ( keepScoringMatrix ) SW = sw;
        final int[] btrack = new int[(sizeOfRefToAlign+1)*(sizeOfAltToAlign+1)];

        calculateMatrix(refToAlign, altToAlign, sw, btrack, OVERHANG_STRATEGY.INDEL);

        if ( DEBUG_MODE ) {
            System.out.println(new String(refToAlign) + " vs. " + new String(altToAlign));
            debugMatrix(sw, sizeOfRefToAlign+1, sizeOfAltToAlign+1);
            System.out.println("----");
            debugMatrix(btrack, sizeOfRefToAlign + 1, sizeOfAltToAlign + 1);
            System.out.println();
        }

        alignmentResult = calculateCigar(forwardEdgeMatch, reverseEdgeMatch, sizeOfRefToAlign, sizeOfAltToAlign, sw, btrack);
    }

    private void debugMatrix(final double[] matrix, final int dim1, final int dim2) {
        for ( int i = 0; i < dim1; i++ ) {
            for ( int j = 0; j < dim2; j++ )
                System.out.print(String.format("%.1f ", matrix[i * dim2 + j]));
            System.out.println();
        }
    }

    private void debugMatrix(final int[] matrix, final int dim1, final int dim2) {
        for ( int i = 0; i < dim1; i++ ) {
            for ( int j = 0; j < dim2; j++ )
                System.out.print(matrix[i*dim2 + j] + " ");
            System.out.println();
        }
    }

    /**
     * Creates a CIGAR for the case where the prefix/suffix match combination encompasses an entire sequence
     *
     * @param reference            the reference sequence
     * @param alternate            the alternate sequence
     * @param matchingPrefix       the prefix match size
     * @param matchingSuffix       the suffix match size
     * @return non-null CIGAR
     */
    private Cigar makeCigarForStrictPrefixAndSuffix(final byte[] reference, final byte[] alternate, final int matchingPrefix, final int matchingSuffix) {

        final List<CigarElement> result = new ArrayList<CigarElement>();

        // edge case: no D or I element
        if ( reference.length == alternate.length ) {
            result.add(makeElement(State.MATCH, matchingPrefix + matchingSuffix));
        } else {
            // add the first M element
            if ( matchingPrefix > 0 )
                result.add(makeElement(State.MATCH, matchingPrefix));

            // add the D or I element
            if ( alternate.length > reference.length )
                result.add(makeElement(State.INSERTION, alternate.length - reference.length));
            else // if ( reference.length > alternate.length )
                result.add(makeElement(State.DELETION, reference.length - alternate.length));

            // add the last M element
            if ( matchingSuffix > 0 )
                result.add(makeElement(State.MATCH, matchingSuffix));
        }

        return new Cigar(result);
    }

    /**
     * Calculates the CIGAR for the alignment from the back track matrix
     *
     * @param matchingPrefix       the prefix match size
     * @param matchingSuffix       the suffix match size
     * @param refLength            length of the reference sequence
     * @param altLength            length of the alternate sequence
     * @param sw                   the Smith-Waterman matrix to use
     * @param btrack               the back track matrix to use
     * @return non-null SWPairwiseAlignmentResult object
     */
    protected SWPairwiseAlignmentResult calculateCigar(final int matchingPrefix, final int matchingSuffix,
                                                       final int refLength, final int altLength,
                                                       final double[] sw, final int[] btrack) {

        final SWPairwiseAlignmentResult SW_result = calculateCigar(refLength, altLength, sw, btrack, OVERHANG_STRATEGY.INDEL);

        final LinkedList<CigarElement> lce = new LinkedList<CigarElement>(SW_result.cigar.getCigarElements());
        if ( matchingPrefix > 0 )
            lce.addFirst(makeElement(State.MATCH, matchingPrefix));
        if ( matchingSuffix > 0 )
            lce.addLast(makeElement(State.MATCH, matchingSuffix));

        return new SWPairwiseAlignmentResult(AlignmentUtils.consolidateCigar(new Cigar(lce)), 0);
    }
}