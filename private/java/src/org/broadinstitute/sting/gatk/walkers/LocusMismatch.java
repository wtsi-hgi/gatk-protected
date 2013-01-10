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

package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.Collection;

/**
 * Walker to calculate the number of mismatches, their base counts, and their quality sums at confidence ref sites" 
 */
@By(DataSource.REFERENCE)
public class LocusMismatch extends LocusWalker<String,Integer> implements TreeReducible<Integer> {
    @Output
    PrintStream out;

    //@Argument(fullName="confidentRefThreshold",doc="Set the lod score that defines confidence in ref, defaults to 4", required=false)
    //int confidentRefThreshold = 5;
    @Argument(fullName="maxNumMismatches",doc="Set the maximum number of mismatches at a locus before choosing not to use it in calculation. Defaults to 1.", required=false)
    int maxNumMismatches = 100;
    @Argument(fullName="minMappingQuality", doc ="Set the alignment quality below which to ignore reads; defaults to 30", required = false)
    int minMappingQuality = 1;
    @Argument(fullName="minDepth",doc="Set the minimum number of reads at a locus before choosing to use it in calculation. Defaults to 20.", required=false)
    int minDepth = 10;
    @Argument(fullName="maxDepth",doc="Set the minimum number of reads at a locus before choosing to use it in calculation. Defaults to 20.", required=false)
    int maxDepth = 100;
    @Argument(fullName="minBaseQuality", doc = "Set the base quality score below which to ignore bases in the pileup, defaults to 20", required = false)
    int minQualityScore = 1;
    @Argument(fullName="maxBaseQuality", doc = "Set the base quality score below which to ignore bases in the pileup, defaults to no restriction", required = false)
    int maxQualityScore = 99;
    @Argument(fullName="minMismatches", doc = "Minimum number of mismatches at a locus before a site is displayed", required = false)
    int minMismatches = 1;

    @Argument(fullName="skip", doc = "Only display every skip eligable sites.  Defaults to all sites", required = false)
    int skip = 1;

    private UnifiedGenotyperEngine ug;

    public void initialize() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES;
        ug = new UnifiedGenotyperEngine(getToolkit(), uac);

        // print the header
        out.printf("loc ref genotype genotypeQ depth nMM qSumMM A C G T%n");
    }

    public String map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        String result = null;

        ReadBackedPileup pileup = context.getBasePileup();
        if ( locusIsUsable(tracker, ref, pileup, context) ) {
            Genotype g = getGenotype(tracker, ref, context);
            if ( g != null )
                result = errorCounts( ref, pileup, g );
        }

        return result;
    }

    public Integer reduce( String map, Integer reduce  ) {
        if ( map != null && (reduce % skip == 0) )
            out.println(map);

        //if (reduce % skip == 0) System.out.printf("Keeping %d%n", reduce);

        return reduce + (map != null ? 1 : 0);
    }

    public Integer treeReduce( Integer reduce1, Integer reduce2 ) {
        return reduce1 + reduce2;
    }

    public Integer reduceInit() {
        return 1;
    }

    private String errorCounts( ReferenceContext ref, ReadBackedPileup pileup, Genotype g ) {
        int[] baseCounts = { 0, 0, 0, 0 };
        int usableDepth = 0;
        int nMismatches = 0;
        int qSumMismatches = 0;

        for ( PileupElement e : pileup ) {
            if ( useRead(e) ) {
                //System.out.printf("Using %s%n", e.getRead().getReadName());
                baseCounts[e.getBaseIndex()] += 1;
                usableDepth++;
                if ( ! BaseUtils.basesAreEqual(e.getBase(), ref.getBase()) ) {
                    nMismatches++;
                    qSumMismatches += e.getQual();
                }
            }
        }

        if ( nMismatches < maxNumMismatches && nMismatches >= minMismatches && usableDepth >= minDepth ) {
            StringBuffer baseCountString = new StringBuffer();
            for ( byte b : BaseUtils.BASES ) {
                baseCountString.append(baseCounts[BaseUtils.simpleBaseToBaseIndex(b)]);
                baseCountString.append(" ");
            }
            return String.format("%s %c %10s %5.2f %d %d %d %s",
                    pileup.getLocation(), ref.getBaseAsChar(),
                    getGenotypeClass(g), -10 * g.getLog10PError(),
                    usableDepth, nMismatches, qSumMismatches, baseCountString.toString());
        }

        return null;
    }

    private String getGenotypeClass(Genotype g) {
        if ( g.isHomRef() ) return "HOM-REF";
        else if ( g.isHet() ) return "HET";
        else if ( g.isHom() ) return "HOM-NONREF";
        else throw new ReviewedStingException("Unexpected genotype in getGenotypeClass " + g);
    }

    public boolean useRead( PileupElement e ) {
        if ( e.getRead().getMappingQuality() <= minMappingQuality ) {
            return false;
        } else if ( ! BaseUtils.isRegularBase( e.getBase() ) ) {
            return false;
        } else if ( e.getQual() <= minQualityScore || e.getQual() > maxQualityScore ) {
            return false;
        } else {
            return true;
        }
    }

    private boolean locusIsUsable( RefMetaDataTracker tracker, ReferenceContext ref, ReadBackedPileup pileup, AlignmentContext context ) {
        return BaseUtils.isRegularBase(ref.getBase()) &&
                pileup.getNumberOfElements() >= minDepth && pileup.getNumberOfElements() < maxDepth &&
                notCoveredByVariations(tracker, ref) &&
                pileupContainsNoNs(pileup);
//        pileupContainsNoNs(pileup) &&
//        baseIsConfidentRef(tracker,ref,context);
    }

    private boolean notCoveredByVariations( RefMetaDataTracker tracker, ReferenceContext ref ) {
        Collection<VariantContext> vcs = tracker.getValues(VariantContext.class);
        // TODO: check this logic. I think it's the best approximation of what was here before, but it's a different system
        if (vcs != null && vcs.size() > 0 ) {
                return false;
        }

        return true;
    }

    private boolean pileupContainsNoNs(ReadBackedPileup pileup) {
        for ( byte c : pileup.getBases() ) {
            if ( c == 'N' ) {
                return false;
            }
        }

        return true;
    }

    private Genotype getGenotype( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        VariantCallContext calls = ug.calculateLikelihoodsAndGenotypes(tracker,ref,context).get(0);
        if ( calls == null || calls.getNSamples() == 0 || !calls.isSNP() )
            return null;
        else {
            return calls.getGenotype(0);
        }
    }

//    private boolean baseIsConfidentRef( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
//        Pair<VariationCall, List<Genotype>> calls = ug.map(tracker,ref,context);
//        if ( calls == null || calls.first == null)
//            return false;
//        else {
//            VariationCall var = calls.getFirst();
//            return var.isReference() && var.getLog10PError() > confidentRefThreshold;
//            //return  ( var.isReference() > 0 && !calls.second.get(0).isVariant(ref.getBase()) && calls.second.get(0).getLog10PError() > confidentRefThreshold );
//        }
//    }

    public void onTraversalDone(Integer result) {
        ;
    }
}
